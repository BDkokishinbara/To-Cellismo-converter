"""
シングルセル解析ファイル変換ウェブアプリケーション

CSV / h5ad / RDS / MEX のいずれかを h5mu 形式に変換する Flask アプリ。
変換と並行して UMAP 可視化をバックグラウンドジョブで実行する機能も持つ。
"""
from flask import Flask, render_template, request, send_file, jsonify
import os
from pathlib import Path
from werkzeug.utils import secure_filename
import uuid
from datetime import datetime
from dotenv import load_dotenv

from converters.csv_converter import csv_to_h5mu, validate_csv
from converters.rds_converter import rds_to_h5mu, check_rpy2_available
from converters.mex_converter import mex_to_h5mu, validate_mex_directory
from converters.h5ad_converter import h5ad_to_h5mu
from converters.umap_visualizer import umap_from_h5mu
import shutil
import zipfile
import tarfile
import traceback
import threading
import time
from collections import OrderedDict

# .env ファイルから環境変数を読み込む（あれば）
load_dotenv()

app = Flask(__name__)

# 環境変数から各種設定を取得（.env または OS の環境変数で上書き可能）
app.config['SECRET_KEY'] = os.getenv('SECRET_KEY', 'dev-secret-key-change-in-production')
app.config['UPLOAD_FOLDER'] = os.getenv('UPLOAD_FOLDER', 'uploads')
app.config['OUTPUT_FOLDER'] = os.getenv('OUTPUT_FOLDER', 'outputs')
# シングルセル h5ad/h5mu は GB 級が一般的なので 5 GB を既定にする。
# 必要なら MAX_CONTENT_LENGTH 環境変数（バイト単位）で上書き可能。
app.config['MAX_CONTENT_LENGTH'] = int(os.getenv('MAX_CONTENT_LENGTH', 5 * 1024 * 1024 * 1024))

# アップロード可能な拡張子（.tar.gz は別途処理）
ALLOWED_EXTENSIONS = {'csv', 'rds', 'mtx', 'tsv', 'gz', 'zip', 'tar', 'h5ad'}


def allowed_file(filename):
    """指定されたファイル名の拡張子が許可されているか判定する。"""
    # .tar.gz は 2 段拡張子なので個別に判定
    if filename.lower().endswith('.tar.gz'):
        return True
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS


def generate_unique_filename(prefix='file'):
    """タイムスタンプ + UUID で衝突しにくいファイル名を生成する。"""
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    unique_id = str(uuid.uuid4())[:8]
    return f"{prefix}_{timestamp}_{unique_id}"


def find_mex_directory(root_dir):
    """MEX 形式の必須ファイル群（matrix.mtx, features/genes.tsv, barcodes.tsv）
    が含まれるディレクトリを探す。

    Parameters
    ----------
    root_dir : str
        検索の起点ディレクトリ。

    Returns
    -------
    str or None
        見つかれば該当ディレクトリのパス。無ければ None。
    """
    from converters.mex_converter import find_file

    root_path = Path(root_dir)

    # まずは起点ディレクトリ直下をチェック
    if find_file(root_path, ['matrix.mtx', 'matrix.mtx.gz']):
        return str(root_path)

    # 見つからなければサブディレクトリを再帰的に探索
    for dirpath, dirnames, filenames in os.walk(root_dir):
        dir_path = Path(dirpath)
        if find_file(dir_path, ['matrix.mtx', 'matrix.mtx.gz']):
            return str(dir_path)

    return None


def _save_upload(file_storage):
    """アップロードされた FileStorage を UPLOAD_FOLDER に一意な名前で保存する。

    Returns
    -------
    dict
        original_name / upload_path / file_extension を含む情報辞書。
    """
    filename = secure_filename(file_storage.filename)
    unique_name = generate_unique_filename('upload')
    if filename.lower().endswith('.tar.gz'):
        ext_for_path = 'tar.gz'
        file_extension = 'gz'
    elif '.' in filename:
        file_extension = filename.rsplit('.', 1)[1].lower()
        ext_for_path = file_extension
    else:
        file_extension = ''
        ext_for_path = ''

    upload_path = os.path.join(app.config['UPLOAD_FOLDER'], f"{unique_name}.{ext_for_path}" if ext_for_path else unique_name)
    file_storage.save(upload_path)
    return {
        'original_name': filename,
        'upload_path': upload_path,
        'file_extension': file_extension,
    }


def _output_name_for(original_filename):
    """出力 h5mu のファイル名を生成する: ``convert_<basename>.h5mu`` 形式。
    ``.tar.gz`` や最後の拡張子を取り除いてベース名を作る。"""
    base = original_filename
    if base.lower().endswith('.tar.gz'):
        base = base[:-7]
    elif '.' in base:
        base = base.rsplit('.', 1)[0]
    return f"convert_{base}.h5mu"


def _infer_file_type(filename):
    """一括変換の自動モード時に、拡張子からファイル種別を推定する。"""
    name = filename.lower()
    if name.endswith('.csv'):
        return 'csv'
    if name.endswith('.h5ad'):
        return 'h5ad'
    if name.endswith('.rds'):
        return 'rds'
    if name.endswith('.zip') or name.endswith('.tar.gz') or name.endswith('.tar'):
        return 'mex'
    return None


def _run_conversion(upload_path, output_path, file_type, filename, file_extension, options):
    """単一ファイルの変換を実行する。失敗時は日本語メッセージ付きの例外を投げる。

    Returns
    -------
    dict
        変換に関する追加情報（自動転置判定の結果など）。
        特に報告すべき情報がなければ空の dict。
    """
    extra_info = {}
    if file_type == 'csv':
        result = csv_to_h5mu(
            upload_path,
            output_path,
            transpose=options.get('transpose', 'auto'),
            has_header=options.get('has_header', True),
            has_index=options.get('has_index', True),
        )
        if isinstance(result, dict):
            extra_info = {
                'auto_detected_transpose': result.get('auto_detected', False),
                'transpose_applied': result.get('transpose_applied', False),
                'detect_info': result.get('detect_info', {}),
            }

    elif file_type == 'h5ad':
        h5ad_to_h5mu(upload_path, output_path)

    elif file_type == 'rds':
        if not check_rpy2_available():
            raise RuntimeError(
                'RDS変換にはrpy2が必要です。インストール方法については README.md を参照してください。'
            )
        rds_to_h5mu(upload_path, output_path, object_type='auto')

    elif file_type == 'mex':
        extract_dir = upload_path + '_extracted'
        if filename.lower().endswith('.tar.gz'):
            os.makedirs(extract_dir, exist_ok=True)
            try:
                with tarfile.open(upload_path, 'r:gz') as tar_ref:
                    tar_ref.extractall(extract_dir)
                mex_dir = find_mex_directory(extract_dir)
                if not mex_dir:
                    raise RuntimeError(
                        'TAR.GZファイル内にMEXファイル（matrix.mtx、features.tsv、barcodes.tsv）が見つかりません。'
                    )
                mex_to_h5mu(mex_dir, output_path)
            finally:
                if os.path.isdir(extract_dir):
                    shutil.rmtree(extract_dir, ignore_errors=True)

        elif file_extension == 'zip':
            os.makedirs(extract_dir, exist_ok=True)
            try:
                with zipfile.ZipFile(upload_path, 'r') as zip_ref:
                    zip_ref.extractall(extract_dir)
                mex_dir = find_mex_directory(extract_dir)
                if not mex_dir:
                    raise RuntimeError(
                        'ZIPファイル内にMEXファイル（matrix.mtx、features.tsv、barcodes.tsv）が見つかりません。'
                    )
                mex_to_h5mu(mex_dir, output_path)
            finally:
                if os.path.isdir(extract_dir):
                    shutil.rmtree(extract_dir, ignore_errors=True)

        elif file_extension == 'gz':
            raise RuntimeError(
                'MEXファイルは、matrix.mtx、features.tsv、barcodes.tsvを含むZIPまたはTAR.GZファイルとしてアップロードしてください。'
            )
        else:
            raise RuntimeError('MEXファイルはZIPまたはTAR.GZ形式でアップロードしてください。')

    else:
        raise RuntimeError('不明なファイル形式です')

    return extra_info


# ── バックグラウンド UMAP ジョブの管理 ────────────────────────────────────
# 大きなデータでは UMAP 計算に数分かかるため、ワーカースレッドで実行し、
# 現在のステージを /umap_progress/<job_id> 経由で UI に公開する。
# ユーザーが「フリーズした」と誤解しないようリアルタイム更新を可能にする。

_UMAP_JOBS_LOCK = threading.Lock()
_UMAP_JOBS: 'OrderedDict[str, dict]' = OrderedDict()
_UMAP_JOBS_MAX = 50  # メモリリーク防止: 保持するジョブの最大件数


def _job_set(job_id, **fields):
    """既存ジョブのフィールドを安全に更新する。"""
    with _UMAP_JOBS_LOCK:
        job = _UMAP_JOBS.get(job_id)
        if job is None:
            return
        job.update(fields)
        job['updated_at'] = time.time()


def _job_create():
    """新しいジョブ ID を発行し、初期状態で登録する。"""
    job_id = str(uuid.uuid4())[:12]
    with _UMAP_JOBS_LOCK:
        if len(_UMAP_JOBS) >= _UMAP_JOBS_MAX:
            # 上限に達したら一番古いジョブから捨てる
            _UMAP_JOBS.popitem(last=False)
        _UMAP_JOBS[job_id] = {
            'status': 'running',
            'stage': '開始中...',
            'started_at': time.time(),
            'updated_at': time.time(),
        }
    return job_id


def _job_get(job_id):
    """ジョブ情報のスナップショットを返す（lock で防御）。"""
    with _UMAP_JOBS_LOCK:
        job = _UMAP_JOBS.get(job_id)
        if job is None:
            return None
        return dict(job)


def _run_umap_job(job_id, h5mu_path, output_png_path):
    """別スレッドで UMAP 計算を実行する本体。"""
    def progress_cb(msg):
        # umap_visualizer が呼ぶコールバック。現在のステージ名を保存する。
        _job_set(job_id, stage=msg)

    try:
        info = umap_from_h5mu(h5mu_path, output_png_path, progress_cb=progress_cb)
        umap_filename = os.path.basename(output_png_path)
        _job_set(
            job_id,
            status='done',
            stage='完了',
            result={
                'download_url': f'/download/{umap_filename}',
                'preview_url': f'/preview/{umap_filename}',
                'modality': info.get('modality'),
                'n_cells': info.get('n_cells'),
                'n_genes': info.get('n_genes'),
                'color_by': info.get('color_by'),
            },
        )
    except Exception as e:
        traceback.print_exc()
        _job_set(job_id, status='error', stage='エラー', error=str(e))


def _start_umap_job(h5mu_path, output_png_path):
    """UMAP 計算をバックグラウンドスレッドで起動し、ジョブ ID を返す。"""
    job_id = _job_create()
    t = threading.Thread(
        target=_run_umap_job,
        args=(job_id, h5mu_path, output_png_path),
        daemon=True,  # メインプロセス終了時に道連れにする
    )
    t.start()
    return job_id


@app.errorhandler(413)
def _too_large(_):
    """MAX_CONTENT_LENGTH を超えた時に分かりやすい日本語メッセージを返す。"""
    limit_mb = app.config['MAX_CONTENT_LENGTH'] / (1024 * 1024)
    return jsonify({
        'error': f'ファイルが大きすぎます（上限 {limit_mb:.0f} MB）。'
                 f'環境変数 MAX_CONTENT_LENGTH をバイト単位で指定して上限を変更できます。'
    }), 413


@app.route('/')
def index():
    """トップページ。"""
    return render_template('index.html')


@app.route('/convert', methods=['POST'])
def convert():
    """単一ファイルの変換リクエストを処理する。"""
    upload_path = None
    try:
        if 'file' not in request.files:
            return jsonify({'error': 'ファイルが選択されていません'}), 400

        file = request.files['file']
        if file.filename == '':
            return jsonify({'error': 'ファイルが選択されていません'}), 400

        if not allowed_file(file.filename):
            return jsonify({'error': '対応していないファイル形式です'}), 400

        file_type = request.form.get('file_type', 'csv')
        if file_type in (None, '', 'auto'):
            file_type = _infer_file_type(file.filename) or 'csv'

        # transpose は 'auto'（既定）/ 'true' / 'false' のいずれかを受け付ける。
        # 'auto' のときは csv_converter 側で自動判定する。
        transpose_raw = request.form.get('transpose', 'auto')
        if transpose_raw in ('true', 'false'):
            transpose_val = transpose_raw == 'true'
        else:
            transpose_val = 'auto'

        options = {
            'transpose': transpose_val,
            'has_header': request.form.get('has_header', 'true') == 'true',
            'has_index': request.form.get('has_index', 'true') == 'true',
        }
        make_umap = request.form.get('make_umap', 'false') == 'true'

        saved = _save_upload(file)
        upload_path = saved['upload_path']
        output_filename = _output_name_for(saved['original_name'])
        output_path = os.path.join(app.config['OUTPUT_FOLDER'], output_filename)

        extra_info = _run_conversion(
            upload_path=upload_path,
            output_path=output_path,
            file_type=file_type,
            filename=saved['original_name'],
            file_extension=saved['file_extension'],
            options=options,
        ) or {}

        response = {
            'success': True,
            'message': '変換が完了しました',
            'download_url': f'/download/{output_filename}',
        }

        if extra_info.get('auto_detected_transpose'):
            response['detection'] = {
                'auto_detected': True,
                'transpose_applied': extra_info.get('transpose_applied', False),
                'reason': extra_info.get('detect_info', {}).get('reason', ''),
                'scores': {
                    k: v for k, v in extra_info.get('detect_info', {}).items()
                    if k in ('row_gene_score', 'col_gene_score', 'row_barcode_score', 'col_barcode_score', 'shape')
                },
            }

        # オプション: UMAP をバックグラウンドジョブとして起動する。
        # （大きなデータでは PCA → neighbors → UMAP → Leiden → plot で
        # 数分かかるため、UI 側がポーリングして進捗を表示できるよう
        # job_id を返す設計にしている。）
        if make_umap:
            try:
                umap_filename = output_filename.replace('.h5mu', '_umap.png')
                umap_path = os.path.join(app.config['OUTPUT_FOLDER'], umap_filename)
                job_id = _start_umap_job(output_path, umap_path)
                response['umap_job'] = {
                    'job_id': job_id,
                    'progress_url': f'/umap_progress/{job_id}',
                }
                response['message'] = '変換が完了しました。UMAP をバックグラウンドで生成中です...'
            except Exception as ue:
                traceback.print_exc()
                response['umap_job'] = {
                    'error': f'UMAP ジョブの起動に失敗しました: {str(ue)}',
                }

        return jsonify(response)

    except Exception as e:
        traceback.print_exc()
        return jsonify({'error': f'変換エラー: {str(e)}'}), 500
    finally:
        if upload_path and os.path.exists(upload_path):
            try:
                os.remove(upload_path)
            except OSError:
                pass


@app.route('/umap_progress/<job_id>')
def umap_progress(job_id):
    """UI のポーリング用エンドポイント: UMAP ジョブの現在状態を返す。"""
    job = _job_get(job_id)
    if job is None:
        return jsonify({'error': 'ジョブが見つかりません'}), 404
    elapsed = time.time() - job.get('started_at', time.time())
    return jsonify({
        'success': True,
        'status': job.get('status'),
        'stage': job.get('stage'),
        'elapsed_sec': round(elapsed, 1),
        'result': job.get('result'),
        'error': job.get('error'),
    })


@app.route('/preview/<filename>')
def preview(filename):
    """インライン画像プレビュー用（UMAP をダウンロードせず画面表示する用途）。"""
    try:
        safe = secure_filename(filename)
        filepath = os.path.join(app.config['OUTPUT_FOLDER'], safe)
        if not os.path.exists(filepath):
            return jsonify({'error': 'ファイルが見つかりません'}), 404
        return send_file(filepath, mimetype='image/png')
    except Exception as e:
        return jsonify({'error': f'プレビューエラー: {str(e)}'}), 500


@app.route('/convert_batch', methods=['POST'])
def convert_batch():
    """複数ファイルを一括変換し、ZIP にまとめてダウンロードリンクを返す。"""
    file_type_form = request.form.get('file_type', 'auto')

    transpose_raw = request.form.get('transpose', 'auto')
    if transpose_raw in ('true', 'false'):
        transpose_val = transpose_raw == 'true'
    else:
        transpose_val = 'auto'

    options = {
        'transpose': transpose_val,
        'has_header': request.form.get('has_header', 'true') == 'true',
        'has_index': request.form.get('has_index', 'true') == 'true',
    }

    files = request.files.getlist('files')
    if not files or all(f.filename == '' for f in files):
        return jsonify({'error': 'ファイルが選択されていません'}), 400

    results = []
    produced_outputs = []
    upload_paths = []

    try:
        for f in files:
            if f.filename == '':
                continue
            if not allowed_file(f.filename):
                results.append({
                    'filename': f.filename,
                    'success': False,
                    'error': '対応していないファイル形式です',
                })
                continue

            saved = _save_upload(f)
            upload_paths.append(saved['upload_path'])

            ftype = file_type_form
            if ftype in (None, '', 'auto'):
                ftype = _infer_file_type(saved['original_name'])
                if ftype is None:
                    results.append({
                        'filename': saved['original_name'],
                        'success': False,
                        'error': 'ファイル種別を判定できませんでした',
                    })
                    continue

            output_filename = _output_name_for(saved['original_name'])
            output_path = os.path.join(app.config['OUTPUT_FOLDER'], output_filename)

            try:
                _run_conversion(
                    upload_path=saved['upload_path'],
                    output_path=output_path,
                    file_type=ftype,
                    filename=saved['original_name'],
                    file_extension=saved['file_extension'],
                    options=options,
                )
                produced_outputs.append((output_filename, output_path))
                results.append({
                    'filename': saved['original_name'],
                    'success': True,
                    'output': output_filename,
                    'file_type': ftype,
                })
            except Exception as conv_err:
                traceback.print_exc()
                results.append({
                    'filename': saved['original_name'],
                    'success': False,
                    'error': str(conv_err),
                    'file_type': ftype,
                })

        if not produced_outputs:
            return jsonify({
                'success': False,
                'error': 'すべての変換に失敗しました',
                'results': results,
            }), 500

        zip_name = f"batch_{datetime.now().strftime('%Y%m%d_%H%M%S')}_{str(uuid.uuid4())[:8]}.zip"
        zip_path = os.path.join(app.config['OUTPUT_FOLDER'], zip_name)
        with zipfile.ZipFile(zip_path, 'w', compression=zipfile.ZIP_DEFLATED) as zf:
            for out_name, out_path in produced_outputs:
                zf.write(out_path, arcname=out_name)

        # 個別の h5mu ファイルは OUTPUT_FOLDER に残しておく
        # （ZIP はバッチ全体の一括ダウンロード用、個別ファイルは将来の
        # 個別操作用に保持する）。

        return jsonify({
            'success': True,
            'message': f'{len(produced_outputs)} / {len(results)} 件の変換が完了しました',
            'download_url': f'/download/{zip_name}',
            'results': results,
        })

    finally:
        for p in upload_paths:
            if os.path.exists(p):
                try:
                    os.remove(p)
                except OSError:
                    pass


@app.route('/download/<filename>')
def download(filename):
    """変換結果ファイルのダウンロードエンドポイント。"""
    try:
        filepath = os.path.join(app.config['OUTPUT_FOLDER'], secure_filename(filename))
        if not os.path.exists(filepath):
            return jsonify({'error': 'ファイルが見つかりません'}), 404

        return send_file(
            filepath,
            as_attachment=True,
            download_name=filename,
            mimetype='application/octet-stream'
        )
    except Exception as e:
        return jsonify({'error': f'ダウンロードエラー: {str(e)}'}), 500


@app.route('/validate', methods=['POST'])
def validate():
    """アップロードファイルの簡易バリデーション（現状 CSV のみ実装）。"""
    try:
        if 'file' not in request.files:
            return jsonify({'error': 'ファイルが選択されていません'}), 400

        file = request.files['file']
        if file.filename == '':
            return jsonify({'error': 'ファイルが選択されていません'}), 400

        # 一時保存して中身を確認する
        filename = secure_filename(file.filename)
        unique_name = generate_unique_filename('validate')
        file_extension = filename.rsplit('.', 1)[1].lower()
        temp_path = os.path.join(app.config['UPLOAD_FOLDER'], f"{unique_name}.{file_extension}")
        file.save(temp_path)

        # ファイル種別ごとに対応するバリデータを呼ぶ
        file_type = request.form.get('file_type', 'csv')

        if file_type == 'csv':
            info = validate_csv(temp_path)
        else:
            info = {'message': 'Validation not yet implemented for this file type'}

        # 一時ファイルを削除
        os.remove(temp_path)

        return jsonify({'success': True, 'info': info})

    except Exception as e:
        return jsonify({'error': f'検証エラー: {str(e)}'}), 500


# アプリ起動時に必要なディレクトリを作成
os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)
os.makedirs(app.config['OUTPUT_FOLDER'], exist_ok=True)


if __name__ == '__main__':
    # 環境変数から設定を取得
    # 注意: macOS Monterey 以降は AirPlay Receiver がポート 5000 を占有しているため、
    # 既定は 5050 を使う。PORT 環境変数で上書き可能。
    port = int(os.getenv('PORT', 5050))
    host = os.getenv('HOST', '0.0.0.0')
    debug = os.getenv('FLASK_ENV', 'development') == 'development'

    # Flask 開発サーバーを起動
    print("Starting Single-cell File Converter Web Application...")
    print(f"Access the application at: http://{host}:{port}")
    print(f"Debug mode: {debug}")
    app.run(debug=debug, host=host, port=port)
