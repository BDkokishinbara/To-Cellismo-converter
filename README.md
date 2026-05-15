# 🧬 To-Cellismo Converter

シングルセル解析ファイルを簡単に変換！CSV / h5ad / RDS / MEX から h5mu 形式への変換ツール（**複数ファイル一括変換対応**）

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

---

## 📖 目次

- [クイックスタート](#クイックスタート)
- [このツールについて](#このツールについて)
- [対応ファイル形式](#対応ファイル形式)
- [インストール方法](#インストール方法)
- [使い方](#使い方)
  - [Webアプリを使う](#webアプリを使う)
  - [コマンドラインで使う](#コマンドラインで使う)
- [h5mu形式とは](#h5mu形式とは)
- [変換したファイルの使い方](#変換したファイルの使い方)
- [トラブルシューティング](#トラブルシューティング)
- [よくある質問](#よくある質問)

---

## 🚀 クイックスタート（IT初学者向け・推奨）

**ダブルクリックするだけで動きます。お使いの PC の Python 環境は一切汚しません。**

### ステップ1: ダウンロード

GitHub の「**Code**」ボタン → 「**Download ZIP**」でダウンロードして解凍してください。
（コマンドラインに慣れている方は `git clone https://github.com/BDkokishinbara/To-Cellismo-converter.git` でも OK）

### ステップ2: 起動スクリプトをダブルクリック

| OS | ダブルクリックするファイル |
|---|---|
| **macOS** | `start.command` |
| **Windows** | `start.bat` |
| **Linux** | `start.sh` |

> 💡 **macOS の方へ**: `.sh` ファイルは Cursor や VSCode などのエディタで開いてしまうため、macOS では必ず **`start.command`** をダブルクリックしてください。中身は `start.sh` と同じ動作をします。
>
> 初回ダブルクリック時に「**開発元を確認できないため開けません**」と表示された場合は、`start.command` を **右クリック → 開く → 開く** で許可してください（次回からは普通にダブルクリックで起動できます）。

### それだけ！

スクリプトが自動で:

1. 専用の仮想環境（プロジェクト内 `.cellismo_venv/` フォルダ）を作成
2. 必要なパッケージをそこにだけインストール
3. アプリを起動し、ブラウザで **http://127.0.0.1:5050** を開く

**システム側の Python やパッケージには一切手を加えません。**
不要になったら、ダウンロードした **フォルダごと削除するだけで完全に元通り** になります。

> 💡 **事前に必要なもの**: Python 3.8 以降だけ（https://www.python.org からインストール）。
> Python さえ入っていれば、それ以外は起動スクリプトがすべて自動でセットアップします。

停止するには、ターミナル（または黒い画面）で `Ctrl+C` を押してください。

---

### 🔁 2回目以降の起動

**初回と全く同じ — `start.command`（macOS）／ `start.bat`（Windows）／ `start.sh`（Linux）をダブルクリックするだけです。**

| 項目 | 初回起動 | 2回目以降 |
|---|---|---|
| 仮想環境の作成 | 自動で行う（数十秒） | スキップ（既にある） |
| パッケージのインストール | 自動で行う（1〜数分） | スキップ |
| アプリ起動 | 行う | 行う |
| 起動までの時間 | 1〜数分 | **数秒** |

つまり、2回目以降は **ダブルクリックして数秒待つとブラウザが開く** だけです。
特別な操作は何も必要ありません。

#### 停止のしかた

- ターミナル（黒い画面）で `Ctrl+C`
- または、Windows の場合はそのウィンドウを閉じる

#### 完全にゼロから作り直したいとき

何かおかしくなった、最初からやり直したい、というときは:

1. プロジェクトフォルダの中の **`.cellismo_venv/` フォルダを削除**
2. もう一度 `start.command` / `start.bat` / `start.sh` をダブルクリック

→ 初回と同じように、仮想環境とパッケージが新しく作られます。

> 💡 ホームフォルダや他のプロジェクトの `.venv` には影響しません。このアプリ専用の `.cellismo_venv` フォルダだけが対象です。

#### パッケージを最新版に更新したいとき

```bash
# プロジェクトフォルダで
source .cellismo_venv/bin/activate            # macOS / Linux （Windows は .cellismo_venv\Scripts\activate）
pip install --upgrade -r requirements.txt
```

または、`.cellismo_venv/` フォルダを丸ごと削除してから起動スクリプトを実行し直しても OK です。

---

### ⌨️ 手動で起動したい場合（中〜上級者向け）

仮想環境を自分で管理したい場合は次のように:

```bash
# 仮想環境を作成して有効化
python -m venv .venv
source .venv/bin/activate            # macOS / Linux
# .venv\Scripts\activate            # Windows

# 依存パッケージをインストール（venv 内だけに入ります）
pip install -r requirements.txt

# 起動
python app.py
```

ブラウザで **http://127.0.0.1:5050** にアクセスしてください。停止は `Ctrl+C`。

---

## このツールについて

### 何ができるの？

このツールは、シングルセル解析で使われる様々なファイル形式を、**h5mu（MuData）形式**という統一されたファイル形式に変換します。

### なぜ必要なの？

シングルセル解析では、実験機器や解析ソフトウェアによって異なるファイル形式が使われます：
- CSV形式のカウントマトリックス
- RDS形式のSeuratオブジェクト
- 10x Genomicsのスパースマトリックス

これらを **h5mu形式に統一** することで：
✅ データの管理が簡単になる
✅ PythonでもRでも同じファイルを使える
✅ ファイルサイズが小さくなる
✅ 複数のデータセットを1つのファイルにまとめられる

---

## 対応ファイル形式

### 📊 CSV形式
- **説明**: 表形式のテキストファイル（Excelで開けるもの）
- **使用例**: SeqGeqやカスタム解析で出力されるカウントマトリックス
- **形式**:
  - 細胞×遺伝子（横に細胞、縦に遺伝子）
  - 遺伝子×細胞（横に遺伝子、縦に細胞）
- **特別対応**: SeqGeq形式（メタデータセクション付き）を自動認識

### 🧪 h5ad形式
- **説明**: AnnData (`.h5ad`) — Scanpy などで広く使われているシングルセルデータ形式
- **使用例**: `sc.write('data.h5ad', adata)` で保存したファイル
- **動作**: AnnData を `rna` モダリティとして MuData にラップして h5mu に書き出します
- **注意**: 既存の `obs` / `var` メタデータはそのまま保持されます

### 📦 RDS形式
- **説明**: R言語で保存されたオブジェクト
- **使用例**: Seurat、SingleCellExperimentなどのRパッケージで作成されたデータ
- **注意**: rpy2のインストールが必要（後述）

### 🧬 MEX形式
- **説明**: 10x Genomics社のスパースマトリックス形式
- **構成**: 3つのファイルが必要
  - `matrix.mtx` または `matrix.mtx.gz` - データ本体
  - `features.tsv` または `genes.tsv` - 遺伝子名
  - `barcodes.tsv` - 細胞バーコード
- **アップロード**: ZIP形式でまとめてアップロード

---

## インストール方法

### 🔧 必要なもの

- **Python 3.8以上**（https://www.python.org からインストール）
- **インターネット接続**（初回インストール時のみ）

### おすすめ: 起動スクリプトに任せる

[クイックスタート](#-クイックスタートit初学者向け推奨)の通り、`start.sh`（macOS/Linux）または `start.bat`（Windows）をダブルクリックすれば、仮想環境の作成・パッケージのインストール・アプリ起動まで全自動で行われます。
**特別な準備は何も要りません。**

### 手動でインストールする場合（中〜上級者向け）

仮想環境を作って、その中だけにパッケージを入れます（システムの Python は汚れません）。

```bash
# 1. リポジトリを取得
git clone https://github.com/BDkokishinbara/To-Cellismo-converter.git
cd To-Cellismo-converter

# 2. プロジェクト専用の仮想環境を作成
python -m venv .venv

# 3. 仮想環境を有効化
source .venv/bin/activate            # macOS / Linux
# .venv\Scripts\activate            # Windows

# 4. 仮想環境内に必要なパッケージをインストール
pip install -r requirements.txt
```

> 💡 仮想環境は単なるフォルダ（`.venv/`）です。気に入らなければ削除するだけで元に戻ります。

### オプション: RDSファイル変換を使う場合

RDSファイルを変換したい場合は、追加で R 本体と rpy2 をインストールしてください：

```bash
# Rのインストール（Ubuntuの場合）
sudo apt-get install r-base r-base-dev

# 仮想環境を有効化したうえで（または起動スクリプトで一度起動した後で）
pip install rpy2
```

**注意**: RDSファイルを使わない場合、この手順はスキップできます。

---

## 使い方

### 🌐 Webアプリを使う（推奨）

#### 1. アプリを起動

**かんたん起動（推奨）**: `start.sh`（macOS/Linux）または `start.bat`（Windows）を**ダブルクリック**するだけ。
初回は仮想環境の作成と依存インストールが入るので少し時間がかかります。2回目以降はすぐ起動します。

**手動起動（起動スクリプトが作った仮想環境を直接使う場合）**:

```bash
cd To-Cellismo-converter
source .cellismo_venv/bin/activate    # macOS / Linux （Windows は .cellismo_venv\Scripts\activate）
python app.py
```

起動すると以下のようなメッセージが表示されます：

```
Starting Single-cell File Converter Web Application...
Access the application at: http://0.0.0.0:5050
```

#### 2. ブラウザでアクセス

ブラウザで以下のURLを開きます：

**http://127.0.0.1:5050**

#### 3. ファイルを変換

1. **ファイル形式を選択**: 自動 / CSV / h5ad / RDS / MEX から選ぶ
   - **自動**: 拡張子から判定するモード。複数ファイル一括変換にも最適
2. **（任意）複数ファイルを一括変換する**にチェック
   - 複数ファイルをドラッグ&ドロップでまとめてアップロードできます
   - 異なる形式（CSVとh5adなど）が混在していても、自動モードなら拡張子で振り分けて変換します
3. **ファイルをアップロード**: ドラッグ&ドロップまたはクリックして選択
4. **オプションを設定**（CSVの場合）:
   - **🔄 転置は自動判定**: 遺伝子名（ACTB, CD3D, ENSG... など）が行と列のどちらにあるかを検出して、必要に応じて自動で転置します。判定結果は変換完了時に画面に表示されます
   - **Has header**: 1行目がヘッダー（細胞名）の場合にチェック（通常はオン）
   - **Has index**: 1列目がインデックス（遺伝子名）の場合にチェック（通常はオン）
5. **（任意）「変換と同時に UMAP も作成する」にチェック**
   - scanpy で標準的な前処理（filter → normalize → log1p → HVG → PCA → neighbors → Leiden → UMAP）を実行
   - 数十秒〜数分かかります。一括変換ではオフ推奨
6. **「h5mu に変換」ボタンをクリック**
7. **自動的にダウンロード開始**:
   - 単一ファイル: `convert_元のファイル名.h5mu`
   - 一括変換: `batch_<日時>_<ID>.zip`（中に各 `convert_*.h5mu` が含まれます）
   - UMAP オン時: `convert_元のファイル名_umap.png` も同時に作成され、画面にプレビュー表示されます


#### 4. アプリを停止

ターミナルで `Ctrl+C` を押すと停止します。

---

### 💻 コマンドラインで使う

#### CSVファイルの変換

```python
from converters.csv_converter import csv_to_h5mu

# 基本的な使い方（転置は自動判定）
csv_to_h5mu(
    csv_path='your_data.csv',
    output_path='output.h5mu',
    transpose='auto',     # 'auto'（既定）/ True / False
    has_header=True,
    has_index=True
)
```

返り値は dict で、`auto_detected` / `transpose_applied` / `detect_info` から判定根拠を確認できます。

#### SeqGeq形式のCSVを変換

```python
from converters.csv_converter import csv_to_h5mu

# SeqGeq形式は自動認識されます
csv_to_h5mu(
    csv_path='SeqGeq_demo.csv',
    output_path='SeqGeq_demo.h5mu',
    transpose=False,
    has_header=True,
    has_index=True
)
```

#### RDSファイルの変換

```python
from converters.rds_converter import rds_to_h5mu

rds_to_h5mu(
    rds_path='seurat_object.rds',
    output_path='output.h5mu',
    object_type='auto'  # 'seurat' または 'sce' も指定可能
)
```

#### h5adファイルの変換

```python
from converters.h5ad_converter import h5ad_to_h5mu

h5ad_to_h5mu(
    h5ad_path='data.h5ad',
    output_path='output.h5mu'
)
```

#### MEXファイルの変換

```python
from converters.mex_converter import mex_to_h5mu

mex_to_h5mu(
    mex_directory='path/to/mex_folder',  # matrix.mtx等があるフォルダ
    output_path='output.h5mu'
)
```

---

## h5mu形式とは

### 📚 MuData (h5mu) について

**h5mu**は、**MuData (Multimodal Data)** 形式のファイルです。

#### 特徴

✅ **効率的**: 大きなデータを圧縮して保存
✅ **マルチモーダル**: RNA、ATAC、Proteinなど複数のデータを1つのファイルに
✅ **互換性**: PythonでもRでも使える
✅ **標準化**: シングルセル解析コミュニティで広く使われている

#### 構造

```
h5mu ファイル
├── rna (RNA-seqデータ)
│   ├── X (カウントマトリックス)
│   ├── obs (細胞のメタデータ)
│   └── var (遺伝子のメタデータ)
├── atac (ATAC-seqデータ) ※オプション
└── protein (タンパク質データ) ※オプション
```

---

## 変換したファイルの使い方

### 🐍 Pythonで使う

#### 基本的な読み込み

```python
import mudata as md
import scanpy as sc

# h5muファイルを読み込む
mdata = md.read_h5mu('convert_your_data.h5mu')

# RNAデータを取り出す
adata = mdata['rna']

# データの情報を確認
print(f"細胞数: {adata.n_obs}")
print(f"遺伝子数: {adata.n_vars}")
print(adata)
```

#### Scanpyで解析

```python
import scanpy as sc

# 標準的な前処理
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# 高変動遺伝子の検出
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

# 次元削減
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)

# クラスタリング
sc.tl.leiden(adata)

# 可視化
sc.pl.umap(adata, color='leiden')
```

### 📊 Rで使う

```r
library(MuData)
library(Seurat)

# h5muファイルを読み込む
mdata <- readH5MU('convert_your_data.h5mu')

# Seuratオブジェクトに変換（オプション）
# 注意: この機能は使用するパッケージのバージョンによります
```

---

## トラブルシューティング

### ❌ よくあるエラーと対処法

#### エラー0: pandas/numpy のビルドエラー（`_PyLong_AsByteArray`, `_PyDict_SetItem_KnownHash` などのエラーが出る）

**原因**: お使いの Python が **最新すぎる**（例: Python 3.14）。古いバージョンの pandas/numpy はまだ Python 3.14 用のビルド済みファイル（wheel）を提供しておらず、ソースからビルドしようとして失敗します。

**解決方法**:

1. プロジェクトフォルダ内の **`.cellismo_venv/` フォルダを削除**（失敗した仮想環境を片付ける）。
   起動スクリプトを使っていれば自動で削除されているはずです。
2. もう一度 `start.command` / `start.bat` をダブルクリック。
   起動スクリプトは Python 3.13 / 3.12 / 3.11 を優先して使うので、それらが入っていればそちらが選ばれます。
3. それでもダメな場合は、Python 3.13 など安定版をインストールしてください:
   - macOS: `brew install python@3.13`
   - Windows: https://www.python.org/downloads/ から Python 3.13 をインストール

なお、本プロジェクトの `requirements.txt` はバージョン固定をゆるめてあるため、新しい Python でも対応する新バージョンの pandas/numpy が自動で選ばれるようになっています（古いバージョンに固執しません）。

---

#### エラー1: `ModuleNotFoundError: No module named 'flask'`

**原因**: 必要なパッケージがインストールされていない、または仮想環境が有効化されていない

**解決方法**:

1. **起動スクリプトを使う方（一番簡単）**: `start.command`（macOS）／ `start.bat`（Windows）／ `start.sh`（Linux）をダブルクリックすれば自動でインストールされます。
2. **手動起動の方**: 仮想環境を有効化してからインストールしてください。
   ```bash
   source .cellismo_venv/bin/activate    # macOS / Linux （Windows は .cellismo_venv\Scripts\activate）
   pip install -r requirements.txt
   ```

#### エラー2: `ParserError: Error tokenizing data`

**原因**: CSVファイルの形式が想定と異なる

**解決方法**:
- ファイルの最初の数行を確認してください
- ヘッダー行やインデックス列の設定を変更してみてください
- SeqGeq形式の場合、自動認識されているか確認してください

#### エラー3: `rpy2 not available`

**原因**: RDSファイル変換に必要なrpy2がインストールされていない

**解決方法**:
```bash
# Rをインストール（Ubuntuの場合）
sudo apt-get install r-base r-base-dev

# rpy2をインストール
pip install rpy2
```

#### エラー4: `Permission denied`

**原因**: ファイルの読み書き権限がない

**解決方法**:
```bash
# uploadsとoutputsフォルダに書き込み権限を付与
chmod 755 uploads outputs
```

#### エラー5: ブラウザで「**127.0.0.1 へのアクセスが拒否されました（HTTP ERROR 403）**」と出る

**原因**: macOS Monterey 以降は **AirPlay Receiver** がポート 5000 を占有しているため、アプリが起動できない、または AirPlay 側が代わりに 403 を返してしまうことがあります。

**解決方法**: 本アプリは既定で **ポート 5050** を使うようになっています。`start.command` を最新版で実行すれば、ブラウザは自動的に http://127.0.0.1:5050 を開きます。

もし古いバージョンのスクリプトを使っていた / 別のポートを試したい場合は:

```bash
# 別のポートを使用する（例: 8080）
export PORT=8080
python app.py
```

または、macOS の AirPlay Receiver を無効化することでもポート 5000 が解放されます:
**システム設定 → 一般 → AirDrop と Handoff → AirPlay レシーバー** をオフ。

#### エラー6: `Port XXXX already in use`

**原因**: そのポートが既に他のプロセスに使われている

**解決方法**:
```bash
# 別のポートを使用する
export PORT=8080
python app.py
```

---

## よくある質問

### Q1: どのくらいのファイルサイズまで対応していますか？

**A**: デフォルトでは **5 GB** までです。シングルセルの h5ad/h5mu は GB 級になることが多いため、余裕を持たせています。
変更したい場合は環境変数 `MAX_CONTENT_LENGTH`（バイト単位）で指定できます。例:

```bash
# 10 GB に拡張
export MAX_CONTENT_LENGTH=10737418240
python app.py
```

### Q2: SeqGeq形式とは何ですか？

**A**: SeqGeqソフトウェアが出力するCSV形式で、ファイルの先頭にメタデータセクション（`[Metadata]`と`[Data]`）があります。このツールは自動的に認識して適切に処理します。

### Q3: 変換にどのくらい時間がかかりますか？

**A**: ファイルサイズによりますが、目安は以下の通りです：
- 小規模（〜1万細胞）: 数秒〜1分
- 中規模（1万〜10万細胞）: 1分〜10分
- 大規模（10万細胞以上）: 10分以上

### Q4: 変換したファイルはどこに保存されますか？

**A**:
- Webアプリ: ブラウザのダウンロードフォルダ
- コマンドライン: 指定した`output_path`

### Q5: 元のデータは失われますか？

**A**: いいえ。このツールは元のファイルを読み込んで新しいファイルを作成するだけで、元のファイルは変更されません。

### Q6: 複数のファイルを一度に変換できますか？

**A**: はい。Webアプリで「複数ファイルを一括変換する」にチェックを入れて、複数ファイルをまとめてアップロードしてください。各ファイルが個別に変換され、結果は1つのZIPファイルにまとめてダウンロードされます。形式が混在していても、ファイル形式を「自動」にしておけば拡張子から判定されます。

### Q7: オフラインで使えますか？

**A**: はい。一度必要なパッケージをインストールすれば、インターネット接続なしで使用できます。

### Q8: 商用利用できますか？

**A**: はい。MITライセンスの下で、商用利用を含めて自由に使用できます。

---

## 📝 ファイル構成

```
To-Cellismo-converter/
├── app.py                               # Flask アプリケーション本体
├── requirements.txt                     # 必要な Python パッケージ
├── README.md                            # このファイル
├── start.command                        # macOS 用ダブルクリック起動
├── start.sh                             # Linux 用起動スクリプト
├── start.bat                            # Windows 用起動スクリプト
├── converters/                          # 変換ロジック
│   ├── csv_converter.py                # CSV 変換（SeqGeq 対応 / 自動転置判定）
│   ├── h5ad_converter.py               # h5ad (AnnData) 変換
│   ├── rds_converter.py                # RDS 変換
│   ├── mex_converter.py                # MEX 変換
│   └── umap_visualizer.py              # 変換後 UMAP 可視化
├── templates/                           # HTML テンプレート
│   └── index.html
├── static/                              # CSS / JavaScript
│   ├── style.css
│   └── script.js
├── uploads/                             # アップロードファイル一時保存
└── outputs/                             # 変換結果（h5mu / UMAP PNG）
```

---

## 🔗 関連リンク

### ドキュメント
- [MuData Documentation](https://mudata.readthedocs.io/) - MuData形式の公式ドキュメント
- [Scanpy Documentation](https://scanpy.readthedocs.io/) - Pythonでのシングルセル解析
- [Seurat Website](https://satijalab.org/seurat/) - Rでのシングルセル解析

### リポジトリ
- [GitHub](https://github.com/BDkokishinbara/To-Cellismo-converter)
- [Hugging Face Space](https://huggingface.co/spaces/BDKoki/To-Cellismo-converter)

---

## 🤝 貢献

バグ報告、機能リクエスト、プルリクエストを歓迎します！

---

## 📄 ライセンス

MIT License - 詳細は[LICENSE](LICENSE)ファイルをご覧ください。

---

## 🙏 謝辞

このツールは以下のオープンソースプロジェクトを使用しています：

- **Flask** - Webアプリケーションフレームワーク
- **AnnData** - 注釈付き行列データ構造
- **MuData** - マルチモーダルデータ管理
- **Scanpy** - シングルセル解析ツールキット
- **pandas** - データ操作
- **NumPy** - 数値計算
- **SciPy** - 科学計算

---

## 📮 お問い合わせ

質問や問題がある場合は、[GitHubのIssues](https://github.com/BDkokishinbara/To-Cellismo-converter/issues)で報告してください。

---

**Happy analyzing! 🧬✨**
