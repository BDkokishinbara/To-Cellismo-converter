# 🧬 To-Cellismo Converter

シングルセル解析ファイルを簡単に変換！CSV/RDS/MEX から h5mu 形式への変換ツール

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

## 🚀 クイックスタート

初めて使う方向けの簡単3ステップ！

### ステップ1: ダウンロード
```bash
git clone https://github.com/BDkokishinbara/To-Cellismo-converter.git
cd To-Cellismo-converter
```

### ステップ2: インストール
```bash
pip install -r requirements.txt
```

### ステップ3: 起動

#### 🖱️ 簡単起動（推奨）

**Linux / macOS**:
```bash
./start.sh
```
または、ファイルマネージャーで `start.sh` をダブルクリック

**Windows**:
```
start.bat
```
または、エクスプローラーで `start.bat` をダブルクリック

起動スクリプトを使うと、自動的にブラウザが開きます！

#### ⌨️ 手動起動

```bash
python app.py
```

ブラウザで **http://127.0.0.1:5000** にアクセスして使い始めましょう！

停止するには、ターミナルで `Ctrl+C` を押してください。

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

- **Python 3.8以上**
- **インターネット接続**（初回インストール時）

### ステップ1: リポジトリをダウンロード

```bash
# GitHubからダウンロード
git clone https://github.com/BDkokishinbara/To-Cellismo-converter.git
cd To-Cellismo-converter
```

または、GitHubのページから「Code」→「Download ZIP」でダウンロードして解凍してください。

### ステップ2: 必要なパッケージをインストール

```bash
# 必要なPythonパッケージをインストール
pip install -r requirements.txt
```

### ステップ3（オプション）: RDSファイル変換を使う場合

RDSファイルを変換したい場合は、追加でrpy2をインストールしてください：

```bash
# Rのインストール（Ubuntuの場合）
sudo apt-get install r-base r-base-dev

# rpy2のインストール
pip install rpy2
```

**注意**: RDSファイルを使わない場合、この手順はスキップできます。

---

## 使い方

### 🌐 Webアプリを使う（推奨）

#### 1. アプリを起動

```bash
# プロジェクトのフォルダに移動
cd To-Cellismo-converter

# Flaskアプリを起動
python app.py
```

起動すると以下のようなメッセージが表示されます：

```
Starting Single-cell File Converter Web Application...
Access the application at: http://0.0.0.0:5000
```

#### 2. ブラウザでアクセス

ブラウザで以下のURLを開きます：

**http://127.0.0.1:5000**

#### 3. ファイルを変換

1. **ファイル形式を選択**: CSV / RDS / MEX から選ぶ
2. **ファイルをアップロード**: ドラッグ&ドロップまたはクリックして選択
3. **オプションを設定**（CSVの場合）:
   - **Transpose**: データが遺伝子×細胞の場合にチェック
   - **Has header**: 1行目がヘッダー（細胞名）の場合にチェック（通常はオン）
   - **Has index**: 1列目がインデックス（遺伝子名）の場合にチェック（通常はオン）
4. **「h5mu に変換」ボタンをクリック**
5. **自動的にダウンロード開始**: `convert_元のファイル名.h5mu` という名前で保存されます

#### 4. アプリを停止

ターミナルで `Ctrl+C` を押すと停止します。

---

### 💻 コマンドラインで使う

#### CSVファイルの変換

```python
from converters.csv_converter import csv_to_h5mu

# 基本的な使い方
csv_to_h5mu(
    csv_path='your_data.csv',
    output_path='output.h5mu',
    transpose=False,      # 遺伝子×細胞の場合はTrue
    has_header=True,      # ヘッダー行がある場合
    has_index=True        # インデックス列がある場合
)
```

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

#### エラー1: `ModuleNotFoundError: No module named 'flask'`

**原因**: 必要なパッケージがインストールされていない

**解決方法**:
```bash
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

#### エラー5: `Port 5000 already in use`

**原因**: ポート5000が既に使用されている

**解決方法**:
```bash
# 別のポートを使用する
export PORT=8080
python app.py
```

---

## よくある質問

### Q1: どのくらいのファイルサイズまで対応していますか？

**A**: デフォルトでは500MBまでです。`app.py`の`MAX_CONTENT_LENGTH`を変更することで調整できます。

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

**A**: 現在は1つずつの変換のみ対応しています。複数のファイルを処理したい場合は、Pythonスクリプトでループ処理を書くことができます。

### Q7: オフラインで使えますか？

**A**: はい。一度必要なパッケージをインストールすれば、インターネット接続なしで使用できます。

### Q8: 商用利用できますか？

**A**: はい。MITライセンスの下で、商用利用を含めて自由に使用できます。

---

## 📝 ファイル構成

```
To-Cellismo-converter/
├── app.py                               # Flaskアプリケーション本体
├── requirements.txt                     # 必要なPythonパッケージ
├── README.md                            # このファイル
├── start.sh                             # Linux/macOS用起動スクリプト
├── start.bat                            # Windows用起動スクリプト
├── To-Cellismo-Converter.desktop        # Linuxデスクトップアプリ設定
├── converters/                          # 変換ロジック
│   ├── csv_converter.py                # CSV変換（SeqGeq対応）
│   ├── rds_converter.py                # RDS変換
│   └── mex_converter.py                # MEX変換
├── templates/                           # HTMLテンプレート
│   └── index.html
├── static/                              # CSS/JavaScript
│   ├── style.css
│   └── script.js
├── demodata/                            # デモファイル
│   ├── sample_data.csv                 # サンプルCSV
│   └── GSM4630028_ccRCC1.tar.gz        # サンプルMEX (tar.gz)
├── uploads/                             # アップロードファイル一時保存
└── outputs/                             # 変換後ファイル保存
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
