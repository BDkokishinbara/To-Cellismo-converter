# シングルセル解析ファイル変換Webアプリケーション

CSV、RDS、MEXファイルをh5mu形式に変換するWebアプリケーションです。

## 概要

このツールは、シングルセル解析で使用される様々なファイル形式をMuData形式（h5mu）に変換します。h5muファイルは、マルチモーダルなシングルセルデータを効率的に保存でき、PythonやRで広く利用可能です。

### 対応ファイル形式

- **CSV**: 遺伝子×細胞または細胞×遺伝子のマトリックスファイル
- **RDS**: SeuratオブジェクトまたはSingleCellExperimentオブジェクト
- **MEX**: 10x Genomics形式の疎行列ファイル（ZIP形式）

## インストール

### 必要要件

- Python 3.8以上
- R（RDSファイルの変換に必要）

### Python環境のセットアップ

1. リポジトリのクローンまたはダウンロード

```bash
cd csv_h5mu
```

2. 仮想環境の作成（推奨）

```bash
python -m venv venv
source venv/bin/activate  # Linux/Mac
# または
venv\Scripts\activate  # Windows
```

3. 依存パッケージのインストール

```bash
pip install -r requirements.txt
```

### RDSファイルのサポート（オプション）

RDSファイルの変換には、Rとrpy2が必要です。

#### Rのインストール

**Ubuntu/Debian:**
```bash
sudo apt-get update
sudo apt-get install r-base r-base-dev
```

**macOS (Homebrew):**
```bash
brew install r
```

**Windows:**
[CRAN](https://cran.r-project.org/)からRをダウンロードしてインストール

#### 必要なRパッケージのインストール

Rを起動して以下を実行：

```r
install.packages("Seurat")
install.packages("SingleCellExperiment")
install.packages("BiocManager")
BiocManager::install("SingleCellExperiment")
```

## 使い方

### ローカル開発環境での起動

```bash
python app.py
```

アプリケーションは `http://localhost:5000` で起動します。

ブラウザで上記のURLにアクセスしてください。

### 基本的な使い方

1. **ファイル形式を選択**: CSV、RDS、MEXから選択
2. **ファイルをアップロード**: クリックまたはドラッグ&ドロップ
3. **変換オプションを設定**（CSVの場合）
4. **「h5muに変換」ボタンをクリック**
5. **変換完了後、自動的にダウンロードが開始**

## デプロイ

本番環境へのデプロイ方法については [DEPLOYMENT.md](DEPLOYMENT.md) をご覧ください。

### クイックデプロイオプション

#### Docker を使用
```bash
docker-compose up -d
```

#### VPS/Ubuntu サーバー
```bash
chmod +x deploy.sh
./deploy.sh
```

#### クラウドプラットフォーム
- **Render**: GitHubリポジトリを接続するだけで自動デプロイ
- **Railway**: ワンクリックデプロイ
- **Fly.io**: `flyctl launch` でデプロイ

詳細な手順は [DEPLOYMENT.md](DEPLOYMENT.md) を参照してください。

### 各ファイル形式の詳細

#### CSV形式

遺伝子発現マトリックスをCSV形式で保存したファイル。

**フォーマット例（遺伝子×細胞）:**
```csv
,Cell1,Cell2,Cell3
Gene1,5,10,3
Gene2,2,8,15
Gene3,0,5,7
```

**オプション:**
- **マトリックスを転置**: 遺伝子×細胞形式の場合にチェック
- **1行目をヘッダーとして扱う**: 細胞名が含まれる場合
- **1列目をインデックスとして扱う**: 遺伝子名が含まれる場合

#### RDS形式

RのSeuratオブジェクトまたはSingleCellExperimentオブジェクトを保存したファイル。

**サポートされるオブジェクト:**
- Seurat v3/v4/v5
- SingleCellExperiment

**注意:** RDSファイルの変換には、RとRパッケージのインストールが必要です。

#### MEX形式

10x Genomicsなどで使用される疎行列形式。

**必要なファイル（ZIP圧縮）:**
- `matrix.mtx` または `matrix.mtx.gz`: 疎行列データ
- `features.tsv` または `genes.tsv`: 遺伝子情報
- `barcodes.tsv`: 細胞バーコード情報

**準備方法:**
```bash
# MEXファイルをZIP圧縮
zip mex_data.zip matrix.mtx.gz features.tsv.gz barcodes.tsv.gz
```

または、これらのファイルを含むディレクトリをZIP圧縮：
```bash
zip -r mex_data.zip filtered_feature_bc_matrix/
```

## h5muファイルの利用

### Python

```python
import mudata as md
import scanpy as sc

# h5muファイルを読み込み
mdata = md.read_h5mu('converted.h5mu')

# RNAデータにアクセス
adata = mdata['rna']

# Scanpyで解析
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
```

### R

```r
library(MuData)

# h5muファイルを読み込み
mdata <- readH5MU('converted.h5mu')

# RNAデータにアクセス
rna_data <- mdata[['rna']]
```

## プロジェクト構造

```
csv_h5mu/
├── app.py                    # Flaskアプリケーション
├── converters/               # 変換モジュール
│   ├── __init__.py
│   ├── csv_converter.py      # CSV変換
│   ├── rds_converter.py      # RDS変換
│   └── mex_converter.py      # MEX変換
├── templates/                # HTMLテンプレート
│   └── index.html
├── static/                   # 静的ファイル
│   ├── style.css
│   └── script.js
├── uploads/                  # アップロード一時ディレクトリ
├── outputs/                  # 変換結果出力ディレクトリ
├── requirements.txt          # Python依存関係
└── README.md                 # このファイル
```

## トラブルシューティング

### rpy2のインストールエラー

**エラー:** `ERROR: Failed building wheel for rpy2`

**解決策:**
1. Rがインストールされていることを確認
2. 開発ツールをインストール:
   ```bash
   # Ubuntu/Debian
   sudo apt-get install build-essential libcurl4-openssl-dev libssl-dev

   # macOS
   xcode-select --install
   ```

### メモリエラー

大きなファイルを変換する際にメモリ不足エラーが発生する場合：

1. より少ないセルでテスト
2. システムのメモリを増やす
3. スパース行列形式が適切に使用されているか確認

### ZIPファイルエラー

MEXファイルのZIPアップロードでエラーが発生する場合：

1. ZIP内にmatrix.mtx、features.tsv、barcodes.tsvが含まれているか確認
2. ファイル名が正しいか確認（大文字小文字も含む）
3. gzip圧縮されたファイル（.gz）も対応しています

## 開発者向け情報

### 新しい変換形式の追加

1. `converters/`ディレクトリに新しいコンバーターモジュールを作成
2. `app.py`にインポートして統合
3. `templates/index.html`にUIオプションを追加

### テスト

```bash
# サンプルCSVファイルで基本テスト
python -c "
import pandas as pd
import numpy as np

# サンプルデータ作成
data = pd.DataFrame(
    np.random.randint(0, 100, size=(100, 50)),
    index=[f'Gene{i}' for i in range(100)],
    columns=[f'Cell{i}' for i in range(50)]
)
data.to_csv('sample.csv')
"

# アプリケーション起動
python app.py
```

## ライセンス

このプロジェクトはMITライセンスの下で公開されています。

## 貢献

バグ報告、機能リクエスト、プルリクエストを歓迎します。

## 参考資料

- [MuData Documentation](https://mudata.readthedocs.io/)
- [Scanpy Documentation](https://scanpy.readthedocs.io/)
- [Seurat Documentation](https://satijalab.org/seurat/)
- [10x Genomics MEX Format](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices)

## 謝辞

このツールは以下のオープンソースプロジェクトを使用しています：

- Flask
- AnnData
- MuData
- Scanpy
- pandas
- NumPy
- SciPy
- rpy2
