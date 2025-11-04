---
title: To-Cellismo Converter
emoji: 🧬
colorFrom: blue
colorTo: purple
sdk: docker
pinned: false
license: mit
---

# 🧬 To-Cellismo Converter

シングルセル解析ファイル変換ツール - CSV/RDS/MEX から h5mu 形式への変換

## 📋 概要

このアプリケーションは、シングルセル解析で使用される様々なファイル形式を MuData 形式（h5mu）に変換します。

## ✨ 対応ファイル形式

- **CSV**: 遺伝子×細胞または細胞×遺伝子のマトリックスファイル
- **RDS**: Seurat オブジェクトまたは SingleCellExperiment オブジェクト
- **MEX**: 10x Genomics 形式の疎行列ファイル（ZIP形式）

## 🚀 使い方

1. ファイル形式を選択（CSV/RDS/MEX）
2. ファイルをアップロード（ドラッグ&ドロップまたはクリック）
3. 変換オプションを設定（CSV の場合）
4. 「h5mu に変換」ボタンをクリック
5. 変換完了後、自動的にダウンロードが開始されます

## 📊 h5mu 形式について

h5mu は MuData フォーマットで、マルチモーダルなシングルセル解析データを効率的に保存できます。Python（mudata、scanpy）や R（MuData）で利用可能です。

## 💡 使用例

### Python
```python
import mudata as md
import scanpy as sc

# h5mu ファイルを読み込み
mdata = md.read_h5mu('converted.h5mu')

# RNA データにアクセス
adata = mdata['rna']

# Scanpy で解析
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
```

### R
```r
library(MuData)

# h5mu ファイルを読み込み
mdata <- readH5MU('converted.h5mu')
```

## 📝 ライセンス

MIT License - 商用利用を含め、自由に使用できます。

## 🔗 リンク

- [GitHub リポジトリ](https://github.com/BDkokishinbara/To-Cellismo-converter)
- [MuData Documentation](https://mudata.readthedocs.io/)
- [Scanpy Documentation](https://scanpy.readthedocs.io/)

## 🙏 謝辞

このツールは以下のオープンソースプロジェクトを使用しています：
Flask, AnnData, MuData, Scanpy, pandas, NumPy, SciPy
