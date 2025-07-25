# eRNAbaseメタデータ

SEgene_analyzer_ernaの解析に必要なeRNAbaseメタデータファイルです。

**データソース**: [eRNAbase](https://bio.liclab.net/eRNAbase/index.php) - エンハンサーRNAデータベース  
**データバージョン**: このメタデータは **2025年7月22日時点** のeRNAbaseデータに対応しています。

## ファイル

- **eRNAbase_data.csv**: CSV形式
- **eRNAbase_data.parquet**: Parquet形式

## データ統計

- **総サンプル数**: 1,012
- **ヒトサンプル**: 858、**マウスサンプル**: 154
- **組織タイプ**: 43種類、**細胞タイプ**: 119種類

## 使用方法

```bash
# メタデータファイルを指定して解析実行
erna-analyzer single -m metadata/eRNAbase_data.parquet --chr chr7 --start 1000000 --end 2000000
```

## データ構造

| カラム | 説明 |
|--------|------|
| Sample ID | サンプル識別子 |
| Species | 種名（Homo sapiens / Mus musculus） |
| Tissue Type | 組織タイプ |
| Cell Type | 細胞タイプ |
| Biosample Type | サンプルタイプ |
| Experiment Type | 実験タイプ |