# データディレクトリ

このディレクトリは、SEdbデータファイルやその他のデータセットを保存するためのものです。

## SEgene_analyzerに必要なファイル

[SEdb 2.0](http://www.licpathway.net/sedb/download.php) から以下のファイルをダウンロードし、`SEdb/` サブディレクトリに配置してください：

### SEdb 2.0 ファイル
- `human_sample_information_sedb2.txt` - サンプルメタデータ
- `SE_package_hg38.bed` - ヒト（hg38）のスーパーエンハンサー定義

### ディレクトリ構造
```
data/
└── SEdb/
    ├── human_sample_information_sedb2.txt
    └── SE_package_hg38.bed
```

## ダウンロード手順

1. アクセス: [SEdb 2.0 ダウンロードページ](http://www.licpathway.net/sedb/download.php)
2. ダウンロード: "Sample information" → `human_sample_information_sedb2.txt`
3. ダウンロード: "Package of SE, SE element and TE" → Human (hg38) → `SE_package_hg38.bed`
4. ファイルを `SEdb/` サブディレクトリに配置

## 注意事項

- これらのデータファイルは、サイズとライセンスの考慮により、**gitリポジトリには含まれていません**
- `data/` ディレクトリは `.gitignore` に記載されています
- 各ユーザーは自分でデータファイルのコピーをダウンロードする必要があります