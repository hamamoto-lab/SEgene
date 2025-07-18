# SEgene Analyzer 使用ガイド

## 概要

SEgene Analyzer は、SEdb（Super-enhancer Database）データを使用してスーパーエンハンサー領域を解析するための包括的なツールです。本ドキュメントでは、`sedb-analyzer` コマンドラインインターフェースの詳細な使用方法を説明します。

## 基本的な使用方法

### 必要なファイル

1. **SEdb BEDファイル** (`SE_package_hg38.bed`): スーパーエンハンサー領域情報
2. **サンプル情報ファイル** (`human_sample_information_sedb2.txt`): 各サンプルのメタデータ
3. **対象領域ファイル** (TSVまたはBED形式): 解析する遺伝子領域のリスト

### コマンドライン構文

```bash
sedb-analyzer <command> [options]
```

## 主要コマンド

### 1. prepare-stats - 統計情報の事前計算

データベースの統計情報を事前に計算し、キャッシュファイルに保存します。

```bash
sedb-analyzer prepare-stats \
  -b data/SEdb/SE_package_hg38.bed \
  -s data/SEdb/human_sample_information_sedb2.txt \
  -o cache/sedb_stats.pkl
```

**オプション：**
- `-b, --bed-file`: SEdb BEDファイルへのパス
- `-s, --sample-info`: サンプル情報ファイルへのパス
- `-o, --output`: キャッシュファイルの出力パス

### 2. batch - バッチ解析

複数の領域を一度に解析します。

```bash
sedb-analyzer batch \
  -b data/SEdb/SE_package_hg38.bed \
  -s data/SEdb/human_sample_information_sedb2.txt \
  -r regions.tsv \
  -o results/ \
  --use-cached-stats cache/sedb_stats.pkl
```

**オプション：**
- `-r, --regions-file`: 解析する領域のリスト（TSVまたはBED形式）
- `-o, --output-dir`: 結果の出力ディレクトリ
- `--use-cached-stats`: 事前計算済み統計ファイルを使用
- `--max-rows`: 処理する最大行数（テスト用）

### 3. single - 単一領域解析

特定のゲノム領域を詳細に解析します。

```bash
sedb-analyzer single \
  -b data/SEdb/SE_package_hg38.bed \
  -s data/SEdb/human_sample_information_sedb2.txt \
  --chr chr7 --start 1000000 --end 2000000 \
  -o results/
```

**オプション：**
- `--chr`: 染色体名
- `--start`: 開始位置
- `--end`: 終了位置
- `--region-name`: 領域名（オプション）

### 4. report - データベースレポート生成

SEdbデータベースの包括的なメタデータレポートを生成します。

```bash
sedb-analyzer report \
  -b data/SEdb/SE_package_hg38.bed \
  -s data/SEdb/human_sample_information_sedb2.txt \
  -o report/
```

## 高度な使用例

### 1. 効率的なバッチ処理ワークフロー

```bash
# 1. 統計情報の事前計算
sedb-analyzer prepare-stats \
  -b data/SEdb/SE_package_hg38.bed \
  -s data/SEdb/human_sample_information_sedb2.txt \
  -o cache/sedb_stats.pkl

# 2. 高速バッチ解析
sedb-analyzer batch \
  -b data/SEdb/SE_package_hg38.bed \
  -s data/SEdb/human_sample_information_sedb2.txt \
  -r regions.tsv \
  -o results/ \
  --use-cached-stats cache/sedb_stats.pkl \
  --figure-formats png svg

# 3. 特定領域の詳細解析
sedb-analyzer single \
  -b data/SEdb/SE_package_hg38.bed \
  -s data/SEdb/human_sample_information_sedb2.txt \
  --chr chr7 --start 1000000 --end 2000000 \
  -o results/chr7_detailed/ \
  --use-cached-stats cache/sedb_stats.pkl
```

### 2. カスタム閾値を使用した解析

```bash
sedb-analyzer batch \
  -r regions.tsv \
  -b data/SEdb/SE_package_hg38.bed \
  -s data/SEdb/human_sample_information_sedb2.txt \
  -o results/ \
  --pvalue-threshold 0.01 \
  --fdr-threshold 0.05 \
  --fc-threshold 2.0
```

### 3. 高速処理のための並列実行

```bash
# 領域ファイルを分割
split -l 100 regions.tsv region_batch_

# 並列実行
for batch in region_batch_*; do
  sedb-analyzer batch \
    -r "$batch" \
    -b data/SEdb/SE_package_hg38.bed \
    -s data/SEdb/human_sample_information_sedb2.txt \
    -o "results_${batch##*_}/" \
    --use-cached-stats cache/sedb_stats.pkl &
done
wait
```

## 出力形式

### バッチ解析結果

結果は以下の形式で出力されます：

```
results/
├── batch_analysis_results_20250117_123456.tsv  # メイン結果ファイル
├── data/                                      # 中間データ
├── reports/                                   # 個別レポート
│   ├── region_1_chr1_1000000_2000000/
│   │   ├── enrichment_analysis.tsv
│   │   ├── tissue_enrichment_plot.png
│   │   └── summary_report.html
│   └── region_2_chr2_3000000_4000000/
└── sedb_report/                               # データベースレポート
    ├── database_overview.html
    └── figures/
```

### 結果TSVファイルの内容

メイン結果ファイルには以下の情報が含まれます：

- **region_id**: 領域ID
- **se_unique_sample_count**: 関連するユニークサンプル数
- **linked_sample_count**: メタデータにリンクされたサンプル数
- **significant_tissues**: 有意な組織タイプ
- **significant_tissues_fc**: 対応するFold Change
- **significant_tissues_pvalue**: 対応するP値
- **significant_tissues_fdr**: 対応するFDR調整P値

### HTMLレポート

各領域について以下を含む詳細HTMLレポートが生成されます：

- **組織エンリッチメント解析**: 統計的有意性評価
- **可視化グラフ**: インタラクティブなチャート
- **データテーブル**: 詳細な数値データ
- **サマリー統計**: 重要な結果のハイライト

## トラブルシューティング

### 一般的な問題

1. **ファイルが見つからない**
   ```bash
   # ファイルパスを確認
   ls -la data/SEdb/
   ```

2. **メモリ不足**
   ```bash
   # 少数の領域でテスト
   sedb-analyzer batch --max-rows 10 ...
   ```

3. **処理時間が長い**
   ```bash
   # 統計キャッシュを使用
   sedb-analyzer prepare-stats ...
   sedb-analyzer batch --use-cached-stats ...
   ```

### デバッグオプション

```bash
# 詳細なログ出力
sedb-analyzer batch --verbose ...

# ログファイルに出力
sedb-analyzer batch --log-file analysis.log ...
```

## パフォーマンス最適化

### 推奨設定

1. **システム要件**:
   - Python 3.10以上
   - 16GB以上のメモリ
   - マルチコアCPU

2. **効率的な処理**:
   - 統計キャッシュの使用
   - 適切なバッチサイズ
   - 並列処理の活用

3. **最適化されたコマンド例**:
   ```bash
   # 大規模データセット向け
   sedb-analyzer batch \
     -r regions.tsv \
     -b data/SEdb/SE_package_hg38.bed \
     -s data/SEdb/human_sample_information_sedb2.txt \
     -o results/ \
     --use-cached-stats cache/sedb_stats.pkl \
     --figure-formats png \
     --no-tables
   ```

## 参考資料

- [インストールガイド](segene_analyzer_installation_ja.md)
- [使用例](../examples/README_ja.md)
- [プロジェクト概要](../README_ja.md)

