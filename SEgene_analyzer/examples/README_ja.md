# SEgene Analyzer 使用例

このディレクトリには、SEgene Analyzer の実用的な使用例が含まれています。

## 概要

SEgene Analyzer は、SEdb（Super-enhancer Database）データを使用してスーパーエンハンサー領域を解析するための包括的なツールです。以下の例では、効率的な解析ワークフローを実演します。

## 利用可能な例

### 1. 効率的なバッチ解析 (`efficient_batch_analysis.sh`)

このスクリプトは、大規模なバッチ解析を効率的に実行する方法を示しています。

**特徴:**
- 統計情報の事前計算とキャッシュ
- TSVとBEDファイル形式の両方をサポート
- 単一領域解析を含む包括的なワークフロー
- 最適化されたパラメータ設定

### 2. 基本的な使用例

#### データ準備

```bash
# 必要なディレクトリを作成
mkdir -p data/SEdb cache results

# SEdb 2.0 データファイルをダウンロード
# - human_sample_information_sedb2.txt
# - SE_package_hg38.bed
# を data/SEdb/ ディレクトリに配置
```

#### 統計情報の事前計算

```bash
# データベース統計を事前計算（推奨）
sedb-analyzer prepare-stats \
  -b data/SEdb/SE_package_hg38.bed \
  -s data/SEdb/human_sample_information_sedb2.txt \
  -o cache/sedb_stats.pkl
```

## 詳細な使用例

### 1. 基本的なバッチ解析

```bash
# 領域ファイル（regions.tsv）を準備
# フォーマット: region_id<TAB>chromosome<TAB>start<TAB>end
echo -e "region1\tchr1\t1000000\t2000000" > regions.tsv
echo -e "region2\tchr2\t3000000\t4000000" >> regions.tsv

# バッチ解析の実行
sedb-analyzer batch \
  -b data/SEdb/SE_package_hg38.bed \
  -s data/SEdb/human_sample_information_sedb2.txt \
  -r regions.tsv \
  -o results/ \
  --use-cached-stats cache/sedb_stats.pkl \
  --figure-formats png svg pdf
```

### 2. カスタム閾値を使用した解析

```bash
# より厳密な統計的閾値を使用
sedb-analyzer batch \
  -b data/SEdb/SE_package_hg38.bed \
  -s data/SEdb/human_sample_information_sedb2.txt \
  -r regions.tsv \
  -o results_strict/ \
  --pvalue-threshold 0.001 \
  --fdr-threshold 0.01 \
  --fc-threshold 3.0 \
  --use-cached-stats cache/sedb_stats.pkl
```

### 3. 高速処理のための最適化

```bash
# 図表生成を最小限に抑えた高速解析
sedb-analyzer batch \
  -b data/SEdb/SE_package_hg38.bed \
  -s data/SEdb/human_sample_information_sedb2.txt \
  -r regions.tsv \
  -o results_fast/ \
  --use-cached-stats cache/sedb_stats.pkl \
  --figure-formats png \
  --no-tables
```

### 4. 単一領域の詳細解析

```bash
# 特定の領域を詳細に解析
sedb-analyzer single \
  -b data/SEdb/SE_package_hg38.bed \
  -s data/SEdb/human_sample_information_sedb2.txt \
  --chr chr7 \
  --start 1000000 \
  --end 2000000 \
  --region-name "MYC_enhancer_region" \
  -o results/single_region/ \
  --use-cached-stats cache/sedb_stats.pkl
```

### 5. データベースレポートの生成

```bash
# SEdbデータベースの包括的なメタデータレポートを生成
sedb-analyzer report \
  -b data/SEdb/SE_package_hg38.bed \
  -s data/SEdb/human_sample_information_sedb2.txt \
  -o report/ \
  --top-n 50 \
  --dpi 300
```

## 大規模データセットでの実践的なワークフロー

### 並列処理による高速化

```bash
# 大きな領域ファイルを複数の小さなファイルに分割
split -l 50 large_regions.tsv region_batch_

# 並列実行
for batch in region_batch_*; do
  sedb-analyzer batch \
    -b data/SEdb/SE_package_hg38.bed \
    -s data/SEdb/human_sample_information_sedb2.txt \
    -r "$batch" \
    -o "results_${batch##*_}/" \
    --use-cached-stats cache/sedb_stats.pkl &
done

# 全ての並列処理の完了を待つ
wait

# 結果を統合
cat results_*/batch_analysis_results_*.tsv > combined_results.tsv
```

### メモリ効率的な処理

```bash
# 大規模データセットでのメモリ効率化
sedb-analyzer batch \
  -b data/SEdb/SE_package_hg38.bed \
  -s data/SEdb/human_sample_information_sedb2.txt \
  -r regions.tsv \
  -o results/ \
  --use-cached-stats cache/sedb_stats.pkl \
  --max-rows 100 \
  --no-peak-analysis
```

## 出力ファイルの解釈

### 結果TSVファイル

バッチ解析結果のTSVファイルには以下の情報が含まれます：

- **region_id**: 領域の識別子
- **chromosome, start, end**: ゲノム座標
- **se_unique_sample_count**: 重複を除いたサンプル数
- **linked_sample_count**: メタデータにリンクされたサンプル数
- **significant_tissues**: 統計的に有意な組織タイプ
- **significant_tissues_fc**: 対応するFold Change値
- **significant_tissues_pvalue**: 対応するP値
- **significant_tissues_fdr**: FDR調整後のP値

### HTMLレポート

各領域について以下を含む詳細なHTMLレポートが生成されます：

- **組織エンリッチメント解析**: 統計的有意性の評価
- **可視化**: 詳細なチャートとグラフ
- **データテーブル**: 詳細な数値データ
- **サマリー統計**: 重要な結果のハイライト

## よくある問題とその対処法

### 1. メモリ不足

**症状**: 処理中にメモリエラーが発生する

**対処法**:
```bash
# 少数の領域でテスト
sedb-analyzer batch --max-rows 10 ...

# 統計キャッシュを使用
sedb-analyzer prepare-stats ...
sedb-analyzer batch --use-cached-stats ...
```

### 2. 処理時間が長い

**症状**: 解析に予想以上の時間がかかる

**対処法**:
```bash
# 統計キャッシュを使用
sedb-analyzer prepare-stats ...

# 出力を最小限に抑える
sedb-analyzer batch --figure-formats png --no-tables ...
```

### 3. ファイルが見つからない

**症状**: データファイルが見つからないエラー

**対処法**:
```bash
# ファイルパスを確認
ls -la data/SEdb/
realpath data/SEdb/SE_package_hg38.bed
```

## パフォーマンス最適化のヒント

### 1. 統計キャッシュの活用

```bash
# 一度だけ実行
sedb-analyzer prepare-stats -b data/SEdb/SE_package_hg38.bed -s data/SEdb/human_sample_information_sedb2.txt -o cache/sedb_stats.pkl

# 以後の解析で使用
sedb-analyzer batch --use-cached-stats cache/sedb_stats.pkl ...
```

### 2. 適切なバッチサイズ

```bash
# 小さなバッチサイズで高速処理
sedb-analyzer batch --max-rows 50 ...

# 大きなバッチサイズで効率的な処理
sedb-analyzer batch --max-rows 500 ...
```

### 3. 出力形式の最適化

```bash
# 高速処理用（最小限の出力）
sedb-analyzer batch --figure-formats png --no-tables ...

# 詳細解析用（完全な出力）
sedb-analyzer batch --figure-formats png svg pdf eps ...
```

## 高度な使用例

### カスタム解析パイプライン

```bash
#!/bin/bash

# 1. データベース統計の事前計算
echo "統計情報を事前計算中..."
sedb-analyzer prepare-stats \
  -b data/SEdb/SE_package_hg38.bed \
  -s data/SEdb/human_sample_information_sedb2.txt \
  -o cache/sedb_stats.pkl

# 2. 予備的な解析（少数の領域）
echo "予備解析を実行中..."
sedb-analyzer batch \
  -b data/SEdb/SE_package_hg38.bed \
  -s data/SEdb/human_sample_information_sedb2.txt \
  -r regions.tsv \
  -o results_preliminary/ \
  --use-cached-stats cache/sedb_stats.pkl \
  --max-rows 10

# 3. 完全な解析
echo "完全な解析を実行中..."
sedb-analyzer batch \
  -b data/SEdb/SE_package_hg38.bed \
  -s data/SEdb/human_sample_information_sedb2.txt \
  -r regions.tsv \
  -o results_complete/ \
  --use-cached-stats cache/sedb_stats.pkl \
  --figure-formats png svg

# 4. データベースレポートの生成
echo "データベースレポートを生成中..."
sedb-analyzer report \
  -b data/SEdb/SE_package_hg38.bed \
  -s data/SEdb/human_sample_information_sedb2.txt \
  -o report/

echo "解析完了！"
```

### 品質管理チェック

```bash
# 結果の品質管理
echo "解析結果の品質チェック:"
echo "- 総領域数: $(wc -l < regions.tsv)"
echo "- 処理完了領域数: $(grep -c "region_" results_complete/batch_analysis_results_*.tsv)"
echo "- 有意な結果を持つ領域数: $(grep -c -v "^$" results_complete/batch_analysis_results_*.tsv | grep significant_tissues)"
```

## 関連リソース

- [使用ガイド](../docs/usage_ja.md) - 詳細なコマンドリファレンス
- [インストールガイド](../docs/segene_analyzer_installation_ja.md) - セットアップ手順
- [メインドキュメント](../README_ja.md) - プロジェクト概要

