#!/bin/bash

# 効率的なバッチ解析の例

# データファイルのパス
BED_FILE="data/SEdb/SE_package_hg38.bed"
SAMPLE_INFO="data/SEdb/human_sample_information_sedb2.txt"
REGIONS_TSV="output/20241209_sorted_SE_data.tsv"
REGIONS_BED="output/regions.bed"
CACHE_DIR="cache"
OUTPUT_DIR="results"

# 1. 統計情報を事前計算してキャッシュ
echo "Step 1: Preparing statistics cache..."
sedb-analyzer prepare-stats \
  -b "$BED_FILE" \
  -s "$SAMPLE_INFO" \
  -o "$CACHE_DIR/sedb_stats.pkl" \
  --verbose

# 2. TSVファイルのバッチ解析（キャッシュ使用）
echo "Step 2: Batch analysis with TSV file..."
sedb-analyzer batch \
  -b "$BED_FILE" \
  -s "$SAMPLE_INFO" \
  -r "$REGIONS_TSV" \
  -o "$OUTPUT_DIR/tsv_analysis" \
  --use-cached-stats "$CACHE_DIR/sedb_stats.pkl" \
  --max-rows 100 \
  --figure-formats png svg \
  --verbose

# 3. BEDファイルのバッチ解析（同じキャッシュを再利用）
echo "Step 3: Batch analysis with BED file..."
sedb-analyzer batch \
  -b "$BED_FILE" \
  -s "$SAMPLE_INFO" \
  -r "$REGIONS_BED" \
  -o "$OUTPUT_DIR/bed_analysis" \
  --use-cached-stats "$CACHE_DIR/sedb_stats.pkl" \
  --figure-formats png pdf \
  --verbose

# 4. 単一領域の解析（キャッシュ使用）
echo "Step 4: Single region analysis..."
sedb-analyzer single \
  -b "$BED_FILE" \
  -s "$SAMPLE_INFO" \
  --chr chr7 \
  --start 5000000 \
  --end 6000000 \
  -o "$OUTPUT_DIR/single_region" \
  --use-cached-stats "$CACHE_DIR/sedb_stats.pkl" \
  --verbose