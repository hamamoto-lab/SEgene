#!/bin/bash
#
# eRNAbase Region Analyzer CLI Usage Examples
#
# This script demonstrates various CLI commands for the erna-analyzer tool.
# Make sure to install the package first:
#   cd /path/to/SEgene_analyzer_erna
#   pip install -e .
#

# Set data paths (adjust to your environment)
BED_DIR="data/peak"
METADATA="metadata/eRNAbase_data.parquet"
REGIONS_TSV="data/20241209_sorted_SE_data.tsv"
OUTPUT_BASE="results/cli_test"

echo "=== eRNAbase Region Analyzer CLI Examples ==="
echo

# 1. Show help
echo "1. Display help information:"
echo "erna-analyzer --help"
echo

# 2. Prepare statistics cache
echo "2. Pre-calculate and cache database statistics:"
echo "erna-analyzer prepare-stats \\"
echo "  -b $BED_DIR \\"
echo "  -m $METADATA \\"
echo "  -o $OUTPUT_BASE/cache/erna_stats.pkl \\"
echo "  --species human"
echo

# 3. Single region analysis
echo "3. Analyze a single genomic region:"
echo "erna-analyzer single \\"
echo "  -b $BED_DIR \\"
echo "  -m $METADATA \\"
echo "  --chr chr7 \\"
echo "  --start 1000000 \\"
echo "  --end 2000000 \\"
echo "  --region-name 'Test_Region_Chr7' \\"
echo "  -o $OUTPUT_BASE/single \\"
echo "  --species human \\"
echo "  --figure-formats png svg \\"
echo "  --pvalue-threshold 0.05"
echo

# 4. Batch analysis (without cache)
echo "4. Batch analysis without cache:"
echo "erna-analyzer batch \\"
echo "  -b $BED_DIR \\"
echo "  -m $METADATA \\"
echo "  -r $REGIONS_TSV \\"
echo "  -o $OUTPUT_BASE/batch_no_cache \\"
echo "  --species human \\"
echo "  --max-rows 5 \\"
echo "  --figure-formats png svg \\"
echo "  --fc-threshold 1.5"
echo

# 5. Batch analysis (with cache)
echo "5. Batch analysis with cached statistics:"
echo "erna-analyzer batch \\"
echo "  -b $BED_DIR \\"
echo "  -m $METADATA \\"
echo "  -r $REGIONS_TSV \\"
echo "  -o $OUTPUT_BASE/batch_with_cache \\"
echo "  --use-cached-stats $OUTPUT_BASE/cache/erna_stats.pkl \\"
echo "  --species human \\"
echo "  --max-rows 10 \\"
echo "  --figure-formats png"
echo

# 6. Generate metadata report
echo "6. Generate eRNAbase metadata report:"
echo "erna-analyzer report \\"
echo "  -m $METADATA \\"
echo "  -o $OUTPUT_BASE/report \\"
echo "  --species human \\"
echo "  --top-n 20 \\"
echo "  --figure-formats png svg pdf"
echo

# 7. Verbose mode with logging
echo "7. Run with verbose logging:"
echo "erna-analyzer single \\"
echo "  -b $BED_DIR \\"
echo "  -m $METADATA \\"
echo "  --chr chr1 \\"
echo "  --start 1000000 \\"
echo "  --end 2000000 \\"
echo "  -o $OUTPUT_BASE/verbose \\"
echo "  -v \\"
echo "  --log-file $OUTPUT_BASE/analysis.log"
echo

# 8. Save statistics cache for future use
echo "8. Batch analysis and save statistics cache:"
echo "erna-analyzer batch \\"
echo "  -b $BED_DIR \\"
echo "  -m $METADATA \\"
echo "  -r $REGIONS_TSV \\"
echo "  -o $OUTPUT_BASE/batch_save_cache \\"
echo "  --save-stats-cache $OUTPUT_BASE/cache/new_stats.pkl \\"
echo "  --species human \\"
echo "  --max-rows 3"
echo

echo "=== Notes ==="
echo "- Adjust the data paths according to your environment"
echo "- Use --max-rows to limit the number of regions for testing"
echo "- Cached statistics significantly speed up repeated analyses"
echo "- All commands support -v flag for verbose output"
echo "- Output directories are created automatically"