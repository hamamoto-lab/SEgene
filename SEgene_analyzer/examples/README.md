# SEgene Analyzer Usage Examples

This directory contains practical usage examples for SEgene Analyzer.

## Overview

SEgene Analyzer is a comprehensive tool for analyzing super-enhancer regions using SEdb (Super-enhancer Database) data. The examples below demonstrate efficient analysis workflows.

## Available Examples

### 1. Efficient Batch Analysis (`efficient_batch_analysis.sh`)

This script demonstrates how to perform large-scale batch analysis efficiently.

**Features:**
- Pre-calculation and caching of statistics
- Support for both TSV and BED file formats
- Comprehensive workflow including single region analysis
- Optimized parameter settings

## Usage Instructions

### Preparation

1. **Set up data directories**
```bash
# Create data directories
mkdir -p data/SEdb output cache results

# Place required data files:
# - data/SEdb/SE_package_hg38.bed
# - data/SEdb/human_sample_information_sedb2.txt
# - output/20241209_sorted_SE_data.tsv
# - output/regions.bed
```

2. **Prepare region files**

**TSV file example** (`regions.tsv`):
```tsv
chromosome	start	end	name	description
chr1	1000000	2000000	region1	Test region 1
chr2	5000000	6000000	region2	Test region 2
chr7	27000000	28000000	MYC_region	MYC locus region
chr8	128000000	129000000	MYC_enhancer	MYC enhancer region
```

**BED file example** (`regions.bed`):
```bed
chr1	1000000	2000000	region1
chr2	5000000	6000000	region2
chr7	27000000	28000000	MYC_region
chr8	128000000	129000000	MYC_enhancer
```

### Execution

```bash
# Make the script executable
chmod +x efficient_batch_analysis.sh

# Execute
./efficient_batch_analysis.sh
```

## Detailed Usage Examples

### 1. Pre-calculating Statistics

```bash
# Cache statistics for large-scale analysis
sedb-analyzer prepare-stats \
  -b data/SEdb/SE_package_hg38.bed \
  -s data/SEdb/human_sample_information_sedb2.txt \
  -o cache/sedb_stats.pkl \
  --verbose
```

**Effects of this process:**
- Calculate tissue distribution
- Count cell_id occurrences
- Filter for standard chromosomes
- Save results to Pickle file

### 2. Batch Analysis with TSV File

```bash
# Batch analyze multiple regions from TSV file
sedb-analyzer batch \
  -b data/SEdb/SE_package_hg38.bed \
  -s data/SEdb/human_sample_information_sedb2.txt \
  -r regions.tsv \
  -o results/tsv_analysis \
  --use-cached-stats cache/sedb_stats.pkl \
  --max-rows 100 \
  --figure-formats png svg \
  --peak-threshold 0.2 \
  --pvalue-threshold 0.05 \
  --fdr-threshold 0.05 \
  --fc-threshold 1.5 \
  --verbose
```

**Output files:**
- `results/tsv_analysis/data/batch_analysis_results_YYYYMMDD_HHMMSS.tsv`
- `results/tsv_analysis/reports/batch_analysis_summary_YYYYMMDD_HHMMSS.html`
- `results/tsv_analysis/reports/figures/` (various graphs)

### 3. Batch Analysis with BED File

```bash
# Batch analyze multiple regions from BED file
sedb-analyzer batch \
  -b data/SEdb/SE_package_hg38.bed \
  -s data/SEdb/human_sample_information_sedb2.txt \
  -r regions.bed \
  -o results/bed_analysis \
  --use-cached-stats cache/sedb_stats.pkl \
  --figure-formats png pdf \
  --no-tables \
  --verbose
```

### 4. Single Region Detailed Analysis

```bash
# Analyze a specific region in detail
sedb-analyzer single \
  -b data/SEdb/SE_package_hg38.bed \
  -s data/SEdb/human_sample_information_sedb2.txt \
  --chr chr7 \
  --start 5000000 \
  --end 6000000 \
  --region-name "MYC_promoter_region" \
  -o results/single_region \
  --use-cached-stats cache/sedb_stats.pkl \
  --figure-formats png svg pdf \
  --peak-threshold 0.15 \
  --pvalue-threshold 0.01 \
  --verbose
```

**Output files:**
- `results/single_region/chr7_5000000_6000000_analysis_report.html`
- `results/single_region/chr7_5000000_6000000_tissue_enrichment.png`
- `results/single_region/chr7_5000000_6000000_peak_analysis.png`

### 5. Database-wide Report Generation

```bash
# Generate SEdb database metadata report
sedb-analyzer report \
  -b data/SEdb/SE_package_hg38.bed \
  -s data/SEdb/human_sample_information_sedb2.txt \
  -o results/database_report \
  --top-n 25 \
  --dpi 300 \
  --verbose
```

## Advanced Usage Examples

### Custom Threshold Analysis

```bash
# Analysis with more stringent thresholds
sedb-analyzer batch \
  -b data/SEdb/SE_package_hg38.bed \
  -s data/SEdb/human_sample_information_sedb2.txt \
  -r important_regions.tsv \
  -o results/strict_analysis \
  --use-cached-stats cache/sedb_stats.pkl \
  --peak-threshold 0.5 \
  --pvalue-threshold 0.001 \
  --fdr-threshold 0.01 \
  --fc-threshold 2.0 \
  --figure-formats png svg pdf eps \
  --verbose
```

### High-speed Processing (No Peak Analysis)

```bash
# Disable peak analysis for faster processing
sedb-analyzer batch \
  -b data/SEdb/SE_package_hg38.bed \
  -s data/SEdb/human_sample_information_sedb2.txt \
  -r large_regions.tsv \
  -o results/fast_analysis \
  --use-cached-stats cache/sedb_stats.pkl \
  --no-peak-analysis \
  --no-tables \
  --figure-formats png \
  --max-rows 1000 \
  --verbose
```

### Staged Batch Processing

```bash
# Process large datasets in batches
for i in {1..10}; do
  echo "Processing batch $i..."
  sedb-analyzer batch \
    -b data/SEdb/SE_package_hg38.bed \
    -s data/SEdb/human_sample_information_sedb2.txt \
    -r regions_batch_${i}.tsv \
    -o results/batch_${i} \
    --use-cached-stats cache/sedb_stats.pkl \
    --max-rows 100 \
    --figure-formats png \
    --verbose
done
```

## Interpreting Results

### Batch Analysis Results TSV File

The results TSV file contains the following information:

- **Basic information**: Chromosome, start position, end position, region name
- **Overlap information**: Number of overlapping samples, overlap rate
- **Tissue enrichment**: Significant tissue types, P-values, FDR values
- **Peak information**: Number of detected peaks, peak positions
- **Statistical information**: Fold changes, confidence intervals

### HTML Reports

HTML reports include:

- **Summary information**: Analysis overview and key results
- **Tissue distribution plots**: Frequency of occurrence in each tissue
- **Enrichment results**: Statistically significant tissues
- **Peak analysis results**: Peak positions and intensities
- **Detailed charts**: Comprehensive visualization and analysis

## Performance Optimization

### Recommendations

1. **Utilize caching**: Pre-calculate statistics and use cache
2. **Appropriate batch sizes**: Limit processing with `--max-rows`
3. **Essential features only**: Use `--no-peak-analysis` or `--no-tables` for speed
4. **Select output formats**: Specify only necessary `--figure-formats`

### Parallel Processing

```bash
# Run multiple batches in parallel
sedb-analyzer batch -r batch1.tsv -o results/batch1 --use-cached-stats cache/stats.pkl &
sedb-analyzer batch -r batch2.tsv -o results/batch2 --use-cached-stats cache/stats.pkl &
sedb-analyzer batch -r batch3.tsv -o results/batch3 --use-cached-stats cache/stats.pkl &
wait
```

## Troubleshooting

### Common Issues

1. **Out of memory**: Limit processing with `--max-rows`
2. **Long processing time**: Disable peak analysis with `--no-peak-analysis`
3. **File format errors**: Check input file format
4. **Permission errors**: Check output directory write permissions

### Debugging Tips

```bash
# Output detailed logs
sedb-analyzer batch ... --verbose

# Test with small sample
sedb-analyzer batch ... --max-rows 10

# Check statistics cache
ls -la cache/sedb_stats.pkl
```

## Related Documentation

- [Usage Guide](../docs/usage_en.md)
- [Installation Guide](../docs/segene_analyzer_installation.md)
- [API Reference](../README.md)

## Update History

- **2025-01-17**: Initial version created
- Added comprehensive usage examples
- Added performance optimization guidelines
- Added troubleshooting section