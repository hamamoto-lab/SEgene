# SEgene Analyzer Usage Guide

## Overview

SEgene Analyzer is a comprehensive tool for analyzing super-enhancer regions using SEdb (Super-enhancer Database) data. This document provides detailed instructions on how to use the `sedb-analyzer` command-line interface.

## Basic Usage

### Required Files

1. **SEdb BED file** (`SE_package_hg38.bed`): Super-enhancer region information
2. **Sample information file** (`human_sample_information_sedb2.txt`): Metadata for each sample
3. **Target regions file** (TSV or BED format): List of genomic regions to analyze

### Command-line Syntax

```bash
sedb-analyzer <command> [options]
```

## Main Commands

### 1. Pre-calculate Statistics (`prepare-stats`)

Pre-calculate and cache database-wide statistics before performing large-scale analyses.

```bash
sedb-analyzer prepare-stats \
  -b data/SE_package_hg38.bed \
  -s data/human_sample_information_sedb2.txt \
  -o cache/sedb_stats.pkl \
  --verbose
```

**Key Options:**
- `-b, --bed-file`: Path to SEdb BED file
- `-s, --sample-info`: Path to sample information file
- `-o, --output`: Output cache file (default: `sedb_stats_cache.pkl`)
- `--no-filter-chromosomes`: Disable filtering for standard human chromosomes
- `-v, --verbose`: Enable verbose logging

### 2. Batch Analysis (`batch`)

Analyze multiple genomic regions in batch mode.

```bash
sedb-analyzer batch \
  -b data/SE_package_hg38.bed \
  -s data/human_sample_information_sedb2.txt \
  -r regions.tsv \
  -o results/ \
  --use-cached-stats cache/sedb_stats.pkl \
  --max-rows 100 \
  --figure-formats png svg \
  --verbose
```

**Key Options:**
- `-r, --regions-file`: Target regions file (TSV or BED format)
- `-o, --output-dir`: Output directory (default: `sedb_results`)
- `--use-cached-stats`: Path to cached statistics file
- `--save-stats-cache`: Save calculated statistics to cache file
- `--max-rows`: Maximum number of regions to process
- `--figure-formats`: Figure output formats (png, svg, pdf, eps)
- `--peak-threshold`: Peak detection threshold (default: 0.2)
- `--pvalue-threshold`: P-value threshold (default: 0.05)
- `--fdr-threshold`: FDR threshold (default: 0.05)
- `--fc-threshold`: Fold change threshold (default: 1.5)
- `--no-tables`: Disable CSV table saving
- `--no-peak-analysis`: Disable peak analysis (for faster processing)

### 3. Single Region Analysis (`single`)

Analyze a specific genomic region in detail.

```bash
sedb-analyzer single \
  -b data/SE_package_hg38.bed \
  -s data/human_sample_information_sedb2.txt \
  --chr chr7 \
  --start 5000000 \
  --end 6000000 \
  --region-name "MYC_enhancer" \
  -o results/single_region \
  --use-cached-stats cache/sedb_stats.pkl \
  --verbose
```

**Key Options:**
- `--chr`: Chromosome name (e.g., chr7)
- `--start`: Start position
- `--end`: End position
- `--region-name`: Region name (optional)
- `-o, --output-dir`: Output directory (default: `sedb_results`)
- `--use-cached-stats`: Path to cached statistics file
- `--figure-formats`: Figure output formats
- `--peak-threshold`: Peak detection threshold
- `--pvalue-threshold`: P-value threshold

### 4. Database Report Generation (`report`)

Generate a comprehensive metadata report for the entire SEdb database.

```bash
sedb-analyzer report \
  -b data/SE_package_hg38.bed \
  -s data/human_sample_information_sedb2.txt \
  -o report/ \
  --top-n 20 \
  --dpi 300 \
  --verbose
```

**Key Options:**
- `-o, --output-dir`: Output directory (default: `sedb_report`)
- `--top-n`: Number of top categories to display (default: 20)
- `--dpi`: DPI for raster images (default: 300)

## Processing Flow

1. **Load and filter SEdb database**
2. **Calculate database-wide statistics** (tissue distribution, etc.)
3. **Use these statistics as background for region-specific analysis**

## Input File Formats

### Regions File (TSV format)

```tsv
chromosome	start	end	name
chr1	1000000	2000000	region1
chr2	5000000	6000000	region2
chr7	27000000	28000000	MYC_region
```

### Regions File (BED format)

```bed
chr1	1000000	2000000	region1
chr2	5000000	6000000	region2
chr7	27000000	28000000	MYC_region
```

## Output Files

### Batch Analysis Output Structure

```
results/
├── data/
│   ├── batch_analysis_results_20250117_143022.tsv
│   └── individual_region_reports/
├── reports/
│   ├── batch_analysis_summary_20250117_143022.html
│   └── figures/
└── cache/
    └── sedb_stats.pkl
```

### Single Region Analysis Output Structure

```
results/
├── chr7_5000000_6000000_analysis_report.html
├── chr7_5000000_6000000_tissue_enrichment.png
├── chr7_5000000_6000000_peak_analysis.png
└── chr7_5000000_6000000_data_tables.csv
```

## Practical Usage Examples

### Efficient Batch Analysis Workflow

```bash
# 1. Pre-calculate and cache statistics
sedb-analyzer prepare-stats \
  -b data/SE_package_hg38.bed \
  -s data/human_sample_information_sedb2.txt \
  -o cache/sedb_stats.pkl \
  --verbose

# 2. Batch analysis using cache
sedb-analyzer batch \
  -b data/SE_package_hg38.bed \
  -s data/human_sample_information_sedb2.txt \
  -r regions.tsv \
  -o results/batch_analysis \
  --use-cached-stats cache/sedb_stats.pkl \
  --max-rows 100 \
  --figure-formats png svg \
  --verbose

# 3. Detailed analysis of specific region
sedb-analyzer single \
  -b data/SE_package_hg38.bed \
  -s data/human_sample_information_sedb2.txt \
  --chr chr8 \
  --start 128000000 \
  --end 129000000 \
  --region-name "MYC_locus" \
  -o results/myc_analysis \
  --use-cached-stats cache/sedb_stats.pkl \
  --verbose
```

### Optimization for Fast Processing

```bash
# Disable peak analysis for faster processing
sedb-analyzer batch \
  -b data/SE_package_hg38.bed \
  -s data/human_sample_information_sedb2.txt \
  -r large_regions.tsv \
  -o results/fast_analysis \
  --use-cached-stats cache/sedb_stats.pkl \
  --no-peak-analysis \
  --no-tables \
  --figure-formats png \
  --max-rows 1000 \
  --verbose
```

## Troubleshooting

### Common Issues and Solutions

1. **Out of Memory Error**
   - Use `--max-rows` option to limit the number of regions processed
   - Use `--use-cached-stats` option to utilize pre-calculated statistics

2. **Long Processing Time**
   - Use `--no-peak-analysis` option to disable peak analysis
   - Use `--use-cached-stats` option to utilize cached statistics

3. **File Format Error**
   - Check file encoding (UTF-8 recommended)
   - Verify TSV/BED format column count and order

4. **No Output Files Generated**
   - Check write permissions for output directory
   - Use `--verbose` option for detailed logging

## Performance Optimization

### Recommended Workflow

1. **Small-scale Testing**: Test with a few regions first
2. **Statistics Caching**: Cache statistics before large-scale analysis
3. **Staged Processing**: Process large datasets in batches
4. **Appropriate Options**: Adjust analysis features as needed

### System Requirements

- **Memory**: 8GB+ (16GB+ recommended for large-scale analysis)
- **Storage**: 10x the size of analysis data in free space
- **CPU**: Multi-core recommended (supports parallel processing)

## Related Documentation

- [Installation Guide](segene_analyzer_installation.md)
- [Usage Examples](../examples/README_en.md)
- [Project Overview](../README.md)

