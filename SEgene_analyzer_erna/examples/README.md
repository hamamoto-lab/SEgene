# Examples

This directory contains example scripts demonstrating how to use SEgene_analyzer_erna.

## Available Examples

### 1. Command-Line Interface (`cli_usage.sh`)

Comprehensive examples of using the `erna-analyzer` CLI tool:
- Single region analysis
- Batch processing with/without cache
- Statistics pre-calculation for performance
- Metadata report generation
- Various output formats and thresholds

```bash
# View examples
cat cli_usage.sh

# Run a specific example
erna-analyzer single -b data/peak -m metadata/eRNAbase_data.parquet --chr chr7 --start 1000000 --end 2000000
```

### 2. Basic Usage (`basic_usage.py`)

Demonstrates the fundamental workflow for analyzing a single genomic region:
- Initialize the ERNARegionAnalyzer
- Load eRNAbase metadata
- Load BED files
- Analyze a specific region
- Generate reports

```bash
python basic_usage.py
```

### 3. Batch Analysis (`batch_analysis.py`)

Shows how to analyze multiple genomic regions using a TSV input file:
- Create or use an existing regions TSV file
- Perform batch analysis on multiple regions
- Generate comprehensive reports for all regions

```bash
python batch_analysis.py
```

### 4. Production Workflow (`production_workflow.py`)

Production-ready script with error handling and logging:
- Command-line argument parsing
- Comprehensive error handling
- Progress tracking
- Detailed logging

```bash
python production_workflow.py --bed-dir data/peak --metadata metadata/eRNAbase_data.parquet --regions regions.tsv
```

### 5. Sample Data Generation (`create_sample_data.py`)

Generate test data for development and testing:
- Create synthetic BED files
- Generate metadata in both CSV and Parquet formats
- Create test region lists

```bash
python create_sample_data.py --num-samples 10 --num-regions 5
```

## Input Data Requirements

### eRNAbase Metadata
- Format: CSV or Parquet
- Required columns: sample information, tissue types, cell types, etc.

### BED Files
- Standard BED format with genomic coordinates
- Directory containing multiple BED files for different samples

### Regions TSV (for batch analysis)
- Tab-separated file with columns: `region_name`, `chrom`, `start`, `end`
- Example:
  ```
  region_name	chrom	start	end
  region_1	chr1	1000000	2000000
  region_2	chr2	5000000	6000000
  ```

## Running the Examples

1. **Prepare your data**: Ensure you have eRNAbase metadata and BED files
2. **Update paths**: Modify the file paths in the examples to point to your actual data
3. **Install dependencies**: Make sure all required packages are installed
4. **Run the scripts**: Execute the examples from this directory

## Output

The examples will create output directories containing:
- Overlapping eRNA data (CSV format)
- Statistical analysis results
- Visualization plots (PNG format)
- HTML reports with comprehensive analysis