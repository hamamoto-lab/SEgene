# Installation Guide

## System Requirements

### Operating System
- Linux (Ubuntu 18.04+ recommended)
- macOS (10.14+)
- Windows (WSL2 recommended)

### Software Requirements
- Python 3.11+
- pip 20.0+
- Git
- bedtools

## Installation Steps

### 1. Install bedtools

#### Ubuntu/Debian
```bash
sudo apt-get update
sudo apt-get install bedtools
```

#### macOS (Homebrew)
```bash
brew install bedtools
```

#### CentOS/RHEL
```bash
sudo yum install bedtools
```

#### Conda Environment
```bash
conda install -c bioconda bedtools
```

### 2. Install SEgene_analyzer_erna

#### Method 1: Using Conda/Miniforge (Recommended)

```bash
# Install miniforge (if not already installed)
# Visit: https://github.com/conda-forge/miniforge

# Create conda environment with bedtools
conda create -n segene-analyzer-erna python=3.11
conda activate segene-analyzer-erna

# Install bedtools
conda install -c bioconda bedtools pybedtools

# Clone the repository
git clone https://github.com/your-org/SEgene_analyzer_erna.git
cd SEgene_analyzer_erna

# Install Python dependencies
pip install pandas matplotlib seaborn scipy statsmodels networkx japanize-matplotlib natsort Jinja2

# Install in development mode
pip install -e .
```

#### Method 2: Using pip (Alternative)

```bash
# Clone the repository
git clone https://github.com/your-org/SEgene_analyzer_erna.git
cd SEgene_analyzer_erna

# Install dependencies
pip install -r requirements.txt

# Install in development mode
pip install -e .
```

### 3. Prepare eRNAbase Data

SEgene_analyzer_erna requires eRNAbase database files:

1. **eRNAbase BED files**: Directory containing multiple BED files (one per sample)
2. **Metadata file**: Parquet or CSV file containing sample metadata
3. **Region list file**: TSV file containing genomic regions to analyze

#### Example Data Structure
```
data/
├── peak/                           # BED files directory
│   ├── Sample-01-0001.bed
│   ├── Sample-01-0002.bed
│   └── ...
├── eRNAbase_data.parquet          # Metadata file
└── regions.tsv                    # Analysis target regions
```

### 4. Verify Installation

```bash
# Test import (development mode)
python -c "from erna_analyzer import ERNARegionAnalyzer; print('✅ Installation successful')"

# Or test source import
python -c "from src.erna_analyzer import ERNARegionAnalyzer; print('✅ Installation successful')"
```

## Usage Examples

### Basic Usage
```python
from src.erna_analyzer import ERNARegionAnalyzer
# Or after pip install: from erna_analyzer import ERNARegionAnalyzer

# Initialize analyzer
analyzer = ERNARegionAnalyzer(results_dir='results/erna_analysis')

# Load eRNAbase metadata
analyzer.load_erna_metadata('data/eRNAbase_data.parquet', species='human')

# Analyze regions from TSV file
result_file = analyzer.batch_analyze_regions_from_tsv(
    bed_directory='data/peak',
    metadata_file='data/eRNAbase_data.parquet',
    regions_tsv_file='data/regions.tsv',
    data_output_dir='results/data',
    report_output_dir='results/reports',
    max_rows=10,
    species='human'
)
```

## Troubleshooting

### Common Issues

1. **Import Error: bedtools not found**
   - Solution: Install bedtools using your system package manager or conda

2. **Memory Issues with Large Datasets**
   - Solution: Increase system memory or reduce batch size using `max_rows` parameter

3. **Permission Errors**
   - Solution: Ensure write permissions for output directories

