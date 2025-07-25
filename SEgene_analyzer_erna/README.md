# SEgene_analyzer_erna (Development Version)

*(For the Japanese version of this README, please see [README_ja.md](README_ja.md).)*

**SEgene_analyzer_erna** is an **eRNAbase-specific** genomic region analysis tool designed for analyzing enhancer RNA regions using eRNAbase data with enhanced functionality and modern Python packaging.

> **⚠️ Development Status**: This is a **specialized development version** for eRNAbase analysis. For SEdb-based analysis, please consider using [SEgene_analyzer](https://github.com/hamamoto-lab/SEgene/tree/main/SEgene_analyzer).

[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Project Context

SEgene_analyzer_erna is part of the **SEgene project** ecosystem, specifically designed for eRNAbase data analysis.

### Related Components

- **Parent Project**: [SEgene](https://github.com/hamamoto-lab/SEgene) - Complete super-enhancer analysis platform
- **SEdb Version**: [SEgene_analyzer](https://github.com/hamamoto-lab/SEgene/tree/main/SEgene_analyzer) - SEdb-based region analyzer
- **Data Source**: [eRNAbase](https://bio.liclab.net/eRNAbase/index.php) - Enhancer RNA database

## About eRNAbase

This analysis tool uses data obtained from **eRNAbase** (https://bio.liclab.net/eRNAbase/index.php), a comprehensive database of enhancer RNAs (eRNAs). eRNAbase provides experimentally validated eRNA datasets across multiple species and tissue types.

### Data Requirements

**Required eRNA BED files:**
- Download `download_peak.zip` from the [eRNAbase download page](https://bio.liclab.net/eRNAbase/download.php)
- Extract BED files to a directory (e.g., `data/peak/`)

**Required metadata files:**
- Use the metadata files in `metadata/`:
  - `eRNAbase_data.csv` - CSV format
  - `eRNAbase_data.parquet` - Parquet format

**Data Version:** This tool has been developed and tested against eRNAbase data as of **July 22, 2025**.

## Key Features

### Enhanced eRNAbase Analysis
- **Multiple BED files support**: Handle directory containing multiple eRNAbase BED files
- **Parquet metadata integration**: Efficient metadata processing using Parquet format
- **Species-specific analysis**: Human and other species support with proper filtering

### Core Analysis Capabilities
- **Batch Processing**: Analyze multiple genomic regions from TSV files
- **Single Region Analysis**: Detailed analysis of specific genomic regions
- **Tissue Enrichment Analysis**: Statistical evaluation of tissue-specific eRNA enrichment
- **Database Integration**: eRNAbase integration with comprehensive metadata support

### Output Options
- **Multiple formats**: PNG, SVG, EPS figure formats
- **Comprehensive reports**: HTML reports with detailed analysis results
- **Data export**: TSV, CSV, and Excel format support
- **Flexible organization**: Hierarchical output directory structure

### Developer-Friendly Features
- **Modern Python Packaging**: Python 3.11+ with enhanced type hints
- **Method chaining**: Fluent API design for streamlined workflows
- **Comprehensive logging**: Detailed logging for debugging and monitoring
- **Test suite**: Comprehensive unit tests for quality assurance

## Features

### eRNAbase-Specific Functionality
- **Multi-sample BED processing**: Handle hundreds of eRNAbase BED files efficiently
- **Metadata correlation**: Correlate genomic overlaps with tissue/cell type metadata
- **Species filtering**: Built-in human/species-specific filtering capabilities
- **Statistical enrichment**: Fisher's exact test with FDR correction

### Visualization Capabilities
- **Read distribution plots**: Stacked visualization of eRNA signal distribution
- **Enrichment volcano plots**: Statistical significance visualization
- **Tissue/cell type bar charts**: Distribution analysis with Japanese font support
- **Network analysis**: Comprehensive network visualization for complex relationships

## Installation

### Prerequisites

- **Python**: 3.11+
- **System Tools**: bedtools

### Installation Steps

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
git clone https://github.com/hamamoto-lab/SEgene_analyzer_erna.git
cd SEgene_analyzer_erna

# Install Python dependencies
pip install pandas matplotlib seaborn scipy statsmodels networkx japanize-matplotlib natsort Jinja2

# Install in development mode
pip install -e .
```

#### Method 2: Using pip (Alternative)

```bash
# Clone the repository
git clone https://github.com/hamamoto-lab/SEgene_analyzer_erna.git
cd SEgene_analyzer_erna

# Install dependencies
pip install -r requirements.txt

# Install in development mode
pip install -e .
```

### Data Setup

Prepare required eRNAbase data files:

1. **eRNAbase BED files**: Directory containing multiple BED files (one per sample)
2. **Metadata file**: Parquet or CSV file containing sample information
3. **Region list**: TSV file containing genomic regions to analyze

#### Example Data Structure
```
data/
├── peak/                           # BED files directory
│   ├── Sample-01-0001.bed
│   ├── Sample-01-0002.bed
│   └── ...
├── eRNAbase_data.parquet          # Metadata file
└── target_regions.tsv             # Analysis regions
```

## Usage

### Command-Line Interface (CLI)

After installation, the `erna-analyzer` command becomes available:

```bash
# Show help
erna-analyzer --help

# Analyze single region
erna-analyzer single -b data/peak -m metadata/eRNAbase_data.parquet \
  --chr chr7 --start 1000000 --end 2000000 \
  --species human -o results/

# Batch analysis
erna-analyzer batch -b data/peak -m metadata/eRNAbase_data.parquet \
  -r regions.tsv -o results/ --max-rows 10

# Pre-calculate statistics for faster analysis
erna-analyzer prepare-stats -b data/peak -m metadata/eRNAbase_data.parquet \
  -o cache/stats.pkl

# Generate metadata report
erna-analyzer report -m metadata/eRNAbase_data.parquet -o report/
```

### Python API

```python
from src.erna_analyzer import ERNARegionAnalyzer
# Or after installation: from erna_analyzer import ERNARegionAnalyzer

# Initialize analyzer
analyzer = ERNARegionAnalyzer(results_dir='results/erna_analysis')

# Load eRNAbase metadata
analyzer.load_erna_metadata('metadata/eRNAbase_data.parquet', species='human')

# Batch analysis from TSV file
result_file = analyzer.batch_analyze_regions_from_tsv(
    bed_directory='data/peak',                     # Directory containing BED files
    metadata_file='metadata/eRNAbase_data.parquet',  # eRNAbase metadata
    regions_tsv_file='data/target_regions.tsv',   # Analysis regions
    data_output_dir='results/data',
    report_output_dir='results/reports',
    max_rows=10,                                  # Limit for testing
    species='human',                              # Species filtering
    figure_formats=['png', 'svg'],
    save_tables=True
)
```

### Advanced Usage

For detailed examples and advanced usage patterns, see:
- **[Installation Guide](docs/segene_analyzer_erna_installation.md)**: Detailed installation instructions
- **[Usage Examples](examples/README.md)**: Comprehensive usage examples
- **[Test Strategy](TEST_STRATEGY.md)**: Testing and validation approaches

## Data Requirements vs SEgene_analyzer

| Aspect | SEgene_analyzer_erna | SEgene_analyzer |
|--------|----------------------|------------------|
| **Data Source** | eRNAbase | SEdb |
| **Input Format** | Multiple BED files | Single BED file |
| **Metadata** | Parquet/CSV | TSV |
| **Sample Scale** | Hundreds of samples | Thousands of samples |
| **Use Case** | eRNA-specific analysis | General SE analysis |

## Development

### Development Environment Setup

```bash
# Install in development mode
pip install -e .

# Install additional development packages
pip install pytest black flake8
```

### Code Quality

```bash
# Code formatting
black src/

# Linting
flake8 src/
```

### Testing

```bash
# Run basic tests
python tests/test_erna_analyzer.py

# Test with sample data
python examples/create_sample_data.py --num-samples 10
python examples/production_workflow.py --bed-dir sample_data/peak --metadata sample_data/eRNAbase_data.parquet
```

## Related Projects

- [SEgene](https://github.com/hamamoto-lab/SEgene) - Parent project
- [SEgene_analyzer](https://github.com/hamamoto-lab/SEgene/tree/main/SEgene_analyzer) - SEdb analysis tool (foundation of this tool)

## License

This program is released under the MIT License. For more details, please refer to the LICENSE file.

## Citation

If you use this tool in your research, please cite:

Shinkai, N., Asada, K., Machino, H., Takasawa, K., Takahashi, S., Kouno, N., Komatsu, M., Hamamoto, R., & Kaneko, S. (2025). SEgene identifies links between super enhancers and gene expression across cell types. *npj Systems Biology and Applications*, 11(1), 49. https://doi.org/10.1038/s41540-025-00533-x

For detailed citation information and additional references, please refer to the [CITATION](https://github.com/hamamoto-lab/SEgene/blob/main/CITATION) file.

