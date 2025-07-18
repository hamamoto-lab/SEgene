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
- bedtools (required for pybedtools)

### Hardware Requirements
- Memory: 8GB+ (16GB+ recommended for large-scale analysis)
- Storage: 10GB+ free space

## Installation Steps

### 1. Install bedtools

bedtools is required for pybedtools functionality.

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

### 2. Install SEgene_analyzer

#### Method 1: Using Conda/Miniforge (Recommended)

For bioinformatics workflows, conda/miniforge is strongly recommended:

```bash
# Install miniforge (if not already installed)
# Visit: https://github.com/conda-forge/miniforge

# Create conda environment with bedtools
conda create -n segene-analyzer python=3.11
conda activate segene-analyzer

# Install bedtools via conda (recommended)
conda install -c bioconda bedtools pybedtools

# Clone the repository
git clone https://github.com/hamamoto-lab/SEgene.git
cd SEgene/SEgene_analyzer

# Install Python dependencies
pip install pandas matplotlib seaborn scipy statsmodels networkx japanize-matplotlib

# Install in development mode
pip install -e .
```

#### Method 2: Using pip (Alternative)

```bash
# Clone the repository
git clone https://github.com/hamamoto-lab/SEgene.git
cd SEgene/SEgene_analyzer

# Install dependencies
pip install -r requirements.txt

# Install in development mode
pip install -e .
```

### 3. Download SEdb Data

SEgene_analyzer requires SEdb 2.0 database files:

1. **Visit**: [SEdb 2.0 Download Page](http://www.licpathway.net/sedb/download.php)
2. **Download Required Files**:
   - `human_sample_information_sedb2.txt` (Sample information)
   - `SE_package_hg38.bed` (Super-enhancer package for human hg38)

3. **Create Data Directory**:
   ```bash
   mkdir -p data/SEdb
   mv human_sample_information_sedb2.txt data/SEdb/
   mv SE_package_hg38.bed data/SEdb/
   ```

### 4. Verify Installation

```bash
# Check if sedb-analyzer is available
sedb-analyzer --help

# Test with version info
sedb-analyzer --version

# Verify dependencies
python -c "import pandas, numpy, scipy, pybedtools; print('All dependencies installed successfully')"
```

## Installation Verification

### Basic Functionality Test

Create a simple test to verify the installation:

```bash
# Create test directories
mkdir -p test_output cache

# Test statistics preparation (requires SEdb data)
sedb-analyzer prepare-stats \
  -b data/SEdb/SE_package_hg38.bed \
  -s data/SEdb/human_sample_information_sedb2.txt \
  -o cache/test_stats.pkl

# If successful, you should see output without errors
```

### Test with Sample Data

```bash
# Create a sample regions file
cat > test_regions.tsv << EOF
test_region1	chr1	1000000	2000000
test_region2	chr2	3000000	4000000
EOF

# Run a small batch analysis
sedb-analyzer batch \
  -b data/SEdb/SE_package_hg38.bed \
  -s data/SEdb/human_sample_information_sedb2.txt \
  -r test_regions.tsv \
  -o test_output/ \
  --max-rows 2
```

## Troubleshooting

### Common Issues

#### 1. bedtools Not Found
**Error**: `bedtools not found in PATH`

**Solution**:
```bash
# Verify bedtools installation
which bedtools

# If not found, install using appropriate method above
# For conda users:
conda install -c bioconda bedtools
```

#### 2. Python Version Issues
**Error**: `Python 3.11+ required`

**Solution**:
```bash
# Check Python version
python --version

# If using older version, install Python 3.11+
# Using conda:
conda install python=3.11
```

#### 3. Permission Errors
**Error**: `Permission denied when installing`

**Solution**:
```bash
# Use conda environment
conda create -n segene-analyzer python=3.11
conda activate segene-analyzer
pip install -r requirements.txt
pip install -e .
```

#### 4. Memory Issues
**Error**: `Out of memory during analysis`

**Solution**:
```bash
# Use smaller batch sizes
sedb-analyzer batch --max-rows 50 ...

# Pre-calculate statistics
sedb-analyzer prepare-stats ...
sedb-analyzer batch --use-cached-stats ...
```

#### 5. SEdb Data Download Issues
**Error**: `Cannot download SEdb data`

**Solutions**:
- Check internet connection
- Try downloading at different times (server may be busy)
- Use alternative download methods if available

### Performance Optimization

#### 1. Memory Usage
```bash
# Monitor memory usage
htop
# or
top

# Reduce memory usage
sedb-analyzer batch --max-rows 100 --no-peak-analysis ...
```

#### 2. Processing Speed
```bash
# Use statistics caching
sedb-analyzer prepare-stats ...
sedb-analyzer batch --use-cached-stats ...

# Minimize output
sedb-analyzer batch --figure-formats png --no-tables ...
```

## Advanced Installation

### Development Installation

For developers who want to modify the code:

```bash
# Clone with development dependencies
git clone https://github.com/hamamoto-lab/SEgene.git
cd SEgene/SEgene_analyzer

# Install in development mode
pip install -e .

# Install testing dependencies
pip install pytest pytest-cov

# Run tests
python tests/run_tests.py
```


## Updating SEgene_analyzer

### Update from Git

```bash
cd SEgene/SEgene_analyzer
git pull origin main
pip install -r requirements.txt
pip install -e .
```

### Update Dependencies

```bash
# Update all dependencies
pip install --upgrade -r requirements.txt

# Update specific packages
pip install --upgrade pandas numpy scipy
```

## Uninstallation

### Remove Package

```bash
# Uninstall SEgene_analyzer
pip uninstall segene-analyzer

# Remove virtual environment (if used)
rm -rf segene_env

# Remove data files (optional)
rm -rf data/SEdb
```

### Clean Installation

```bash
# Remove cache and temporary files
rm -rf cache/ results/ output/

# Remove __pycache__ directories
find . -name "__pycache__" -type d -exec rm -rf {} +
```


## Related Documentation

- [Usage Guide](usage.md) - Complete usage instructions
- [Usage Guide](usage_ja.md) (Japanese) - Comprehensive command reference
- [Examples](../examples/README.md) - Real-world workflows
- [Main Documentation](../README.md) - Project overview