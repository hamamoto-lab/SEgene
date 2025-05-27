# cpm_peakprep - edgeR Normalized CPM Guide

This document provides an overview and basic usage instructions for the **edgeR normalized CPM** (calcnorm CPM) feature of cpm_peakprep. For general information and the CPM method overview, please refer to the [CPM Method Detailed Guide](./cpm_README.md).

## Table of Contents

1. [Overview](#overview)
2. [Requirements and Dependencies](#requirements-and-dependencies)
3. [Processing Flow](#processing-flow)
4. [Command Line Arguments](#command-line-arguments)
5. [Output File Formats](#output-file-formats)
6. [Usage Examples](#usage-examples)
7. [Related Documentation](#related-documentation)

## Overview

The edgeR‑normalized CPM (calcnorm CPM) feature uses the **edgeR** package in R to normalize ChIP‑seq peak count data. You can enable this functionality by adding the `--calcnorm` option.  
> **Prerequisites:** R ≥ 4.2, edgeR ≥ 3.40, and rpy2 ≥ 3.5 must be installed.

### Differences Between Standard log2-CPM and calcnorm CPM

1. **Standard log2-CPM (traditional method)**:
   - Calculated based on total mapped reads obtained from samtools flagstat
   - Formula: CPM = (count × 10⁶) / total_mapped_reads
   - Then log₂(CPM + pseudocount) transformation is applied

2. **edgeR normalized CPM (calcnorm CPM)**:
   - Calculated based on reads counted in the relevant regions by featureCounts
   - Uses normalization factors calculated by edgeR's `calcNormFactors()` function
   - Applies scaling factors based on count distribution characteristics of each sample

## Requirements and Dependencies

To use this feature, the following software is required:

1. **R language**: Version 4.2.0 or higher
2. **rpy2**: Version 3.5.0 or higher (Python-R interface library)
3. **edgeR**: Version 3.40.0 or higher (Bioconductor package)

There are multiple ways to install R/edgeR. For example, using conda-forge:

```bash
# Installation example with conda-forge
conda install -c conda-forge r-base=4.2
conda install -c conda-forge r-essentials
conda install -c bioconda bioconductor-edger=3.40
conda install -c conda-forge rpy2

# Alternative method: installation with pip
# pip install rpy2>=3.5.0
```

For installation methods specific to your environment, please refer to the official R/Bioconductor documentation.

## Processing Flow

The edgeR normalized CPM processing follows these steps:

1. **Data Preprocessing**:
   - Loading count data from featureCounts output file
   - Cleaning sample names (extracting meaningful sample names from BAM filenames)

2. **DGEList Creation**:
   - Converting count data to edgeR's `DGEList` object

3. **Filtering** (optional):
   - Selecting peaks that meet specific criteria based on CPM values
   - Parameters: `min_cpm` (minimum CPM value) and `min_samples` (minimum number of samples)

4. **Calculating Normalization Factors**:
   - Computing normalization factors using the selected method
   - User-selectable methods: upperquartile (default), TMM, RLE, none

5. **CPM Calculation** (without log transformation):
   - Computing CPM values based on normalized library sizes

6. **Result Output**:
   - Saving in standard format (PeakID, Chr, Start, End, values for each sample)
   - Optionally saving a complete version with all metadata columns

### About Normalization Methods

The edgeR package provides the following main normalization methods:

- **Upper Quartile** (default): Normalization based on the 75th percentile of the count distribution in each sample
- **TMM**: Trimmed Mean of M-values. Normalization based on the trimmed mean of log ratios between samples
- **RLE**: Relative Log Expression. Normalization using the geometric mean of peaks as a reference
- **none**: No normalization factors applied

For detailed theoretical background on each normalization method, please refer to the edgeR documentation.

### About Filtering

The `min_cpm` and `min_samples` parameters allow you to select only peaks that meet specific criteria:

- A peak is retained if the number of samples with "CPM > min_cpm" is at least "min_samples"
- The default settings are "min_cpm=1.0, min_samples=0", which effectively performs no filtering

## Command Line Arguments

The main command line arguments available for the edgeR normalized CPM feature are as follows:

### Basic Arguments

| Argument | Description | Default Value |
|--------|------|-------------|
| `--calcnorm` | Flag to enable edgeR normalized CPM | `False` |
| `--calcnorm_method` | Normalization method: upperquartile, TMM, RLE, none | `upperquartile` |
| `--calcnorm_output_name` | calcnorm CPM output filename. Empty string to disable saving | `calcnorm.tsv` |

### Filtering-related Arguments

| Argument | Description | Default Value |
|--------|------|-------------|
| `--min_cpm` | Minimum CPM threshold for filtering | `1.0` |
| `--min_samples` | Minimum number of samples that should exceed CPM threshold. 0=no filtering | `0` |

### Sample Name Cleaning Arguments

| Argument | Description | Default Value |
|--------|------|-------------|
| `--remove_extensions` | Automatically remove common BAM file extensions (`.bam`, `.sorted.bam`, etc.) | `False` |
| `--pattern_rules` | JSON file containing additional pattern rules for sample name cleaning | `None` |

### Output Format Arguments

| Argument | Description | Default Value |
|--------|------|-------------|
| `--full_metadata_output` | Save output with all featureCounts metadata columns in addition to standard output | `None` |

## Output File Formats

### Standard Output File (calcnorm.tsv)

The standard output format has the following structure:

```
PeakID  Chr     Start   End     Sample1  Sample2  Sample3
peak1   chr1    999     2000    2.45     3.12     1.87
peak2   chr1    2999    4000    4.21     3.98     4.56
...
```

- **Note**: The `Start` column is output in 0-based coordinates (BED format standard)

### Full Metadata Output File

If `--full_metadata_output` is specified, a complete output with all featureCounts metadata columns is generated:

```
Geneid  Chr     Start   End     Strand  Length  Sample1  Sample2  Sample3
peak1   chr1    999     2000    +       1001    2.45     3.12     1.87
peak2   chr1    2999    4000    +       1001    4.21     3.98     4.56
...
```

This output format is useful when additional metadata like Strand information or peak length is needed.

### JSON File for Sample Name Cleaning

For complex sample name patterns, you can specify a JSON file with the `--pattern_rules` option:

```json
[
  {"pattern": "\\.mLb\\.clN\\.sorted\\.bam$", "replace": ""},
  {"pattern": "^Sample_", "replace": ""},
  {"pattern": "-rep([0-9])", "replace": "_R$1"}
]
```

Each rule consists of a pair of `pattern` (regular expression pattern) and `replace` (replacement string).

## Usage Examples

### Basic Example Saving Both CPM Values

```bash
python cpm_peakprep.py \
    --bam_dir /path/to/bam_files \
    --annotation_file peaks.bed \
    --output_dir results \
    --calcnorm \
    --threads 10
```

This example saves both Standard log2-CPM (`logCPM.tsv`) and calcnorm CPM (`calcnorm.tsv`).

### Example Using TMM Normalization

```bash
python cpm_peakprep.py \
    --bam_dir /path/to/bam_files \
    --annotation_file peaks.bed \
    --output_dir results \
    --calcnorm \
    --calcnorm_method TMM \
    --threads 10
```

### Example with Filtering Applied

```bash
python cpm_peakprep.py \
    --bam_dir /path/to/bam_files \
    --annotation_file peaks.bed \
    --output_dir results \
    --calcnorm \
    --min_cpm 1.0 \
    --min_samples 3 \
    --threads 10
```

In this example, only peaks with CPM values above 1.0 in at least 3 samples will be retained.

### Example Saving Full Metadata Output

```bash
python cpm_peakprep.py \
    --bam_dir /path/to/bam_files \
    --annotation_file peaks.bed \
    --output_dir results \
    --calcnorm \
    --full_metadata_output \
    --threads 10
```

## Related Documentation

- [Main README](./README.md): Overview of the entire tool
- [CPM Method Detailed Guide](./cpm_README.md): Basic usage and detailed information for the CPM method
- [BigWig Method Detailed Documentation](./bigwig_README.md): Detailed information on the BigWig method

## Citation

If you use this tool in your research, please cite:

Shinkai, N., Asada, K., Machino, H., Takasawa, K., Takahashi, S., Kouno, N., Komatsu, M., Hamamoto, R., & Kaneko, S. (2025). SEgene identifies links between super enhancers and gene expression across cell types. *npj Systems Biology and Applications*, 11(1), 49. https://doi.org/10.1038/s41540-025-00533-x

## License

This program is released under the MIT license. For details, see [LICENSE](https://github.com/hamamoto-lab/SEgene/blob/main/LICENSE).