# cpm_peakprep - CPM Method Detailed Guide

This document provides detailed usage instructions and examples specific to the **CPM method** of SEgene_peakprep. For basic overview and common requirements, please refer to the [main README](./README.md).

## Complete Argument List

Here is a comprehensive list of all arguments available for command-line execution:

| Argument | Description | Default value |
|--------|------|-------------|
| **Required Arguments** |||
| `--bam_dir` / `-b` | Directory containing BAM files | - |
| `--annotation_file` / `-a` | Path to the annotation file (BED, SAF, or merge_SE.tsv format) | - |
| `--output_dir` / `-o` | Directory to save results | - |
| **Annotation-related Options** |||
| `--is_mergese_format` | Flag to specify that the annotation file is in merge_SE.tsv format | `False` |
| **File/Sample Selection Options** |||
| `--filename_pattern` | Wildcard pattern for filtering BAM filenames | `*.bam` |
| `--sample_delimiter` | Delimiter string to extract sample names from BAM filenames | `None` |
| **Tool Paths and Thread Count** |||
| `--samtools_path` | Path to samtools executable | `samtools` |
| `--featurecounts_path` | Path to featureCounts executable | `featureCounts` |
| `--threads` / `-T` | Number of threads to use for samtools and featureCounts | `4` |
| `--single_end` | Flag to specify single-end reads (if not specified, paired-end is assumed) | `False` |
| **featureCounts-related Options** |||
| `--fc_basename` | Base name for featureCounts output/log files | `None` |
| `--fc_options` | Additional options to pass to featureCounts (e.g., `--fc_options --minOverlap 10`) | `None` |
| **Standard log2-CPM Calculation Options** |||
| `--add_peakid` | Replace original GeneIDs with sequential PeakIDs | `False` |
| `--id_column` | ID column name in SAF/featureCounts output | `Geneid` |
| `--output_id_name` | ID column name in the final logCPM table | `PeakID` |
| `--log_base` | Base for logarithmic transformation (2 or 10). Values ≤0 disable transformation | `2` |
| `--pseudocount` | Pseudocount to add before logarithmic transformation | `1.0` |
| **Output Filename Options** |||
| `--logcpm_output_name` | Filename for the log2-CPM table | `logCPM.tsv` |
| **Logging Options** |||
| `--log_level` | Logging level (DEBUG/INFO/WARNING/ERROR/CRITICAL) | `INFO` |
| `--script_log_file` | Filename for script execution log (empty string to disable) | `pipeline_run.log` |

## Dependencies and Version Requirements

### Python
- **Python**: ≥3.10
- **pandas**: ≥1.5
- **numpy**: ≥1.23
- **natsort**: ≥8.3

### For Standard CPM Calculation Only
- **samtools**: ≥1.13
- **featureCounts** (Subread): ≥2.0.0

### Additional for edgeR normalized CPM (calcnorm)
- **R language**: ≥4.2.0
- **edgeR**: ≥3.40.0
- **rpy2**: ≥3.5.0

## edgeR normalized CPM (calcnorm) Mode

Adding the `--calcnorm` flag applies edgeR's `calcNormFactors()` function to the raw count data from featureCounts, calculating CPM values based on normalized library sizes. This feature allows for normalization that corrects for compositional biases between samples.

By default, the script outputs both **Standard log2-CPM** (`logCPM.tsv`) and **calcnorm CPM** (`calcnorm.tsv`). The output of each can be controlled with `--logcpm_output_name` and `--calcnorm_output_name`.

| Argument | Description | Default value |
|--------|------|-------------|
| `--calcnorm` | Enable edgeR normalized CPM | `False` |
| `--calcnorm_method` | Normalization method (upperquartile, TMM, RLE, none) | `upperquartile` |
| `--min_cpm` | Minimum CPM threshold for filtering | `1.0` |
| `--min_samples` | Minimum number of samples that should exceed CPM threshold (0=no filtering) | `0` |
| `--calcnorm_output_name` | calcnorm CPM output filename (empty string to disable saving) | `calcnorm.tsv` |
| `--remove_extensions` | Remove common BAM file extensions (.bam, .sorted.bam, etc.) | `False` |
| `--pattern_rules` | JSON file containing additional pattern rules for sample name cleaning | `None` |
| `--full_metadata_output` | Save output with all featureCounts metadata columns in addition to standard output | `None` |

### Choosing a Normalization Method

The `--calcnorm_method` option allows you to select one of the normalization methods implemented in the edgeR package:

- **upperquartile** (default): Upper quartile normalization. Based on the 75th percentile of non-zero count values, suitable for sparse datasets like ChIP-seq data.
- **TMM**: Trimmed Mean of M-values. Based on trimmed mean of log ratios between samples, commonly used in RNA-seq data analysis.
- **RLE**: Relative Log Expression. Based on relative expression levels using the geometric mean of each gene as a reference.
- **none**: No normalization applied, only computes CPM considering library size differences.

For detailed theoretical background on each normalization method, please refer to the edgeR documentation. For more details on calcnorm CPM usage and background, see the [edgeR Normalized CPM Guide](./cpm_calcnorm_README.md).

## Differences in CPM Value Calculation

### Standard log2-CPM vs. edgeR normalized CPM

1. **Standard log2-CPM**:
   - Calculated based on total mapped reads obtained from samtools flagstat
   - Formula: CPM = (count × 10⁶) / total_mapped_reads
   - Then log₂(CPM + pseudocount) transformation is applied

2. **edgeR normalized CPM (calcnorm CPM)**:
   - Calculated based on reads counted in the regions by featureCounts
   - Uses edgeR normalization factors to adjust library sizes between samples
   - Allows for normalization that accounts for differences in sample composition

These differences can result in different values, especially when comparing samples with varying compositions.

## Basic Usage Examples

### Simple Execution Example

```bash
python cpm_peakprep.py \
    --bam_dir /path/to/bam_files \
    --annotation_file peaks.bed \
    --output_dir results \
    --filename_pattern "Sample*.sorted.bam" \
    --sample_delimiter=".sorted.bam" \
    --threads 10
```

### Specifying Additional Options for featureCounts

```bash
python cpm_peakprep.py \
    --bam_dir /path/to/bam_files \
    --annotation_file peaks.saf \
    --output_dir results \
    --fc_options --minOverlap 10 --fracOverlap 0.2 --ignoreDup \
    --threads 8
```

### Custom ID Generation and Log Transformation Settings

```bash
python cpm_peakprep.py \
    --bam_dir /path/to/bam_files \
    --annotation_file peaks.bed \
    --output_dir results \
    --add_peakid \
    --output_id_name "CustomPeakID" \
    --log_base 10 \
    --pseudocount 0.1
```

### Using edgeR normalized CPM

```bash
python cpm_peakprep.py \
    --bam_dir /path/to/bam_files \
    --annotation_file peaks.bed \
    --output_dir results \
    --calcnorm \
    --calcnorm_method TMM \
    --threads 10
```

### Saving calcnorm CPM Only

```bash
python cpm_peakprep.py \
    --bam_dir /path/to/bam_files \
    --annotation_file peaks.bed \
    --output_dir results \
    --calcnorm \
    --logcpm_output_name "" \
    --threads 10
```

### Applying Expression Filtering

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

## Jupyter Notebook Execution Examples

### Simple Usage

In its simplest form, you can import and use the main script:

```python
import sys
sys.path.append('/path/to/SEgene/SEgene_peakprep')

# Run as a module
import cpm_peakprep

# Set parameters
args = [
    "--bam_dir", "/path/to/bam_files",
    "--annotation_file", "peaks.bed",
    "--output_dir", "results",
    "--sample_delimiter", ".sorted.bam",
    "--threads", "8"
]

# Call main function directly
sys.argv = ["cpm_peakprep.py"] + args
cpm_peakprep.main()
```

### Detailed Step-by-Step Execution Example

Here's a detailed example that controls each step individually:

```python
import os
import sys
import logging
import pandas as pd

# Set up path for module import
sys.path.append('/path/to/SEgene/SEgene_peakprep')

# Import specific functions
from cpm_peakprep_utils import (
    get_sample_data,
    run_samtools_flagstat,
    run_featurecounts,
    calculate_logcpm,
    convert_bed_to_saf
)

# Set parameters
bam_dir = "/path/to/bam_files"
annotation_file = "peaks.bed"
output_dir = "results"
filename_pattern = "Sample*.sorted.bam"
sample_delimiter = ".sorted.bam"
threads = 10

# Create output directory
os.makedirs(output_dir, exist_ok=True)

# Configure logger
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger("ChIPseqPipeline")

# Convert BED file to SAF (if needed)
saf_file = annotation_file
if annotation_file.lower().endswith('.bed'):
    saf_file = convert_bed_to_saf(annotation_file, logger)
    if saf_file is None:
        print("Failed to convert BED file.")
        sys.exit(1)

# Step 1: Get BAM file information
sample_dict, bam_list = get_sample_data(
    bam_folder=bam_dir,
    logger=logger,
    sample_name_delimiter=sample_delimiter,
    filename_pattern=filename_pattern
)

# Step 2: Get total read counts
flagstat_dir = os.path.join(output_dir, "flagstat")
total_reads_map = run_samtools_flagstat(
    bam_files=bam_list,
    output_dir=flagstat_dir,
    logger=logger,
    threads=threads
)

# Step 3: Count reads (featureCounts)
fc_output = os.path.join(output_dir, f"{os.path.splitext(os.path.basename(saf_file))[0]}_featureCounts.txt")
fc_log = os.path.join(output_dir, f"{os.path.splitext(os.path.basename(saf_file))[0]}_featureCounts.log")
counts_file_path = run_featurecounts(
    bam_files=bam_list,
    saf_file=saf_file,
    logger=logger,
    output_file=fc_output,
    log_file=fc_log,
    threads=threads,
    is_paired_end=True
)

# Load counts file and Step 4: Calculate LogCPM
counts_df = pd.read_csv(counts_file_path, sep='\t', comment='#', header=0)
if counts_df.columns[0].startswith('#'):
    counts_df = pd.read_csv(counts_file_path, sep='\t', comment='#', header=1)

bam_to_sample_map = {v: k for k, v in sample_dict.items()}
logcpm_df = calculate_logcpm(
    counts_df=counts_df,
    total_reads_dict=total_reads_map,
    logger=logger,
    bam_path_to_sample_name=bam_to_sample_map,
    log_transform_base=2,
    pseudocount=1.0
)

# Save logCPM results
logcpm_output_file = os.path.join(output_dir, "logCPM.tsv")
logcpm_df.to_csv(logcpm_output_file, sep='\t', index=False, float_format='%.4f')
print(f"logCPM table saved: {logcpm_output_file}")
```

## Intermediate Files and Output Formats

### flagstat Output

From the samtools flagstat results saved in the `flagstat/` directory, the number of mapped reads is extracted and used for CPM calculation.

### featureCounts Output

Example featureCounts output:
```
# Program:featureCounts v2.0.6; Command:...
Geneid  Chr     Start   End     Strand  Length  bam_files/sample1.bam  bam_files/sample2.bam
peak1   chr1    1000    2000    .       1001    150     220
peak2   chr1    3000    4000    .       1001    342     289
...
```

### log2-CPM Table

Example final log2-CPM table:
```
PeakID  Chr     Start   End     Sample1  Sample2  Sample3
peak1   chr1    999     2000    2.45     3.12     1.87
peak2   chr1    2999    4000    4.21     3.98     4.56
...
```

> **Note**: The `Start` column in all output files is 0-based (BED format standard). Even if the input is in SAF format (1-based), the output is converted to 0-based.

## Related Documentation

- [Main README](./README.md): Overview of the entire tool
- [edgeR Normalized CPM Guide](./cpm_calcnorm_README.md): Detailed information on edgeR normalized CPM
- [BigWig Method Detailed Documentation](./bigwig_README.md): Detailed information on the BigWig method

## Citation

If you use this tool in your research, please cite:
**(Paper in preparation)**

## License

This program is released under the MIT license. For details, see [LICENSE](https://github.com/hamamoto-lab/SEgene/blob/main/LICENSE)