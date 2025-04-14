# SEgene_peakprep - CPM Method Detailed Guide

This document provides detailed usage instructions and examples specific to the **CPM method** of SEgene_peakprep. For basic overview and common requirements, please refer to the [main README](./README.md).

## Complete Argument List

Here is a comprehensive list of all arguments available for command-line execution:

| Argument | Short form | Default value | Description |
|--------|--------|------------|------|
| **Required Arguments** ||||
| `--bam_dir` | `-b` | - | Directory containing BAM files |
| `--annotation_file` | `-a` | - | Path to the annotation file (BED, SAF, or merge_SE.tsv format) |
| `--output_dir` | `-o` | - | Directory to save results |
| **Annotation-related Options** ||||
| `--is_mergese_format` | - | `False` | Flag to specify that the annotation file is in merge_SE.tsv format |
| **File/Sample Selection Options** ||||
| `--filename_pattern` | - | `*.bam` | Wildcard pattern for filtering BAM filenames |
| `--sample_delimiter` | - | `None` | Delimiter string to extract sample names from BAM filenames |
| **Tool Paths and Thread Count** ||||
| `--samtools_path` | - | `samtools` | Path to samtools executable |
| `--featurecounts_path` | - | `featureCounts` | Path to featureCounts executable |
| `--threads` | `-T` | `4` | Number of threads to use for samtools and featureCounts |
| `--single_end` | - | `False` | Flag to specify single-end reads (if not specified, paired-end is assumed) |
| **featureCounts-related Options** ||||
| `--fc_basename` | - | `None` | Base name for featureCounts output/log files |
| `--fc_options` | - | `None` | Additional options to pass to featureCounts (e.g., `--fc_options --minOverlap 10`) |
| **logCPM Calculation Options** ||||
| `--add_peakid` | - | `False` | Replace original GeneIDs with sequential PeakIDs |
| `--id_column` | - | `Geneid` | ID column name in SAF/featureCounts output |
| `--output_id_name` | - | `PeakID` | ID column name in the final logCPM table |
| `--log_base` | - | `2` | Base for logarithmic transformation (2 or 10). Values ≤0 disable transformation |
| `--pseudocount` | - | `1.0` | Pseudocount to add before logarithmic transformation |
| **Output Filename Options** ||||
| `--logcpm_output_name` | - | `logCPM.tsv` | Filename for the logCPM table |
| **Logging Options** ||||
| `--log_level` | - | `INFO` | Logging level (DEBUG/INFO/WARNING/ERROR/CRITICAL) |
| `--script_log_file` | - | `pipeline_run.log` | Filename for script execution log (empty string to disable) |

## Advanced Use Case Examples

### Specifying Additional Options for featureCounts

Example specifying specific overlap criteria and counting methods:

```bash
python cpm_peakprep.py \
    --bam_dir /path/to/bam_files \
    --annotation_file peaks.saf \
    --output_dir results \
    --fc_options --minOverlap 10 --fracOverlap 0.2 --ignoreDup \
    --threads 8
```

### Custom ID Generation and Log Transformation Settings

Example generating sequential PeakIDs and custom log transformation:

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

# Clean up temporary files
if 'saf_file' in locals() and saf_file != annotation_file and os.path.exists(saf_file):
    try:
        os.remove(saf_file)
        print(f"Temporary SAF file deleted: {saf_file}")
    except Exception as e:
        print(f"Failed to delete temporary file: {e}")
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

### logCPM Table

Example final logCPM table:
```
PeakID  Chr     Start   End     Sample1  Sample2  Sample3
peak1   chr1    999     2000    2.45     3.12     1.87
peak2   chr1    2999    4000    4.21     3.98     4.56
...
```

> **Note**: The `Start` column in all output files is 0-based (BED format standard). Even if the input is in SAF format (1-based), the output is converted to 0-based.

## Citation

If you use this tool in your research, please cite:
**(Paper in preparation)**

## License

This program is released under the MIT license. For details, see [LICENSE](./LICENSE).