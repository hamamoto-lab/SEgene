# SEgene_peakprep

*(For the Japanese version of this README, please see [README_ja.md](./README_ja.md).)*

**SEgene_peakprep** is a Python pipeline for normalizing read counts in specified genomic regions from sequence data (BAM files) and calculating CPM (Counts Per Million) values. This tool is not dependent on any specific reference genome and can be used with data mapped to any genome, such as hg38.

## Development Status

> **Note**: This tool is **currently under development**, and the program structure and usage may change significantly between versions.

## Overview

SEgene_peakprep provides the following features:

- Efficient processing of multiple BAM files and automatic extraction of sample information
- Calculation of total mapped reads using `samtools flagstat`
- Counting of reads in specified genomic regions (BED, SAF, or merge_SV.tsv format) using `featureCounts`
- Calculation of CPM values and log transformation (default is log2(CPM+1))
- Formatting and saving of results

**Position in the SEgene Project:**

This pipeline serves as the initial data preparation step in the [SEgene project](https://github.com/hamamoto-lab/SEgene):

1. **Input**: Mapped sequence data (BAM files) and annotation files (BED/SAF/merge_SV.tsv format)
2. **Output**: CPM/logCPM value tables (TSV files)
3. **Flow**: This output serves as the main input for the `peak_to_gene_links` module, which is combined with RNA-seq TPM data for gene regulatory network analysis

## Requirements

- **Python**: Python 3.7 or higher
- **Required libraries**:
  - pandas
  - numpy
  - natsort
  
  *These libraries can be installed using conda or pip:*
  ```bash
  # Using conda
  conda install pandas numpy natsort
  
  # Using pip
  pip install pandas numpy natsort
  ```
  
- **External tools**:
  - samtools (recommended: v1.13 or higher)
  - featureCounts (part of the Subread package, recommended: v2.0.0 or higher)
  
  *These tools should be in your PATH or specified via options at runtime:*
  ```bash
  # Installation example using conda
  conda install -c bioconda samtools subread
  ```
  
  For detailed installation instructions, refer to the official sites for [samtools](http://www.htslib.org/) and [Subread/featureCounts](http://subread.sourceforge.net/).

## Installation

```bash
# Clone the repository
git clone https://github.com/hamamoto-lab/SEgene.git
cd SEgene/SEgene_peakprep
```

## Preparing Input Data

### Preparing Annotation Files
This pipeline accepts the following annotation file formats as input:

1. **BED format**: Standard BED format (0-based start coordinates). The pipeline will automatically convert it to SAF format internally.

2. **SAF format** (Simplified Annotation Format): A tab-delimited text file containing the following columns:
   - `GeneID` - A unique identifier for each region
   - `Chr` - Chromosome name
   - `Start` - Start position (1-based)
   - `End` - End position
   - `Strand` - Strand (+, -, or .)

   Example:
   ```
   GeneID  Chr     Start   End     Strand
   peak1   chr1    1000    2000    .
   peak2   chr1    3000    4000    .
   ```

3. **merge_SV.tsv format**: A tab-delimited file with region IDs in the first column in the format `chr_start_end` (e.g., `chr1_1000_2000`). When using this format, you must specify the `--is_mergesv_format` flag.

The annotation file is specified using the `--annotation_file` argument. When a BED file or merge_SV file is specified, it will be automatically converted to SAF format at runtime.

### Preparing BAM Files
BAM files should be generated using standard sequence mapping pipelines. It is recommended that the BAM files be indexed (with .bai files).

## Usage

The pipeline can be run from the command line or a Jupyter notebook.

### Command Line Execution

**Execution command using a BED file (Basic usage):**

```bash
python cpm_peakprep.py \
    --bam_dir /path/to/bam_files \
    --annotation_file peaks.bed \
    --output_dir results \
    --filename_pattern "Sample*.sorted.bam" \
    --sample_delimiter=".sorted.bam" \
    --threads 10
```

**Execution command using a SAF file:**

```bash
python cpm_peakprep.py \
    --bam_dir /path/to/bam_files \
    --annotation_file peaks.saf \
    --output_dir results \
    --filename_pattern "Sample*.sorted.bam" \
    --sample_delimiter=".sorted.bam" \
    --threads 10
```

**Execution command using a merge_SV.tsv format file:**

```bash
python cpm_peakprep.py \
    --bam_dir /path/to/bam_files \
    --annotation_file merged_se.tsv \
    --is_mergesv_format \
    --output_dir results \
    --filename_pattern "Sample*.sorted.bam" \
    --sample_delimiter=".sorted.bam" \
    --threads 10
```

> **Note**: When an argument value starts with a hyphen, you must use an equals sign (`=`) to specify it (e.g., `--sample_delimiter="-H3K27ac.sorted.bam"`).

#### Argument Descriptions

**Required Arguments:**

| Argument Name           | Short Form | Description                                           |
| :---------------------- | :--------- | :---------------------------------------------------- |
| `--bam_dir`             | `-b`       | Directory containing BAM files.                       |
| `--annotation_file`     | `-a`       | Path to the annotation file (BED, SAF, or merge_SV.tsv format). |
| `--output_dir`          | `-o`       | Directory to save results.                            |

**Annotation Options:**

| Argument Name           | Short Form | Default Value | Description                                           |
| :---------------------- | :--------- | :------------ | :---------------------------------------------------- |
| `--is_mergesv_format`   |            | `False`       | Flag to specify that the annotation file is in merge_SV.tsv format. |

**File/Sample Selection Options:**

| Argument Name           | Short Form | Default Value | Description                                           |
| :---------------------- | :--------- | :------------ | :---------------------------------------------------- |
| `--filename_pattern`    |            | `*.bam`       | Wildcard pattern for filtering BAM filenames.         |
| `--sample_delimiter`    |            | `None`        | Delimiter string to extract sample names from BAM filenames. |

**Tool Paths and Thread Options:**

| Argument Name           | Short Form | Default Value | Description                                           |
| :---------------------- | :--------- | :------------ | :---------------------------------------------------- |
| `--samtools_path`       |            | `samtools`    | Path to the samtools executable.                      |
| `--featurecounts_path`  |            | `featureCounts` | Path to the featureCounts executable.                |
| `--threads`             | `-T`       | `4`           | Number of threads to use for samtools and featureCounts. |
| `--single_end`          |            | `False`       | Flag to specify single-end reads (default is paired-end). |

**featureCounts Options:**

| Argument Name           | Short Form | Default Value | Description                                           |
| :---------------------- | :--------- | :------------ | :---------------------------------------------------- |
| `--fc_basename`         |            | `None`        | Base name for featureCounts output/log files (default: `ANNOTATION_fc`). |
| `--fc_options`          |            | `None`        | Additional options for featureCounts (e.g., `--fc_options --minOverlap 10`). |

**LogCPM Calculation Options:**

| Argument Name           | Short Form | Default Value | Description                                           |
| :---------------------- | :--------- | :------------ | :---------------------------------------------------- |
| `--add_peakid`          |            | `False`       | Flag to replace original GeneIDs with sequential PeakIDs. |
| `--id_column`           |            | `Geneid`      | ID column name in SAF/featureCounts input.            |
| `--output_id_name`      |            | `PeakID`      | ID column name in the final logCPM table.             |
| `--log_base`            |            | `2`           | Log base (2 or 10). Set to ≤0 for no transformation. |
| `--pseudocount`         |            | `1.0`         | Pseudocount to add before log transformation.         |

**Output File Options:**

| Argument Name           | Short Form | Default Value | Description                                           |
| :---------------------- | :--------- | :------------ | :---------------------------------------------------- |
| `--logcpm_output_name`  |            | `logCPM.tsv`  | Filename for the logCPM table.                        |

**Logging Options:**

| Argument Name           | Short Form | Default Value | Description                                           |
| :---------------------- | :--------- | :------------ | :---------------------------------------------------- |
| `--log_level`           |            | `INFO`        | Logging level (DEBUG/INFO/WARNING/ERROR/CRITICAL).    |
| `--script_log_file`     |            | `pipeline_run.log` | Filename for script execution log (empty string to disable). |

*(Default values might change. Check `python cpm_peakprep.py --help` for the latest information)*

### Execution in Jupyter Notebook

Here's a basic example of running the pipeline in a Jupyter notebook:

```python
import os
import sys
import logging
import pandas as pd

# Set up path to import modules
sys.path.append('/path/to/SEgene/SEgene_peakprep')

# Import specific functions
from cpm_peakprep_utils import (
    get_sample_data,
    run_samtools_flagstat,
    run_featurecounts,
    calculate_logcpm,
    convert_bed_to_saf
)

# Parameter settings
bam_dir = "/path/to/bam_files"
annotation_file = "peaks.bed"  # BED, SAF, or merge_SV.tsv file
output_dir = "results"
filename_pattern = "Sample*.sorted.bam"
sample_delimiter = ".sorted.bam"
threads = 10

# Create output directory
os.makedirs(output_dir, exist_ok=True)

# Set up logger
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger("PipelineLogger")

# Convert BED file to SAF if needed
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

# Step 3: Count reads using featureCounts
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

# Read counts file and Step 4: Calculate LogCPM
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
print(f"Saved logCPM table to: {logcpm_output_file}")

# Clean up temporary files
if 'saf_file' in locals() and saf_file != annotation_file and os.path.exists(saf_file):
    try:
        os.remove(saf_file)
        print(f"Removed temporary SAF file: {saf_file}")
    except Exception as e:
        print(f"Failed to delete temporary file: {e}")
```

## Output Format

The pipeline generates the following files:

1. **logCPM table** (default: `logCPM.tsv`):
   - Log-transformed CPM values for each peak/region in samples
   - Columns: PeakID, Chr, Start (0-based, BED format standard), End, logCPM values for each sample

2. **Intermediate files**:
   - **flagstat output** (`flagstat/` directory):
     - samtools flagstat results for each BAM file (total mapped reads and other statistics)
   - **featureCounts output** (e.g., `annotation_base_fc_featureCounts.txt`):
     - Raw read counts for each region
   - **Execution log** (e.g., `pipeline_run.log`):
     - Detailed log of pipeline execution

Example (logCPM table):
```
PeakID  Chr     Start   End     Sample1  Sample2  Sample3
peak1   chr1    999     2000    2.45     3.12     1.87
peak2   chr1    2999    4000    4.21     3.98     4.56
...
```

> **Note**: The `Start` column in all output files uses 0-based coordinates (BED format standard). Even if the input is in SAF format (1-based), the output will be converted to 0-based.

## Position in the SEgene Workflow

This pipeline functions as part of the SEgene project workflow as follows:

1. **Sequence data preprocessing**: Generate BAM files using standard mapping pipelines
2. **CPM/logCPM value calculation** (this pipeline)
3. **peak-to-gene links construction** (`peak_to_gene_links` module)
   - Input: logCPM files generated by this pipeline
   - Input: RNA-seq TPM data
4. **Super-enhancer identification and analysis** (`SE_to_gene_links`)
5. **Optional: Region evaluation based on public data** (`SEgene_RegionAnalyzer`)

## Troubleshooting

Common issues and solutions:

1. **samtools/featureCounts not found**:
   - Check if they are in your PATH or specify the full path using `--samtools_path` and `--featurecounts_path`

2. **BAM files not found**:
   - Check the `--filename_pattern` argument

3. **Sample name extraction issues**:
   - Verify the `--sample_delimiter` value
   - Use the format `--delimiter="value"` when the delimiter contains hyphens

4. **Annotation file format issues**:
   - For BED files, make sure the extension is `.bed`
   - For merge_SV format, specify the `--is_mergesv_format` flag
   - For SAF files, ensure the format is correct (tab-delimited)

5. **featureCounts returns errors**:
   - Check if BAM files are indexed
   - Verify paired-end/single-end settings (`--single_end` flag)

## Citation

When using this tool for research, please cite:
**(Paper is currently in preparation)**

## License

This program is distributed under the MIT license. For details, see [LICENSE](./LICENSE).