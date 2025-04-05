# SEgene_peakprep

**SEgene_peakprep** is a Python pipeline for normalizing read counts from ChIP-seq data (BAM files) in specified genomic regions and calculating Counts Per Million (CPM) values. This tool is designed for data mapped to the **hg38 reference genome**.

## Overview

SEgene_peakprep provides the following functionalities:

- Efficient processing of multiple BAM files with automatic sample information extraction
- Calculation of total mapped reads using `samtools flagstat`
- Read counting in specified genomic regions (SAF format) using `featureCounts`
- CPM value calculation with optional log transformation (log2, log10)
- Formatting and saving of results

## Development Status

> **Note**: This tool is **currently under development**, and the program file structure and usage methods may change significantly between versions.

**Position within the SEgene project:**

This pipeline functions as the initial data preparation module within the [SEgene project](https://github.com/hamamoto-lab/SEgene) framework:

1. **Input**: Mapped ChIP-seq data (BAM files) and regions of interest (SAF format)
2. **Output**: CPM/logCPM value table (TSV file)
3. **Flow**: This output serves as a key input for the `peak_to_gene_links` module, where it is combined with RNA-seq TPM data for gene regulatory network analysis

## Requirements

- **Python**: Python 3.7 or higher
- **Required libraries**:
  - pandas
  - numpy
  - natsort
  
  *Note: Please install the required libraries using conda or pip*
  
- **External tools**:
  - samtools (for flagstat functionality) - Ensure it's in your PATH or specify the path at runtime
  - featureCounts (part of the Subread package) - Ensure it's in your PATH or specify the path at runtime

## Installation

To use this pipeline, install the required Python libraries (pandas, numpy, natsort) using conda or pip.

Clone this repository or download the scripts:

```bash
git clone https://github.com/hamamoto-lab/SEgene.git
cd SEgene/SEgene_peakprep
```

## Input Data Preparation

### SAF File Preparation
This pipeline uses Simplified Annotation Format (SAF) files as input. SAF files are tab-delimited text files containing the following columns:

1. `GeneID` - Unique identifier for each region
2. `Chr` - Chromosome name
3. `Start` - Start position (1-based)
4. `End` - End position
5. `Strand` - Strand (+, -, or .)

Example:
```
GeneID  Chr     Start   End     Strand
peak1   chr1    1000    2000    .
peak2   chr1    3000    4000    .
```

If you need to convert BED files to SAF format, you can use a script like the following:

```python
import pandas as pd

# Read BED file
bed_df = pd.read_csv("your_file.bed", sep='\t', header=None)

# Extract necessary columns and convert to SAF format
saf_df = pd.DataFrame({
    'GeneID': [f"peak_{i+1}" for i in range(len(bed_df))],
    'Chr': bed_df[0],  # Chromosome column
    'Start': bed_df[1] + 1,  # Add 1 because BED is 0-based
    'End': bed_df[2],
    'Strand': '.'  # Use '.' if strand information is unavailable
})

# Save as SAF file
saf_df.to_csv("output.saf", sep='\t', index=False)
```

### BAM File Preparation
BAM files should be generated using standard ChIP-seq mapping pipelines (such as nf-core/chipseq). It is recommended that BAM files be indexed (with .bai files).

## Usage

The pipeline can be executed from the command line or from a Jupyter notebook.

### Command Line Execution

Basic execution command:

```bash
python cpm_peakprep_pipeline.py \
    --bam_dir /path/to/bam_files \
    --saf_file peaks.saf \
    --output_dir results \
    --filename_pattern "Sample*.sorted.bam" \
    --sample_delimiter ".sorted.bam" \
    --threads 10
```

#### Argument Description

| Argument | Description |
|------|------|
| `--bam_dir` | Directory containing BAM files (**required**) |
| `--saf_file` | Path to the annotation file in SAF format (**required**) |
| `--output_dir` | Directory to save results (**required**) |
| `--filename_pattern` | Wildcard pattern to filter BAM filenames (e.g., `*.bam`, `T*H3K27ac*.bam`) (optional, default: `*.bam`) |
| `--sample_delimiter` | Delimiter string to extract sample names from BAM filenames (optional)<br>Example: If filename is `T9-H3K27ac_R1.mLb.clN.sorted.bam` and delimiter is `-H3K27ac_R1.mLb.clN.sorted.bam`, the sample name will be `T9` |
| `--samtools_path` | Path to samtools executable (optional, default: `samtools`) |
| `--featurecounts_path` | Path to featureCounts executable (optional, default: `featureCounts`) |
| `--threads` | Number of threads for samtools and featureCounts (optional, default: 4) |
| `--single_end` | Flag to specify single-end reads. If not specified, paired-end is assumed (optional, default: paired-end) |
| `--fc_output_name` | featureCounts output filename (optional, default: automatically generated from SAF filename) |
| `--fc_log_name` | featureCounts log filename (optional) |
| `--fc_options` | Additional options to pass to featureCounts (optional) |
| `--add_peakid` | Flag to replace original GeneIDs with sequential PeakIDs (PeakID_XXXX) (optional) |
| `--id_column` | Name of the ID column in featureCounts output (optional, default: `Geneid`) |
| `--output_id_name` | Name of the ID column in the final logCPM table (optional, default: `PeakID`) |
| `--log_base` | Base for log transformation (2 or 10). Set to 0 or negative for no transformation (optional, default: 2) |
| `--pseudocount` | Pseudocount added before log transformation (optional, default: 1.0) |
| `--logcpm_output_name` | Filename for the final logCPM table (optional, default: `logCPM_table.tsv`) |
| `--log_level` | Logging level (DEBUG/INFO/WARNING/ERROR/CRITICAL) (optional, default: INFO) |

### Jupyter Notebook Execution

Here is a basic example of executing the pipeline in a Jupyter notebook:

```python
import os
import logging
import pandas as pd
# Import utility functions
from cpm_peakprep_utils import (
    get_sample_data,
    run_samtools_flagstat,
    run_featurecounts,
    calculate_logcpm
)

# Parameter settings
bam_dir = "/path/to/bam_files"
saf_file = "peaks.saf"
output_dir = "results"
filename_pattern = "Sample*.sorted.bam"
sample_delimiter = ".sorted.bam"
threads = 4

# Create output directory
os.makedirs(output_dir, exist_ok=True)

# Logger configuration
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger("ChIPseqPipeline")

# Step 1: Get BAM file information
sample_dict, bam_list = get_sample_data(
    bam_folder=bam_dir,
    logger=logger,
    sample_name_delimiter=sample_delimiter,
    filename_pattern=filename_pattern
)

# Step 2: Get total read counts
flagstat_dir = os.path.join(output_dir, "flagstat_output")
total_reads_map = run_samtools_flagstat(
    bam_files=bam_list,
    output_dir=flagstat_dir,
    logger=logger,
    threads=threads
)

# Step 3: Count reads (featureCounts)
fc_output_path = os.path.join(output_dir, f"{os.path.splitext(os.path.basename(saf_file))[0]}_featureCounts.txt")
fc_log_path = os.path.join(output_dir, f"{os.path.splitext(os.path.basename(saf_file))[0]}_featureCounts.log")
counts_file_path = run_featurecounts(
    bam_files=bam_list,
    saf_file=saf_file,
    logger=logger,
    output_file=fc_output_path,
    log_file=fc_log_path,
    threads=threads,
    is_paired_end=True
)

# Step 4: LogCPM calculation
bam_to_sample_map = {v: k for k, v in sample_dict.items()}
logcpm_df = calculate_logcpm(
    counts_file=counts_file_path,
    total_reads_dict=total_reads_map,
    logger=logger,
    bam_path_to_sample_name=bam_to_sample_map,
    log_transform_base=2,
    pseudocount=1.0
)

# Save results
logcpm_output_file = os.path.join(output_dir, "logCPM_table.tsv")
logcpm_df.to_csv(logcpm_output_file, sep='\t', index=False, float_format='%.4f')
print(f"Final logCPM table saved to: {logcpm_output_file}")
```

## Output Format

The pipeline generates the following files:

1. **flagstat output** (`flagstat_output/` directory):
   - samtools flagstat results for each BAM file (including total mapped reads information)

2. **featureCounts output** (e.g., `peaks_featureCounts.txt`):
   - Raw read counts for each region
   - Columns: GeneID, Chr, Start, End, Strand, Length, counts for each BAM file

3. **logCPM table** (e.g., `logCPM_table.tsv`):
   - Final log-transformed CPM values
   - Columns: PeakID, Chr, Start, End, logCPM values for each sample

Example (logCPM table):
```
PeakID  Chr     Start   End     T1      T2      T3
peak1   chr1    1000    2000    2.45    3.12    1.87
peak2   chr1    3000    4000    4.21    3.98    4.56
...
```

## Pipeline Workflow

1. **BAM File Information Retrieval**:
   - Search for BAM files matching criteria in the specified directory
   - Extract sample names from filenames

2. **Total Mapped Reads Retrieval**:
   - Calculate total mapped reads for each BAM file using `samtools flagstat`
   - Results are used for CPM normalization

3. **Counting Reads by Region**:
   - Count reads in specified regions (SAF file) from each BAM file using `featureCounts`
   - Compatible with paired-end or single-end reads

4. **CPM Calculation and Log Transformation**:
   - Normalize each read count by total mapped reads (CPM)
   - Optionally apply log transformation (log2 or log10)
   - Convert sample names to user-friendly format

5. **Result Output**:
   - Save the final table in TSV format

## Position in SEgene Workflow

This pipeline functions as part of the following workflow in the SEgene project:

1. **ChIP-seq data preprocessing**: Generate BAM files using standard mapping pipelines (e.g., nf-core/chipseq)
2. **CPM/logCPM value calculation** (this pipeline)
3. **Construction of peak-to-gene links** (`peak_to_gene_links` module)
   - Input: logCPM file generated by this pipeline
   - Input: RNA-seq TPM data
4. **Super-enhancer identification and analysis** (`SE_to_gene_links`)
5. **Optional: Public data-based region evaluation** (`SEgene_RegionAnalyzer`)

## Troubleshooting

Common issues and solutions:

1. **samtools/featureCounts not found**:
   - Verify they are added to your PATH, or specify the complete path with `--samtools_path` and `--featurecounts_path`

2. **BAM files not found**:
   - Check the `--filename_pattern` argument, verify the pattern is appropriate

3. **featureCounts returns an error**:
   - Check if BAM files are indexed
   - Verify paired-end/single-end settings are correct (`--single_end` flag)
   - Check if the SAF file format is correct

## License

This program is released under the MIT license. For details, please see [LICENSE](./LICENSE).

## Citation

If you use this tool in your research, please cite:
**(Paper currently in preparation.)**