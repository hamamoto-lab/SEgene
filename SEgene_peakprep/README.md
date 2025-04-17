# SEgene_peakprep

*(For the Japanese version of this README, please see [README_ja.md](./README_ja.md).)*

**SEgene_peakprep** is a Python pipeline for quantifying and normalizing signal values in specified genomic regions from sequence data (BAM files), and formatting them for downstream SEgene analysis. This tool is not dependent on any specific reference genome and can be used with data mapped to genomes such as hg38.

## Development Status

> **Note:** This tool is **currently under development**, and the program file structure and usage may change between versions.

## Overview

SEgene_peakprep serves as the initial data preparation step in the [SEgene project](https://github.com/hamamoto-lab/SEgene):

1. **Input**: Mapped sequence data (BAM files) and annotation files (BED files, etc.)
2. **Output**: Normalized count value tables
3. **Flow**: This output serves as the primary input for the `peak_to_gene_links` module, which is combined with RNA-seq TPM data for gene regulatory network analysis

This pipeline offers two implementation methods:

1. **CPM Method**: Uses `featureCounts` to directly calculate Counts Per Million (CPM) values
2. **BigWig Method**: Uses `deeptools` to convert BAM to bigWig files, then quantifies peaks using `multiBigwigSummary`

Both methods support BED, SAF, and merge_SE format annotation files. You can choose the appropriate method based on your needs.  
**Note:** The merge_SE format is specific to the SEgene project and is used as output and intermediate files in the pipeline.

## Common Requirements and Dependencies

The following common requirements must be met:

### Python Environment
- Python 3.10 or higher

### Python Libraries (install with pip or conda)
| Library    | Purpose                           | Recommended Version | Required For             |
|------------|-----------------------------------|--------------------:|-------------------------|
| pandas     | Table manipulation                | ≥1.5                | All functionality        |
| numpy      | Numerical computation             | ≥1.23               | All functionality        |
| natsort    | Natural sorting                   | ≥8.3                | All functionality        |
| pyranges   | Genomic region operations         | ≥0.14               | BigWig method only      |

```bash
# Installation with pip
pip install pandas numpy natsort
pip install pyranges  # Required for BigWig method
```

```bash
# Installation with conda
conda install pandas numpy natsort
conda install -c bioconda pyranges  # Required for BigWig method
```

### External Tools (install with conda or binaries)
| Tool Name           | Purpose                          | Recommended Version | Required For            |
|---------------------|----------------------------------|--------------------:|------------------------|
| samtools            | BAM manipulation                 | ≥1.13               | CPM method             |
| featureCounts       | CPM counting                     | ≥2.0.0              | CPM method             |
| deeptools           | BigWig generation & diff analysis| ≥3.5.0              | BigWig method          |
| - bamCoverage       | BAM to bigWig conversion         | (included in deeptools) | BigWig single sample analysis |
| - multiBigwigSummary| bigWig signal aggregation        | (included in deeptools) | BigWig single sample analysis |
| - bamCompare        | Sample/control comparison        | (included in deeptools) | BigWig differential analysis pipeline |

```bash
# Installation with conda
conda install -c bioconda samtools subread  # CPM method
conda install -c bioconda deeptools         # BigWig method
```

### Common Input Files

Both implementation methods require the following input files:

1. **BAM Files**: Mapped sequence data files (recommended: indexed .bai files)
2. **Annotation Files**: Files that define genomic regions, supporting these formats:
   - **BED Format**: Standard BED format (0-based start coordinates)
     ```
     # Example: peaks.bed (tab-delimited)
     chr1    1000    2000    peak1    0    +
     chr1    3000    4000    peak2    0    +
     chr2    5000    6000    peak3    0    -
     ```
   
   - **SAF Format** (Simplified Annotation Format): Tab-delimited text file (**1-based start coordinates**)
     ```
     # Example: peaks.saf (tab-delimited)
     GeneID  Chr     Start   End     Strand
     peak1   chr1    1001    2000    +
     peak2   chr1    3001    4000    +
     peak3   chr2    5001    6000    -
     ```
     **Note:** SAF format coordinates are 1-based (Start=1001 is 1000 in 0-based), but all outputs are unified to 0-based
   
   - **merge_SE.tsv Format**: First column is region ID in the format `chr_start_end` (e.g., `chr1_1000_2000`)
     ```
     # Example: merged_se.tsv (tab-delimited)
     se_data            feature
     chr1_1000_2000     feature1
     chr1_3000_4000     feature2
     chr2_5000_6000     feature3
     ```
     **Notes:** 
      - In merge_SE format, coordinates are interpreted as 0-based
      - The merge_SE format is specific to the SEgene project and is used as output and intermediate files in the SEgene pipeline
      - Coordinates are 0-based (chr1_1000_2000's Start is 1000 in 0-based)
      - This format allows seamless integration with other SEgene components

## CPM Method (cpm_peakprep)

The CPM method uses `featureCounts` to directly count reads in specified genomic regions (BED/SAF/merge_SE.tsv format) and calculates CPM values and their log-transformed values.

### Basic Usage

```bash
python cpm_peakprep.py \
    --bam_dir /path/to/bam_files \
    --annotation_file peaks.bed \
    --output_dir results \
    --filename_pattern "Sample*.sorted.bam" \
    --sample_delimiter=".sorted.bam" \
    --threads 10
```

### Key Arguments

| Argument | Description | Default Value |
|--------|------|-------------|
| **Required Arguments** |||
| `--bam_dir` / `-b` | Directory containing BAM files | - |
| `--annotation_file` / `-a` | Annotation file (BED, SAF, merge_SE format) | - |
| `--output_dir` / `-o` | Directory to save results | - |
| **Annotation Related** |||
| `--is_mergese_format` | Specify if the annotation file is in merge_SE.tsv format | False |
| **File/Sample Selection** |||
| `--filename_pattern` | Wildcard pattern for BAM filenames | "*.bam" |
| `--sample_delimiter` | Delimiter string for extracting sample names | None |
| **Tool Paths & Threads** |||
| `--samtools_path` | Path to samtools executable | "samtools" |
| `--featurecounts_path` | Path to featureCounts executable | "featureCounts" |
| `--threads` / `-T` | Number of threads | 4 |
| `--single_end` | Run in single-end mode (default is paired-end) | False |

To see all options, run:
```bash
python cpm_peakprep.py --help
```

For more detailed configuration options and output explanations (including Jupyter notebook examples), refer to the [CPM Method Detailed Documentation](./cpm_README.md).

## BigWig Method (bigwig_peakprep)

The BigWig method uses `deeptools` to process BAM files and quantify signal values in genomic regions. It provides **two types of pipelines**: **single sample analysis** and **sample-control differential analysis**.

### Pipeline Structure

The BigWig method offers two main execution paths:

1. **Single Sample Analysis Pipeline** (`bigwig_peakprep.py`)
   - Generates normalized bigWig files from BAM files and quantifies signal values in specified regions
   - Internally executes the following scripts in sequence:
     1. `bigwig_peakprep_bamconvert.py`: Converts BAM files to bigWig files
     2. `bigwig_peakprep_summary.py`: Generates count tables for peak regions from bigWig files

2. **Differential Analysis Pipeline** (`bigwig_peakprep_bamdiff.py`) - **New Feature**
   - Performs log2 ratio analysis of ChIP-seq samples and Input controls
   - Internally executes the following scripts in sequence:
     1. `bigwig_peakprep_bamdiff_generatediff.py`: Generates log2 ratio bigWig files from sample and control
     2. `bigwig_peakprep_summary.py`: Generates region count tables from the differential bigWig files

**Utility Modules**:
- `bigwig_peakprep_utils.py`: Provides basic functions and common utilities
- `bigwig_peakprep_bamdiff_utils.py`: Provides specialized functions for differential analysis

### Which Pipeline to Choose

- Use the **Single Sample Analysis Pipeline** (`bigwig_peakprep.py`) when:
  - You want to analyze and compare each sample independently
  - You need normalized signal values from individual BAM files

- Use the **Differential Analysis Pipeline** (`bigwig_peakprep_bamdiff.py`) when:
  - You want to compare ChIP-seq and Input controls, removing background noise
  - You have sample and control pairs (e.g., T1-H3K27ac and T1-Input)

By default, bamCoverage uses **RPGC (Reads Per Genomic Content)** normalization, but you can change to other normalization methods (CPM, BPM, RPKM, RPGC, None) with the `--normalize_using` option.

### Basic Usage

#### For Single Sample Analysis (bigwig_peakprep.py)

```bash
python bigwig_peakprep.py \
    --bam_dir="/path/to/bam_files" \
    --annotation_file="/path/to/peaks.bed" \
    --output_dir="/path/to/output" \
    --tools_dir="/path/to/deeptools/bin" \
    --sample_delimiter="-H3K27ac.sorted.bam" \
    --threads=10
```

#### For Sample-Control Differential Analysis (bigwig_peakprep_bamdiff.py)

```bash
python bigwig_peakprep_bamdiff.py \
    --sample_bam_dir="/path/to/chipseq_bams" \
    --control_bam_dir="/path/to/input_bams" \
    --annotation_file="/path/to/peaks.bed" \
    --output_dir="/path/to/output" \
    --tools_dir="/path/to/deeptools/bin" \
    --sample_delimiter="-H3K27ac.sorted.bam" \
    --threads=10
```

### Key Arguments

#### Single Sample Analysis (bigwig_peakprep.py)

| Argument | Description | Default Value |
|--------|------|-------------|
| **Required Arguments** |||
| `--bam_dir` | Directory containing BAM files | - |
| `--annotation_file` | Annotation file (BED, SAF, merge_SE format) | - |
| `--output_dir` | Directory to save results | - |
| **Tool Paths** |||
| `--tools_dir` | Directory containing deeptools executables | None |
| `--bamcoverage_path` | Path to bamCoverage executable | "bamCoverage" |
| `--multibigwigsummary_path` | Path to multiBigwigSummary executable | "multiBigwigSummary" |
| **Execution Step Control** |||
| `--run_bamconvert_only` | Run only BAM to bigWig conversion | False |
| `--run_summary_only` | Run only bigWig summary processing | False |
| **Normalization Related** |||
| `--normalize_using` | bamCoverage normalization method (RPGC, CPM, BPM, RPKM, None) | "RPGC" |
| `--effective_genome_size` | Effective genome size for RPGC normalization | 2913022398 |

#### Sample-Control Differential Analysis (bigwig_peakprep_bamdiff.py)

| Argument | Description | Default Value |
|--------|------|-------------|
| **Required Arguments** |||
| `--sample_bam_dir` | Directory containing ChIP-seq sample BAM files | - |
| `--control_bam_dir` | Directory containing Input control (Input/IgG) BAM files | - |
| `--annotation_file` | Annotation file (BED, SAF, merge_SE format) | - |
| `--output_dir` | Directory to save results | - |
| **File Selection** |||
| `--sample_pattern` | Glob pattern for sample BAM files | "*.bam" |
| `--control_pattern` | Glob pattern for control BAM files | "\*Input\*.bam" |
| **Processing Related** |||
| `--operation` | bamCompare operation method (log2, ratio, subtract, add, mean, etc.) | "log2" |
| `--bin_size` | bamCompare bin size (in bases) | 50 |
| `--pseudocount` | Pseudocount value to avoid division by zero | 1.0 |
| `--run_diff_only` | Run only differential bigWig generation (skip summary processing) | False |

To see all options, run:
```bash
python bigwig_peakprep.py --help
python bigwig_peakprep_bamdiff.py --help
```

For more detailed configuration options and output explanations, refer to the [BigWig Method Detailed Documentation](./bigwig_README.md).

## Output Formats

> **Important Note:** Even if SAF format (1-based coordinates) is used as the input annotation file, **all output table coordinates (Start column) are unified to 0-based (BED format).**

### CPM Method Output

1. **logCPM Table** (default: `logCPM.tsv`):
   - Table containing log-transformed CPM values for each peak/region
   - Columns: PeakID, Chr, Start (0-based, BED format), End, logCPM values for each sample

2. **Intermediate Files**:
   - flagstat output files (total mapped read counts)
   - featureCounts output files (raw read counts)

### BigWig Method Output

#### For Single Sample Analysis (bigwig_peakprep.py)

1. **bigWig Files** (Step 1 output):
   - Normalized bigWig files generated from each BAM file
   - Filename format: `{sample_name}.{normalization_method}.bigwig`
   - Example: `sample1.RPGC.bigwig`

2. **bam_to_bigwig_details.tsv**:
   - Step 1 detailed output file (input for Step 2)
   - Comment lines contain parameter information and execution metadata
   - Columns: `Sample_name`, `BAM_filename`, `BigWig_filename`, `BAM_fullpath`, `BigWig_fullpath`

3. **Count Tables** (Step 2 output):
   - `multibigwig_counts.tsv`: Raw count table
   - `multibigwig_counts_log2.tsv`: log2-transformed count table

#### For Sample-Control Differential Analysis (bigwig_peakprep_bamdiff.py)

1. **Differential bigWig Files**:
   - bigWig files representing log2 ratios between sample and control
   - Filename format: `{sample_name}_vs_{control_name}.log2ratio.bigwig`
   - Example: `T1-H3K27ac_vs_T1-Input.log2ratio.bigwig`

2. **bamdiff_to_bigwig_details.tsv**:
   - TSV file containing details of the differential analysis
   - Columns: `Sample_name`, `BigWig_filename`, `BigWig_fullpath`
   - Important: This file includes `already_log2_transformed: true` metadata to prevent additional log2 transformation in summary processing

3. **Conversion Plan**:
   - `conversion_plan.yaml`: Execution plan for differential analysis
   - `conversion_plan.tsv`: TSV format of the conversion plan (tabular format)

4. **Count Tables**:
   - `multibigwig_counts.tsv`: Table containing log2 ratio values for each peak region
   - Note: No additional log2 transformation is performed as the values are already log2-transformed

Example output table format:
```
PeakID  Chr     Start   End     Sample1  Sample2  Sample3
peak1   chr1    1000    2000    2.45     3.12     1.87
peak2   chr1    3000    4000    4.21     3.98     4.56
...
```

### Example Output Directory Structure

#### Single Sample Analysis Output Directory Structure:
```
output_dir/
├── bigwig/
│   ├── sample1.RPGC.bigwig
│   ├── sample2.RPGC.bigwig
│   └── ...
├── summary/
│   ├── bigwig_summary.npz
│   └── bigwig_summary.tab
├── bam_to_bigwig_details.tsv
├── multibigwig_counts.tsv
├── multibigwig_counts_log2.tsv
└── pipeline_run.log
```

#### Differential Analysis Output Directory Structure:
```
output_dir/
├── bigwig/
│   ├── T1-H3K27ac_vs_T1-Input.log2ratio.bigwig
│   ├── T2-H3K27ac_vs_T2-Input.log2ratio.bigwig
│   └── ...
├── summary/
│   ├── bigwig_summary.npz
│   └── bigwig_summary.tab
├── bamdiff_to_bigwig_details.tsv
├── bamdiff_to_bigwig_details.yaml
├── conversion_plan.yaml
├── conversion_plan.tsv
├── multibigwig_counts.tsv
└── pipeline_run.log
```

## Troubleshooting

Common issues and solutions:

1. **Tool not found error**:
   - Verify that samtools, featureCounts, or deeptools are installed
   - Specify the full path using the corresponding option
   ```bash
   --tools_dir="/path/to/deeptools/bin"
   --samtools_path="/path/to/samtools"
   --featurecounts_path="/path/to/featureCounts"
   ```

2. **BAM files not found**:
   - Check the `--filename_pattern` argument, `--sample_pattern`, `--control_pattern`
   - Verify file permissions and read access

3. **Sample name extraction issues**:
   - Set the `--sample_delimiter` value appropriately and use an equals sign (`=`)
   - For values with hyphens, always use the equals sign: `--sample_delimiter="-H3K27ac.sorted.bam"`

4. **Out of memory error**:
   - Run on a machine with more memory
   - Reduce the number of files processed or process them in batches

5. **Slow processing**:
   - Increase the number of threads with the `--threads` option
   - Use indexed BAM files
   - Increase the `--bin_size` value for faster processing, though at the cost of lower resolution

6. **Sample-control matching failure**:
   - If you see warning message "Could not create sample-control pairs" in the log, matching has failed
   - Explicitly specify file name patterns with `--sample_pattern` and `--control_pattern`
   - Ensure file names have common identifiers (T1, T2, etc.)
   - For using a single control for all samples, ensure the control directory contains only one file

## Citation

If you use this tool in your research, please cite:
**(Paper in preparation)**

## License

This program is released under the MIT license. For details, see [LICENSE](./LICENSE)