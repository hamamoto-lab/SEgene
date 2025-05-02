# SEgene_peakprep - BigWig Method Detailed Guide

This document provides detailed usage instructions and examples specific to the **BigWig method** of SEgene_peakprep. For basic overview and common requirements, please refer to the [main README](./README.md).

## Pipeline Structure and Processing Flow

The BigWig method consists of two main pipelines that can be used depending on your processing needs:

### 1. Single Sample Analysis Pipeline

**Main script**: `bigwig_peakprep.py`

This pipeline generates normalized bigWig files from individual BAM files and quantifies signal values in specific regions (such as peaks). Internally, it processes data through the following steps:

1. `bigwig_peakprep_bamconvert.py`: Converts BAM files to normalized bigWig files
   - Uses bamCoverage to generate bigWig files for each sample
   - Output: bigWig files and a conversion details file (`bam_to_bigwig_details.tsv`)

2. `bigwig_peakprep_summary.py`: Aggregates signal values from the generated bigWig files for specified regions
   - Uses multiBigwigSummary to extract signal values from each region
   - Output: Raw count table and log2-transformed table

**Usage example**:
```bash
python bigwig_peakprep.py \
    --bam_dir="/path/to/bam_files" \
    --annotation_file="/path/to/peaks.bed" \
    --output_dir="/path/to/output"
```

### 2. Differential Analysis Pipeline

**Main script**: `bigwig_peakprep_bamdiff.py`

This pipeline performs comparisons (primarily log2 ratios) between ChIP-seq samples and Input controls, calculating background-corrected signal values. Internally, it processes data through the following steps:

1. `bigwig_peakprep_bamdiff_generatediff.py`: Pairs samples with controls and performs differential comparison
   - Uses bamCompare to generate log2 ratio bigWig files for each sample/control pair
   - Output: Differential bigWig files and a conversion details file (`bamdiff_to_bigwig_details.tsv`)

2. `bigwig_peakprep_summary.py`: Aggregates signal values from the differential bigWig files
   - Note: Differential bigWig files are already log2 ratios, so no additional log2 transformation is performed

**Usage example**:
```bash
python bigwig_peakprep_bamdiff.py \
    --sample_bam_dir="/path/to/chipseq_bams" \
    --control_bam_dir="/path/to/input_bams" \
    --annotation_file="/path/to/peaks.bed" \
    --output_dir="/path/to/output"
```

### Utility Modules

- `bigwig_peakprep_utils.py`: Provides basic functions and common utilities
  - BAM to bigWig conversion (bamCoverage execution)
  - BED file processing using PyRanges
  - multiBigwigSummary execution and result processing
  - log2 transformation and various data format conversion functions

- `bigwig_peakprep_bamdiff_utils.py`: Provides specialized functions for differential analysis
  - Sample-control pairing
  - bamCompare execution functions
  - YAML conversion plan management
  - Result file generation functions

### Choosing Which Script to Use

- Use **Single Sample Analysis** (`bigwig_peakprep.py`) when:
  - You want to analyze each BAM file individually and obtain normalized signal values

- Use **ChIP-seq/Input Differential Analysis** (`bigwig_peakprep_bamdiff.py`) when:
  - You want to obtain background-corrected signal values (log2 ratios)
  - You have pairs of samples and controls to analyze together

## Output Directory Structure

The typical output directory structure will be as follows:

### For Single Sample Analysis (bigwig_peakprep.py)

```
output_dir/
├── bigwig/                              # bigWig file storage directory
│   ├── sample1.RPGC.bigwig              # Normalized bigWig files
│   ├── sample2.RPGC.bigwig
│   ├── ...
│   └── sampleN.RPGC.bigwig
├── summary/                             # multiBigwigSummary output directory
│   ├── bigwig_summary.npz               # Binary data (NumPy format)
│   ├── bigwig_summary.tab               # Tab-delimited text data
│   └── bigwig_summary.log               # Execution log
├── conversion_logs/                     # Conversion logs directory (if applicable)
│   └── peaks_converted.bed              # Converted BED file
├── bam_to_bigwig_details.tsv            # BAM→bigWig conversion details
├── multibigwig_counts.tsv               # Raw count table
├── multibigwig_counts_log2.tsv          # log2-transformed count table
└── pipeline_run.log                     # Pipeline execution log
```

### For Differential Analysis (bigwig_peakprep_bamdiff.py)

```
output_dir/
├── bigwig/                              # bigWig file storage directory
│   ├── sample1_vs_input1.log2ratio.bigwig  # log2 ratio bigWig files
│   ├── sample2_vs_input2.log2ratio.bigwig
│   ├── ...
│   └── sampleN_vs_inputN.log2ratio.bigwig
├── summary/                             # multiBigwigSummary output directory
│   ├── bigwig_summary.npz               # Binary data (NumPy format)
│   ├── bigwig_summary.tab               # Tab-delimited text data
│   └── bigwig_summary.log               # Execution log
├── conversion_logs/                     # Conversion logs directory (if applicable)
│   └── peaks_converted.bed              # Converted BED file
├── conversion_plan.yaml                 # Conversion plan YAML
├── conversion_plan.tsv                  # Conversion plan TSV format
├── bamdiff_to_bigwig_details.tsv        # Differential bigWig generation details
├── bamdiff_to_bigwig_details.yaml       # Details in YAML format
├── multibigwig_counts.tsv               # Differential value count table
└── pipeline_run.log                     # Pipeline execution log
```

## Argument List

Here is a comprehensive list of all arguments available for command-line execution, organized by the two main scripts:

### Single Sample Analysis (bigwig_peakprep.py) Arguments

| Argument | Description | Default Value |
|--------|------|-------------|
| **Required Arguments** |||
| `--bam_dir` | Directory containing BAM files (required for BAM→bigWig conversion) | - |
| `--annotation_file` | Path to the annotation file (BED, SAF, or merge_SE format) | - |
| `--output_dir` | Directory to save results | - |
| `--bigwig_details_file` | Path to the bam_to_bigwig_details.tsv file, required when running summary only (`--run_summary_only`) | None |
| **Tool Path Related** |||
| `--tools_dir` | Directory containing deeptools executables | None |
| `--bamcoverage_path` | Path to bamCoverage executable (ignored if `--tools_dir` is specified) | "bamCoverage" |
| `--multibigwigsummary_path` | Path to multiBigwigSummary executable (ignored if `--tools_dir` is specified) | "multiBigwigSummary" |
| **File/Sample Selection** |||
| `--filename_pattern` | Wildcard pattern for BAM filenames | "*.bam" |
| `--sample_delimiter` | Delimiter string to extract sample names from BAM filenames | None |
| **bamCoverage Related** |||
| `--bigwig_dir` | Directory to save bigWig files | output_dir/bigwig |
| `--bin_size` | bin size in bases for bamCoverage | 50 |
| `--normalize_using` | Normalization method for bamCoverage (RPGC, CPM, etc.) | "RPGC" |
| `--effective_genome_size` | Effective genome size for RPGC normalization | 2913022398 |
| **multiBigwigSummary Related** |||
| `--summary_dir` | Directory to save multiBigwigSummary files | output_dir/summary |
| `--summary_basename` | Base name for multiBigwigSummary output files | "bigwig_summary" |
| **Data Table Related** |||
| `--force_chr_start_end_ids` | Use chr_start_end format as PeakID even if BED file has name column | False |
| `--pseudocount` | Pseudocount for log2 transformation | 1.0 |
| `--raw_counts_name` | Filename for the raw counts table | "multibigwig_counts.tsv" |
| `--log2_counts_name` | Filename for the log2-transformed counts table | "multibigwig_counts_log2.tsv" |
| **Execution Control** |||
| `--threads` / `-T` | Number of threads for deeptools commands | 4 |
| `--force_overwrite` | Force overwrite existing bigWig files (skips them otherwise) | False |
| `--run_bamconvert_only` | Run only BAM→bigWig conversion | False |
| `--run_summary_only` | Run only bigWig summary processing | False |

### Differential Analysis (bigwig_peakprep_bamdiff.py) Arguments

| Argument | Description | Default Value |
|--------|------|-------------|
| **Required Arguments** |||
| `--sample_bam_dir` | Directory containing ChIP-seq sample BAM files | - |
| `--control_bam_dir` | Directory containing Input control (Input/IgG) BAM files | - |
| `--annotation_file` | Path to the annotation file (BED, SAF, or merge_SE format) | - |
| `--output_dir` | Directory to save results | - |
| **Tool Path Related** |||
| `--tools_dir` | Directory containing deeptools executables | None |
| `--bamcompare_path` | Path to bamCompare executable | "bamCompare" |
| **File Selection** |||
| `--sample_pattern` | Glob pattern for sample BAM files | "*.bam" |
| `--control_pattern` | Glob pattern for control BAM files | "\*Input\*.bam" |
| `--sample_delimiter` | Delimiter for sample name extraction | "_" |
| **Differential Analysis Related** |||
| `--operation` | bamCompare operation type (log2, ratio, subtract, add, mean, etc.) | "log2" |
| `--bin_size` | bamCompare bin size (in bases) | 50 |
| `--pseudocount` | Pseudocount value to avoid division by zero | 1.0 |
| `--effective_genome_size` | Effective genome size | 2913022398 |
| **Execution Control** |||
| `--threads` / `-T` | Number of threads | 4 |
| `--force_overwrite` | Force overwrite existing bigWig files | False |
| `--run_diff_only` | Run only differential bigWig generation (skip summary processing) | False |

## Detailed Usage Examples for Single Sample Analysis Pipeline (bigwig_peakprep.py)

### Basic Usage

```bash
python bigwig_peakprep.py \
    --bam_dir="/path/to/bam_files" \
    --annotation_file="/path/to/peaks.bed" \
    --output_dir="/path/to/output" \
    --tools_dir="/path/to/deeptools/bin" \
    --sample_delimiter="-H3K27ac.sorted.bam" \
    --normalize_using="RPGC" \
    --bin_size=50 \
    --threads=10
```

### Advanced Data Table Settings

```bash
python bigwig_peakprep.py \
    --bam_dir="/path/to/bam_files" \
    --annotation_file="/path/to/peaks.bed" \
    --output_dir="/path/to/output" \
    --tools_dir="/path/to/deeptools/bin" \
    --force_chr_start_end_ids \
    --pseudocount=0.1 \
    --raw_counts_name="my_raw_counts.tsv" \
    --log2_counts_name="my_log2_counts.tsv" \
    --threads=10
```

### Two-Stage Processing: BAM Conversion Only

```bash
python bigwig_peakprep.py \
    --bam_dir="/path/to/bam_files" \
    --annotation_file="/path/to/peaks.bed" \
    --output_dir="/path/to/output" \
    --tools_dir="/path/to/deeptools/bin" \
    --normalize_using="CPM" \
    --threads=10 \
    --run_bamconvert_only
```

### Two-Stage Processing: Summary Processing Only

```bash
python bigwig_peakprep.py \
    --annotation_file="/path/to/peaks.bed" \
    --output_dir="/path/to/output" \
    --tools_dir="/path/to/deeptools/bin" \
    --bigwig_details_file="/path/to/output/bam_to_bigwig_details.tsv" \
    --threads=10 \
    --run_summary_only
```

### Direct Execution Using Individual Step Scripts

Individual scripts are useful when finer control is needed or when you want to run only specific processing steps. Each script has its own options. You can check the detailed options by running each script with the `--help` flag.

BAM→bigWig conversion (Step 1):
```bash
python bigwig_peakprep_bamconvert.py \
    --bam_dir="/path/to/bam_files" \
    --output_dir="/path/to/output" \
    --tools_dir="/path/to/deeptools/bin" \
    --sample_delimiter="-H3K27ac.sorted.bam" \
    --threads=10
```

bigWig summary processing (Step 2):
```bash
python bigwig_peakprep_summary.py \
    --bigwig_details_file="/path/to/output/bam_to_bigwig_details.tsv" \
    --annotation_file="/path/to/peaks.bed" \
    --output_dir="/path/to/output" \
    --tools_dir="/path/to/deeptools/bin" \
    --threads=10
```

## Detailed Usage Examples for Differential Analysis Pipeline (bigwig_peakprep_bamdiff.py)

The differential analysis pipeline compares ChIP-seq samples to Input controls, correcting for background signal to calculate binding site intensities. It primarily performs log2 ratio calculations, allowing for more accurate binding site analysis.

### Basic Usage

```bash
python bigwig_peakprep_bamdiff.py \
    --sample_bam_dir="/path/to/chipseq_bams" \
    --control_bam_dir="/path/to/input_bams" \
    --annotation_file="/path/to/peaks.bed" \
    --output_dir="/path/to/output" \
    --tools_dir="/path/to/deeptools/bin" \
    --threads=10
```

### Automatic Sample and Control Matching

`bigwig_peakprep_bamdiff.py` automatically matches samples (ChIP-seq) with controls (Input) based on filename patterns. For example:

- Samples: `T1-H3K27ac.bam`, `T2-H3K27ac.bam`
- Controls: `T1-Input.bam`, `T2-Input.bam`

In this case, T1 sample will be paired with T1 control, and T2 sample with T2 control automatically. It's also possible to apply a single control to all samples.

### Specifying Sample and Control Using Specific Patterns

```bash
python bigwig_peakprep_bamdiff.py \
    --sample_bam_dir="/path/to/bams" \
    --control_bam_dir="/path/to/bams" \
    --sample_pattern="*-H3K27ac.bam" \
    --control_pattern="*-Input.bam" \
    --annotation_file="/path/to/peaks.bed" \
    --output_dir="/path/to/output" \
    --threads=10
```

### Differential bigWig Generation Only (Skip Summary Processing)

```bash
python bigwig_peakprep_bamdiff.py \
    --sample_bam_dir="/path/to/chipseq_bams" \
    --control_bam_dir="/path/to/input_bams" \
    --annotation_file="/path/to/peaks.bed" \
    --output_dir="/path/to/output" \
    --run_diff_only \
    --threads=10
```

### Changing Operation Type (Other than log2)

```bash
python bigwig_peakprep_bamdiff.py \
    --sample_bam_dir="/path/to/chipseq_bams" \
    --control_bam_dir="/path/to/input_bams" \
    --annotation_file="/path/to/peaks.bed" \
    --output_dir="/path/to/output" \
    --operation="subtract" \
    --threads=10
```

## Output File Details

### Single Sample Analysis Pipeline (bigwig_peakprep.py) Output

1. **bigWig files** (in the `bigwig/` directory):
   - Normalized bigWig files generated for each BAM file
   - Filename format: `{sample_name}.{normalization_method}.bigwig`
   - Example: `sample1.RPGC.bigwig`

2. **bam_to_bigwig_details.tsv**:
   - Detailed output file from Step 1 (input for Step 2)
   - Comment lines contain parameter information and execution metadata
   - Columns: `Sample_name`, `BAM_filename`, `BigWig_filename`, `BAM_fullpath`, `BigWig_fullpath`

3. **multiBigwigSummary output** (in the `summary/` directory):
   - `{summary_basename}.npz`: Binary data in NumPy format
   - `{summary_basename}.tab`: Tab-delimited text data

4. **Count Tables**:
   - `multibigwig_counts.tsv`: Raw count table (with PeakID information)
   - `multibigwig_counts_log2.tsv`: log2-transformed count table (with PeakID information)

### Differential Analysis Pipeline (bigwig_peakprep_bamdiff.py) Output

1. **Differential bigWig files** (in the `bigwig/` directory):
   - bigWig files representing log2 ratios between sample and control
   - Filename format: `{sample_name}_vs_{control_name}.log2ratio.bigwig`

2. **bamdiff_to_bigwig_details.tsv**:
   - TSV file containing details of the differential analysis
   - Columns: `Sample_name`, `BigWig_filename`, `BigWig_fullpath`
   - Important: This file includes `already_log2_transformed: true` metadata to prevent additional log2 transformation during summary processing

3. **Count Tables** (when summary processing is run):
   - `multibigwig_counts.tsv`: Table containing log2 ratio values for each peak region
   - Note: Values are already log2-transformed, so no additional transformation is performed (files specified by log2_counts_name option may not be generated)

Example count table format:
```
PeakID  Chr     Start   End     Sample1  Sample2  Sample3
peak1   chr1    999     2000    2.45     3.12     1.87
peak2   chr1    2999    4000    4.21     3.98     4.56
...
```

> **Note**: The `Start` column in output files is 0-based (BED format standard).

## Advanced Feature: YAML Conversion Plans

The differential analysis pipeline manages sample and control conversion plans in YAML format. This feature allows complex analysis settings to be easily saved and reused.

### YAML Conversion Plan Structure

The conversion plan file (`conversion_plan.yaml`) contains the following information:

```yaml
meta_info:
  generated_date: "2023-01-01 00:00:00"
  pipeline: "bigwig_peakprep_bamdiff pipeline"
parameters:
  operation: "log2"
  bin_size: 50
  pseudocount: 1.0
  effective_genome_size: 2913022398
  threads: 4
samples:
  T1-H3K27ac:
    bam_file: "T1-H3K27ac.bam"
    bam_path: "/path/to/T1-H3K27ac.bam"
controls:
  T1-Input:
    bam_file: "T1-Input.bam"
    bam_path: "/path/to/T1-Input.bam"
conversion_plan:
  - pair_id: "T1-H3K27ac_vs_T1-Input"
    sample_name: "T1-H3K27ac"
    control_name: "T1-Input"
    sample_bam: "T1-H3K27ac.bam"
    control_bam: "T1-Input.bam"
    output_bigwig: "T1-H3K27ac_vs_T1-Input.log2ratio.bigwig"
    sample_bam_path: "/path/to/T1-H3K27ac.bam"
    control_bam_path: "/path/to/T1-Input.bam"
    output_bigwig_path: "/path/to/output/bigwig/T1-H3K27ac_vs_T1-Input.log2ratio.bigwig"
```

### YAML Processing Results

Processing results are also saved in YAML format (`bamdiff_to_bigwig_details.yaml`), containing detailed information about the success/failure of each job:

```yaml
results:
  summary:
    total_jobs: 5
    completed_jobs: 4
    skipped_jobs: 1
    failed_jobs: 0
  job_results:
    - pair_id: "T1-H3K27ac_vs_T1-Input"
      status: "completed"
      start_time: "2023-01-01T00:00:00"
      end_time: "2023-01-01T00:05:00"
      duration_seconds: 300
      output_file: "/path/to/output/bigwig/T1-H3K27ac_vs_T1-Input.log2ratio.bigwig"
  successful_files: ["/path/to/output/bigwig/T1-H3K27ac_vs_T1-Input.log2ratio.bigwig", ...]
  failed_pairs: []
```

These conversion plan and result YAML files help with reproducing analysis settings and reprocessing specific samples in subsequent analyses.

## Normalization Method Selection

### Single Sample Analysis (bigwig_peakprep.py) Normalization

You can select from the following normalization methods using the `--normalize_using` option:

- **RPGC** (default): Reads Per Genomic Content. Normalized based on the specified genome size
- **CPM**: Counts Per Million mapped reads
- **BPM**: Bins Per Million mapped reads
- **RPKM**: Reads Per Kilobase per Million mapped reads
- **None**: No normalization

For example, to use CPM normalization:
```bash
python bigwig_peakprep.py \
    --bam_dir="/path/to/bam_files" \
    --annotation_file="/path/to/peaks.bed" \
    --output_dir="/path/to/output" \
    --normalize_using="CPM" \
    --threads=10
```

### Differential Analysis (bigwig_peakprep_bamdiff.py) Operation Types

You can select from the following operation types using the `--operation` option:

- **log2** (default): log2(sample/control) ratio (pseudocount specified by `--pseudocount`)
- **ratio**: sample/control ratio
- **subtract**: sample - control subtraction
- **add**: sample + control addition
- **mean**: (sample + control)/2 average
- **reciprocal_ratio**: control/sample inverse ratio
- **first**: output sample only
- **second**: output control only

For example, to perform a simple subtraction:
```bash
python bigwig_peakprep_bamdiff.py \
    --sample_bam_dir="/path/to/chipseq_bams" \
    --control_bam_dir="/path/to/input_bams" \
    --annotation_file="/path/to/peaks.bed" \
    --output_dir="/path/to/output" \
    --operation="subtract" \
    --threads=10
```

## Troubleshooting

Common issues and solutions:

1. **Errors with argument values containing hyphens**:
   - Use equals sign format for arguments: `--sample_delimiter="-H3K27ac.sorted.bam"`
   - Do not use space-delimited format: ~~`--sample_delimiter "-H3K27ac.sorted.bam"`~~

2. **bamCoverage/multiBigwigSummary not found errors**:
   - Specify the correct deeptools bin directory with the `--tools_dir` option
   - Alternatively, specify full paths with `--bamcoverage_path` and `--multibigwigsummary_path`

3. **Sample-control automatic matching problems**:
   - Ensure sample and control filenames have common identifiers (T1, T2, etc.)
   - Use `--sample_pattern` and `--control_pattern` to explicitly specify patterns if needed

4. **Out of memory errors**:
   - When processing a large number of files or large genomic regions, run on a machine with more memory
   - Consider reducing the number of files processed or process them in batches

5. **Slow processing**:
   - Increase the number of threads with the `--threads` option
   - Use indexed BAM files
   - Increasing the `--bin_size` value will speed up processing but decrease resolution

### Differential Analysis Specific Issues

6. **Sample-control matching failure**:
   - If you see warning message "Could not create sample-control pairs" in the log, matching has failed
   - Use `--sample_pattern` and `--control_pattern` to explicitly specify filename patterns
   - Ensure filenames have common identifiers (T1, T2, etc.)
   - If using a single control for all samples, ensure the control directory contains only one file

7. **Errors after differential bigWig generation**:
   - If differential bigWig files exist but summary processing fails, binary format might be corrupted
   - Use `--force_overwrite` option to regenerate differential bigWig files
   - Check the log (`bamCompare.log`) for files that failed to generate

8. **Verify correct script is being used**:
   - For sample and control differential analysis, ensure you're using `bigwig_peakprep_bamdiff.py` not `bigwig_peakprep.py`
   - Conversely, for single sample analysis, use `bigwig_peakprep.py`

## Citation

If you use this tool in your research, please cite:
**(Paper in preparation)**

## License

This program is released under the MIT license. For details, see [LICENSE](https://github.com/hamamoto-lab/SEgene/blob/main/LICENSE)