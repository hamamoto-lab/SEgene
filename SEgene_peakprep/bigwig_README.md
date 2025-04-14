# SEgene_peakprep - BigWig Method Detailed Guide

This document provides detailed usage instructions and examples specific to the **BigWig method** of SEgene_peakprep. For basic overview and common requirements, please refer to the [main README](./README.md).

## BigWig Method Overview

The BigWig method uses `deeptools`'s `bamCoverage` to convert BAM files to normalized bigWig files, then uses `multiBigwigSummary` to extract count values from specified regions. This approach offers the following advantages:

- Visualization of coverage data across the entire genome
- Support for various normalization methods (RPGC, CPM, BPM, RPKM)
- Efficient processing and comparison of multi-sample data

## Pipeline Structure

The BigWig method consists of three scripts:

1. **bigwig_peakprep.py**: Main wrapper script (controls the entire process)
2. **bigwig_peakprep_bamconvert.py**: Handles BAM to bigWig file conversion
3. **bigwig_peakprep_summary.py**: Generates multi-sample count tables from bigWig files

This three-stage processing flow generates normalized bigWig files from BAM files and aggregates and quantifies signal values in specified genomic regions. Each step can be executed individually or as a series of processes through the wrapper script.

## Complete Argument List

Here is a comprehensive list of all arguments available for command-line execution:

| Argument | Description | Default Value |
|--------|------|-------------|
| **Required Arguments** |||
| `--bam_dir` | Directory containing BAM files (required for BAM→bigWig conversion) | - |
| `--annotation_file` | Path to the annotation file (BED, SAF, or merge_SE format) | - |
| `--output_dir` | Directory to save results | - |
| `--bigwig_details_file` | Path to bam_to_bigwig_details.tsv file, required when running summary only (`--run_summary_only`) | None |
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
| **Logging Settings** |||
| `--log_level` | Logging level (DEBUG/INFO/WARNING/ERROR/CRITICAL) | "INFO" |
| `--script_log_file` | Filename for script execution log (empty string to disable)        | "pipeline_run.log" |

## Detailed Usage Examples

### Wrapper Script for Complete Processing

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

## Output File Details

### Step 1 (BAM to bigWig Conversion) Output

1. **bigWig files** (in the `bigwig/` directory):
   - Normalized bigWig files generated for each BAM file
   - Filename format: `{sample_name}.{normalization_method}.bigwig`
   - Example: `sample1.RPGC.bigwig`

2. **bam_to_bigwig_details.tsv**:
   - Detailed output file from Step 1 (input for Step 2)
   - Comment lines contain parameter information and execution metadata
   - Columns: `Sample_name`, `BAM_filename`, `BigWig_filename`, `BAM_fullpath`, `BigWig_fullpath`

### Step 2 (bigWig Summary and Count Table Generation) Output

1. **multiBigwigSummary output** (in the `summary/` directory):
   - `{summary_basename}.npz`: Binary data in NumPy format
   - `{summary_basename}.tab`: Tab-delimited text data

2. **Count Tables**:
   - `multibigwig_counts.tsv`: Raw count table (with PeakID information)
   - `multibigwig_counts_log2.tsv`: log2-transformed count table (with PeakID information)

Example count table format:
```
PeakID  Chr     Start   End     Sample1  Sample2  Sample3
peak1   chr1    999     2000    2.45     3.12     1.87
peak2   chr1    2999    4000    4.21     3.98     4.56
...
```

> **Note**: The `Start` column in output files is 0-based (BED format standard).

## Normalization Method Selection

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

## Troubleshooting

Common issues and solutions:

1. **Errors with argument values containing hyphens**:
   - Use equals sign format for arguments: `--sample_delimiter="-H3K27ac.sorted.bam"`
   - Do not use space-delimited format: ~~`--sample_delimiter "-H3K27ac.sorted.bam"`~~

2. **bamCoverage/multiBigwigSummary not found errors**:
   - Specify the correct deeptools bin directory with the `--tools_dir` option
   - Alternatively, specify full paths with `--bamcoverage_path` and `--multibigwigsummary_path`

3. **Sample name extraction issues**:
   - Set the `--sample_delimiter` value appropriately and use an equals sign (`=`)
   - To check BAM filename to sample name mapping, set the log level to `DEBUG`

4. **Out of memory errors**:
   - When processing a large number of BAM files or large genomic regions, run on a machine with more memory
   - Consider reducing the number of files processed or process them in batches

5. **Slow processing**:
   - Increase the number of threads with the `--threads` option
   - Use indexed BAM files
   - Increasing the `--bin_size` value will speed up processing but decrease resolution

## Citation

If you use this tool in your research, please cite:
**(Paper in preparation)**

## License

This program is released under the MIT license. For details, see [LICENSE](./LICENSE).