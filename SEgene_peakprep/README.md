# SEgene_peakprep

*(For the Japanese version of this README, please see [README_ja.md](./README_ja.md).)*

**SEgene_peakprep** is a Python pipeline for quantifying and normalizing signal values in specified genomic regions from sequence data (BAM files), and formatting them for downstream SEgene analysis. This tool is not dependent on any specific reference genome and can be used with data mapped to genomes such as hg38.

## Development Status

> **Note:** This tool is **currently under development**, and the program file structure and usage may change between versions.

## Overview

SEgene_peakprep serves as the initial data preparation step in the [SEgene project](https://github.com/hamamoto-lab/SEgene):

1.  **Input**: Mapped sequence data (BAM files) and annotation files (BED files, etc.)
2.  **Output**: Normalized count value tables
3.  **Flow**: This output serves as the primary input for the `peak_to_gene_links` module, which is combined with RNA-seq TPM data for gene regulatory network analysis

This pipeline offers two implementation methods:

1.  **CPM Method**: Uses `featureCounts` to directly calculate Counts Per Million (CPM) values
2.  **BigWig Method**: Uses `deeptools` to convert BAM to bigWig files, then quantifies peaks using `multiBigwigSummary`

Both methods support BED, SAF, and `merge_SE.tsv` format annotation files. You can choose the appropriate method based on your needs.

## Common Requirements

- **Python**: Python 3.10 or higher
- **Basic Libraries**:
  - pandas
  - numpy
  - natsort
  
  ```bash
  # Using conda
  conda install pandas numpy natsort
  
  # Using pip
  pip install pandas numpy natsort
  ```

### Common Input Files

Both implementation methods require the following input files:

1. **BAM Files**: Mapped sequence data files (recommended: indexed .bai files)
2. **Annotation Files**: Files that define genomic regions, supporting these formats:
   - **BED Format**: Standard BED format (0-based start coordinates)
     ```
     # Example: peaks.bed (tab-delimited)
     chr1    999     2000    peak1    0    +
     chr1    2999    4000    peak2    0    +
     chr2    4999    6000    peak3    0    -
     ```
   
   - **SAF Format** (Simplified Annotation Format): Tab-delimited text file (**1-based start coordinates**)
     ```
     # Example: peaks.saf (tab-delimited)
     GeneID  Chr     Start   End     Strand
     peak1   chr1    1000    2000    +
     peak2   chr1    3000    4000    +
     peak3   chr2    5000    6000    -
     ```
     **Note:** SAF format coordinates are 1-based (Start=1000 is 999 in 0-based), but all outputs are unified to 0-based
   
   - **merge_SE.tsv Format**: First column is region ID in the format `chr_start_end` (e.g., `chr1_1000_2000`)
     ```
     # Example: merged_se.tsv (tab-delimited)
     se_data            feature
     chr1_999_2000      feature1
     chr1_2999_4000     feature2
     chr2_4999_6000     feature3
     ```
     **Note:** In merge_SE format, coordinates are interpreted as 0-based (chr1_999_2000's Start is 999 in 0-based)

## CPM Method (cpm_peakprep)

The CPM method uses `featureCounts` to directly count reads in specified genomic regions (BED/SAF/merge_SE.tsv format) and calculates CPM values and their log-transformed values.

### Additional Requirements

- **External Tools**:
  - samtools (recommended: v1.13 or higher, tested with v1.15.1)
  - featureCounts (part of the Subread package, recommended: v2.0.0 or higher, tested with v2.0.6)
  
  ```bash
  # Installation example with conda
  conda install -c bioconda samtools subread
  ```

### Detailed Guide

For detailed configuration options and usage of the CPM method, please refer to the following:

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

The BigWig method uses `deeptools` `bamCoverage` to convert BAM files to normalized bigWig files, then extracts count values from specified regions using `multiBigwigSummary`.

### Additional Requirements

- **Additional Libraries**:
  - pyranges
  
  ```bash
  # Using conda
  conda install -c bioconda pyranges
  
  # Using pip
  pip install pyranges
  ```

- **External Tools**:
  - `deeptools` (recommended: v3.5.0 or higher, tested with v3.5.5)
  
  ```bash
  # Installation example with conda
  conda install -c bioconda deeptools
  ```

### Pipeline Structure

The BigWig method consists of three scripts:

1. **bigwig_peakprep.py**: Main wrapper script (controls the entire process)
2. **bigwig_peakprep_bamconvert.py**: Handles BAM to bigWig file conversion
3. **bigwig_peakprep_summary.py**: Generates multi-sample count tables from bigWig files using `multiBigwigSummary`

This three-stage processing flow generates normalized bigWig files from BAM files and aggregates and quantifies signal values in specified genomic regions. Each step can be executed individually or as a series of processes through the wrapper script.

By default, `bamCoverage` uses **RPGC (Reads Per Genomic Content)** normalization, but you can change to other normalization methods (CPM, BPM, RPKM, RPGC, None) with the `--normalize_using` option.

### Basic Usage

```bash
python bigwig_peakprep.py \
    --bam_dir="/path/to/bam_files" \
    --annotation_file="/path/to/peaks.bed" \
    --output_dir="/path/to/output" \
    --tools_dir="/path/to/deeptools/bin" \
    --sample_delimiter="-H3K27ac.sorted.bam" \
    --threads=10
```

### Key Arguments

| Argument | Description | Default Value |
|--------|------|-------------|
| **Required Arguments** |||
| `--bam_dir` | Directory containing BAM files | - |
| `--annotation_file` | Annotation file (BED, SAF, `merge_SE.tsv` format) | - |
| `--output_dir` | Directory to save results | - |
| **Tool Paths** |||
| `--tools_dir` | Directory containing `deeptools` executables | None |
| `--bamcoverage_path` | Path to `bamCoverage` executable | "bamCoverage" |
| `--multibigwigsummary_path` | Path to `multiBigwigSummary` executable | "multiBigwigSummary" |
| **Execution Step Control** |||
| `--run_bamconvert_only` | Run only BAM to bigWig conversion | False |
| `--run_summary_only` | Run only bigWig summary processing | False |
| **Normalization Related** |||
| `--normalize_using` | `bamCoverage` normalization method (RPGC, CPM, BPM, RPKM, None) | "RPGC" |
| `--effective_genome_size` | Effective genome size for RPGC normalization | 2913022398 |

To see all options, run:
```bash
python bigwig_peakprep.py --help
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

1. **bigWig Files** (Step 1 output):
   - Normalized bigWig files generated from each BAM file
   - Filename format: `{sample_name}.{normalization_method}.bigwig`

2. **Count Tables** (Step 2 output):
   - `multibigwig_counts.tsv`: Raw count table
   - `multibigwig_counts_log2.tsv`: log2-transformed count table

Example output table format:
```
PeakID  Chr     Start   End     Sample1  Sample2  Sample3
peak1   chr1    999     2000    2.45     3.12     1.87
peak2   chr1    2999    4000    4.21     3.98     4.56
...
```

## Troubleshooting

Common issues and solutions:

1. **Tool not found error**:
   - Verify that samtools, featureCounts, or deeptools are installed
   - Specify the full path using the corresponding option

2. **BAM files not found**:
   - Check the `--filename_pattern` argument

3. **Sample name extraction issues**:
   - Set the `--sample_delimiter` value appropriately and use an equals sign (`=`)
   - For values with hyphens, always use the equals sign: `--sample_delimiter="-H3K27ac.sorted.bam"`

4. **Out of memory error**:
   - Run on a machine with more memory
   - Reduce the number of files processed or process them in batches

5. **Slow processing**:
   - Increase the number of threads with the `--threads` option
   - Use indexed BAM files

## Citation

If you use this tool in your research, please cite:
**(Paper in preparation)**

## License

This program is released under the MIT license. For details, see [LICENSE](./LICENSE).