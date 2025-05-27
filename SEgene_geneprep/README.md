# SEgene_geneprep

*(For the Japanese version of this README, please see [README_ja.md](./README_ja.md).)*

**SEgene_geneprep** is a Python tool designed to extract gene position information from GTF files and merge it with expression data from TPM files. Its primary purpose is to convert salmon output from the nf-core/rnaseq pipeline into a CSV file format suitable for RNAinput in the SEgene analysis pipeline.

## Overview

SEgene_geneprep serves as a utility tool in the [SEgene project](https://github.com/hamamoto-lab/SEgene):

1. **Input**: 
   - GTF files (gene annotations):
     - Can be obtained using the `--save_reference` option in nf-core/rnaseq
     - For external sources, genome annotations from [Illumina iGenomes](https://support.illumina.com/sequencing/sequencing_software/igenome.html) (e.g., GRCh38) can also be used
     - In either case, the file must conform to standard GTF format
   - TPM files (gene expression data): nf-core/rnaseq output files such as `salmon.merged.gene_tpm.tsv`

2. **Output**: Integrated tables containing gene IDs, genomic position information (chromosome, start, end, strand), and expression values

3. **Flow**: This output is used as RNAseq input for SEgene's `peak_to_gene_links` module, leading to subsequent `SE_to_gene_links` analysis

For more details about nf-core/rnaseq, please refer to the official repository (https://github.com/nf-core/rnaseq).

## Requirements and Dependencies

### Python Environment
- Python 3.6 or higher

### Python Libraries
- pandas (version 1.5 or higher recommended)

## Usage

### Basic Usage

```bash
python geneprep.py --gtf path/to/genes.gtf --tpm path/to/salmon.merged.gene_tpm.tsv
```

### All Options

```
usage: geneprep.py [-h] --gtf GTF --tpm TPM [--output OUTPUT] [--output-dir OUTPUT_DIR] 
                   [--unmatched UNMATCHED] [--mismatch MISMATCH] 
                   [--undefined-strand UNDEFINED_STRAND] [--join {inner,left,outer}] 
                   [--log LOG] [--no-create-dirs] [--verbose]

Extract gene position information from GTF and merge with TPM data

optional arguments:
  -h, --help            show this help message and exit
  --gtf GTF             Path to GTF file
  --tpm TPM             Path to TPM file
  --output OUTPUT       Output file name (default: gene_tpm_with_position.csv)
  --output-dir OUTPUT_DIR
                        Output directory for all files (overrides individual file paths)
  --unmatched UNMATCHED Output file for unmatched GTF lines (default: unmatched_gtf_lines.tsv)
  --mismatch MISMATCH   Output file for gene_id/gene_name mismatches (default: id_name_mismatches.tsv)
  --undefined-strand UNDEFINED_STRAND
                        Output file for lines with undefined strand (default: undefined_strand_lines.tsv)
  --join {inner,left,outer}
                        Join method:
                        inner (only genes in both files),
                        left (all genes in GTF),
                        outer (all genes in either file)
                        (default: inner)
  --log LOG             Path to log file (if specified, logs will be saved to this file)
  --no-create-dirs      Do not create output directories if they don't exist
                        (by default, directories are created automatically)
  --verbose, -v         Enable verbose output mode
```

### Join Method (--join) Description

- **inner**: Default setting. Includes only genes present in both files (GTF and TPM).
- **left**: Includes all genes in the GTF file. If a gene has no corresponding data in the TPM file, expression values will be output as missing values (NaN).
- **outer**: Includes all genes present in either file. For genes present in only one file, data from the other file will be treated as missing values.

## Input File Formats

### GTF File

Supports standard GTF (Gene Transfer Format) files used in the nf-core/rnaseq pipeline. Here's an example:

```
chr1    BestRefSeq      exon    11874   12227   .       +       .       gene_id "DDX11L1"; gene_name "DDX11L1"; transcript_id "rna0"; tss_id "TSS31672";
chr1    BestRefSeq      exon    12613   12721   .       +       .       gene_id "DDX11L1"; gene_name "DDX11L1"; transcript_id "rna0"; tss_id "TSS31672";
chr1    BestRefSeq      exon    13221   14409   .       +       .       gene_id "DDX11L1"; gene_name "DDX11L1"; transcript_id "rna0"; tss_id "TSS31672";
```

Position information is integrated at the gene level from each exon.

**How to obtain GTF files**:
- **Via nf-core/rnaseq**: Using the `--save_reference` option during pipeline execution will save the reference GTF file used
- **Via Illumina iGenomes**: Standard genome annotations (e.g., GRCh38, GRCm39) can be obtained from [Illumina iGenomes](https://support.illumina.com/sequencing/sequencing_software/igenome.html)
- **Note**: Regardless of the source, ensure the GTF file **contains exon features**

**Important specifications**:
- Processing target: **exon features only** (CDS, UTR, etc. are excluded)
- Lines with strand information (column 7) as '.' (undefined) are skipped. These lines are recorded separately in the `undefined_strand_lines.tsv` file
- gene_id and gene_name attributes are required

### TPM File

Tab-delimited TPM (Transcripts Per Million) data file output by the nf-core/rnaseq pipeline. Primarily targets `salmon.merged.gene_tpm.tsv`, but other salmon output files with similar format can also be converted. Here's an example:

```
gene_id gene_name       SRX2370497
A1BG    A1BG    0.067802
A1BG-AS1        A1BG-AS1        1.100025
A1CF    A1CF    0
```

The `gene_id` column is required and used as the key for joining with the GTF file.

## Output Files

### Main Output File (default: gene_tpm_with_position.csv)

```
symbol,chrom,strand,start,end,sample1,sample2,...
A1BG,chr19,+,58345178,58353492,0.067802,0.123456,...
A1BG-AS1,chr19,-,58353066,58358438,1.100025,2.345678,...
```

Column descriptions:
- **symbol**: Gene ID (gene_id from GTF file)
- **chrom**: Chromosome name
- **strand**: Transcription direction (+, -) ※ Lines with undefined strand ('.') are not included in the output
- **start**: Start position (0-based)
- **end**: End position
- **sample1, sample2, ...**: Expression levels for each sample from the TPM file

### Other Output Files

- **unmatched_gtf_lines.tsv**: Lines that failed to parse in the GTF file
- **id_name_mismatches.tsv**: List of genes where gene_id and gene_name do not match
- **undefined_strand_lines.tsv**: List of lines with undefined strand ('.')

## Example Executions

### Basic execution using nf-core/rnaseq output

```bash
python geneprep.py --gtf /path/to/nfcore/references/genes.gtf --tpm /path/to/nfcore/results/salmon.merged.gene_tpm.tsv
```

### Specifying a unified output directory

```bash
# Save all output files to the same directory
python geneprep.py --gtf genes.gtf --tpm salmon.merged.gene_tpm.tsv --output-dir /path/to/output
```

### Verbose mode with log file output

```bash
python geneprep.py --gtf genes.gtf --tpm salmon.merged.gene_tpm.tsv --verbose --log processing.log
```

### Execution with outer join (including genes in GTF but not in TPM)

```bash
python geneprep.py --gtf genes.gtf --tpm salmon.merged.gene_tpm.tsv --join outer --output all_genes.csv
```

### Changing the output destination for undefined strand lines

```bash
python geneprep.py --gtf genes.gtf --tpm salmon.merged.gene_tpm.tsv --undefined-strand undefined_strands.txt
```

## Detailed Statistics Display

Using the `--verbose` option will display detailed statistics like the following:

```
Statistics:
Number of genes in GTF: 25631
Number of genes in TPM file: 23486
Number of matched genes: 22978
Number of genes in GTF but not in TPM: 2653
Number of genes in TPM but not in GTF: 508

Detailed Statistics:
Chromosome distribution in GTF (top 10):
  chr1: 2103 genes
  chr2: 1283 genes
  chr3: 1078 genes
  ...

Gene length statistics:
  Mean: 32456.1 bp
  Median: 16782.0 bp
  Min: 102.0 bp
  Max: 2304562.0 bp
```

## Troubleshooting

Common issues and solutions:

1. **File not found error**:
   - Check that the GTF and TPM file paths are correct
   - Verify that file permissions are properly set

2. **Many lines with undefined strand**:
   - Check the `undefined_strand_lines.tsv` file to see which genes were skipped
   - Consider using a different GTF file if necessary

3. **Out of memory error**:
   - When processing large GTF files, run on a machine with more memory

## Notes

- Ensure sufficient memory is available when processing large GTF files
- The TPM file must contain a `gene_id` column
- **Only exon features in the GTF file are processed**. CDS, UTR, and other features are skipped
- **Lines with undefined strand ('.') are excluded from processing** and recorded in a separate file
- Gene position information is integrated from exon information, using the minimum start position and maximum end position
- **Output directories are automatically created by default if they don't exist**. This behavior can be disabled with the `--no-create-dirs` option

## License

This program is released under the MIT License. For more details, please refer to the [LICENSE](https://github.com/hamamoto-lab/SEgene/blob/main/LICENSE) file.

## Citation

If you use this tool in your research, please cite:

Shinkai, N., Asada, K., Machino, H., Takasawa, K., Takahashi, S., Kouno, N., Komatsu, M., Hamamoto, R., & Kaneko, S. (2025). SEgene identifies links between super enhancers and gene expression across cell types. *npj Systems Biology and Applications*, 11(1), 49. https://doi.org/10.1038/s41540-025-00533-x