# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.0] - 2024-02-07

### Added
- **New CLI tools (`cli_tools/`)** for streamlined featureCounts processing, CPM calculation, and SE–gene correlation analysis
  - Scripts such as `generate_gtf.py`, `generate_file_list.sh`, `run_featurecounts_array.sh`, etc.
  - `correlation_analysis.py` computes Pearson correlation between SE (CPM) and genes (TPM)
- **Documentation updates**:
  - Added **English (`README.md`)** and **Japanese (`README_ja.md`)** documentation under `/cli_tools/`
  - Updated the main repository documentation (`README.md`) and `/SE_to_gene_links/` docs to explain how to integrate the new CLI tools

### Changed
- Refined repository structure by separating stand-alone CLI scripts into `cli_tools/`, independent of `SEgene_package`

### Notes
- The newly added documents provide usage instructions and workflow details in both English and Japanese  
- Enables advanced or extended analyses on top of the existing `SE_to_gene_links` functionalities

## [1.0.0] - 2024-01-15

### Added
- Initial release of SEgene, a comprehensive platform for identifying and analyzing Super-Enhancer-to-gene links through statistical approaches
  - `peak_to_gene_links`:
    - Retrieves correlation information between gene expression and enhancer peaks
    - Data analysis and visualization support
  - `SE_to_gene_links`:
    - Evaluates and analyzes super-enhancer to gene links using correlation data
    - Graph theory-based data visualization
    - Interactive analysis with Jupyter Notebook
- Detailed installation instructions and documentation in the GitHub repository

[1.1.0]: https://github.com/hamamoto-lab/SEgene/releases/tag/v1.1.0
[1.0.0]: https://github.com/hamamoto-lab/SEgene/releases/tag/v1.0.0
