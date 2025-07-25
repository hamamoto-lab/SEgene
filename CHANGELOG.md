# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.13.1] - 2025-07-25

### Added
- **SEgene_analyzer_erna metadata documentation**:
  - Added `SEgene_analyzer_erna/metadata/README.md` - Documentation for eRNAbase metadata files
  - Added `SEgene_analyzer_erna/metadata/README_ja.md` - Japanese documentation for eRNAbase metadata files
  - Both files include data source links to eRNAbase and version information (July 22, 2025)
  - Streamlined documentation by removing verbose technical descriptions and focusing on essential information

## [1.13.0] - 2025-07-23

### Added
- **New development version: SEgene_analyzer_erna**:
  - eRNAbase-specific analysis tool derived from SEgene_region_package
  - Command-line interface with 4 specialized commands: `report`, `single`, `batch`, `prepare-stats`
  - pip-installable package with `erna-analyzer` system-wide command
  - Support for multiple BED files and Parquet metadata format specific to eRNAbase
  - Species filtering functionality (human/mouse) for eRNAbase datasets
  - Comprehensive bilingual documentation (English and Japanese)
  - Enhanced chromosome filtering and validation for genomic region analysis
  - WSL-compatible matplotlib backend configuration for headless environments
  - Complete test suite with real eRNAbase data validation (858 samples)
  - Statistical enrichment analysis with Fisher's exact test and FDR correction
  - Multiple output formats: PNG, SVG, PDF, CSV, HTML reports
  - Positioned as eRNAbase-specialized counterpart to SEgene_analyzer

### Changed
- **Documentation updates**:
  - Added SEgene_analyzer_erna to "Development Versions" section in main README files
  - Updated usage links to include the new development version documentation

## [1.12.0] - 2025-07-18

### Added
- **New development version: SEgene_analyzer**:
  - Development version of SEgene_region_analyzer with enhanced CLI interface and modern Python packaging
  - Command-line interface with 4 specialized commands: `prepare-stats`, `batch`, `single`, `report`
  - pip-installable package with `sedb-analyzer` system-wide command
  - Advanced caching system for improved batch processing performance
  - Comprehensive bilingual documentation (English and Japanese)
  - Test suite with 29 test cases covering core functionality
  - Enhanced matplotlib backend configuration for WSL environments
  - Multiple output formats: PNG, SVG, PDF, EPS
  - Positioned as experimental development version alongside stable SEgene_region_analyzer

### Changed
- **Documentation structure updates**:
  - Added "Development Versions" section to main README files
  - Updated Program Structure section to mention development versions availability
  - Enhanced README files with bold emphasis on development program information
  - Improved navigation with direct links to development version documentation

## [1.11.0] - 2025-05-27

### Published
- **🎉 SEgene paper officially published**: 
  - Shinkai, N., Asada, K., Machino, H., Takasawa, K., Takahashi, S., Kouno, N., Komatsu, M., Hamamoto, R., & Kaneko, S. (2025). SEgene identifies links between super enhancers and gene expression across cell types. *npj Systems Biology and Applications*, 11(1), 49. https://doi.org/10.1038/s41540-025-00533-x
  - Test data and supplementary materials released on Figshare: https://doi.org/10.6084/m9.figshare.28171127

### Changed
- **Documentation updated following paper publication**:
  - Updated all README files to include the official paper citation
  - Replaced all "manuscript in preparation" notices with published paper information
  - Updated CITATION file with complete publication details and BibTeX format
  - Standardized citation section titles in Japanese READMEs to "引用 / Citation"
  - Updated Figshare links throughout documentation from placeholders to official DOI

### Added
- **Citation sections**: Added missing citation sections to ensure all components include proper citation information

## [1.10.0] - 2025-05-02

### Added
- **Enhanced SEgene_peakprep with edgeR normalized CPM**:
  - Added edgeR normalized CPM (calcnorm CPM) as the third quantification method in addition to Standard log2-CPM and BigWig methods
  - Implemented integration with R's edgeR package via rpy2 for advanced normalization of ChIP-seq count data
  - Added multiple normalization methods from edgeR: upperquartile (default), TMM, RLE
  - Implemented optional filtering capabilities based on CPM values across samples
  - Added comprehensive documentation in both English and Japanese
  - Created a dedicated edgeR normalized CPM guide (cpm_calcnorm_README.md)

### Changed
- **Updated SEgene_peakprep code and documentation**:
  - Reorganized CPM method to support both standard log2-CPM and edgeR normalized CPM calculations
  - Updated command-line interface with new `--calcnorm` parameter and related options for normalization control
  - Improved sample name handling for better compatibility with complex file naming conventions
  - Updated all relevant documentation (README.md, README_ja.md, cpm_README.md, cpm_README_ja.md) to explain the new capabilities
- **Dependencies**:
  - Added R (≥4.2.0), edgeR (≥3.40.0), and rpy2 (≥3.5.0) as new dependencies for the edgeR normalized CPM functionality
  - Updated installation instructions to explain the setup of R and Bioconductor environments

## [1.9.0] - 2025-04-24

### Added
- **New component: SEgene_geneprep**:
  - Preparation tool for RNA-seq input data to add region information for P2GL input formats
  - Extracts gene position information from GTF files and merges with expression data from TPM files
  - Converts salmon output from nf-core/rnaseq pipeline to CSV format suitable for RNAinput in SEgene analysis
  - Supports multiple input formats including standard GTF and salmon TPM files
  - Compatible with genome annotations from Illumina iGenomes (e.g., GRCh38)
  - Comprehensive documentation in both English and Japanese
- **Workflow visualization**:
  - Added mermaid diagrams to visualize the complete SEgene workflow
  - Clear illustration of the four-step process: P2GL Preparation, P2GL Analysis, SE Analysis, and Region Evaluation

### Changed
- **Program structure reorganization**:
  - Reorganized workflow into four primary components for clearer understanding
  - Positioned SEgene_peakprep and SEgene_geneprep as the initial P2GL data preparation steps
  - Updated documentation structure to reflect the new workflow organization
- **Documentation updates**:
  - Comprehensive revision of main README files (both English and Japanese)
  - Restructured program structure sections for better clarity
  - Enhanced descriptions of component relationships and data flow

## [1.8.0] - 2025-04-17

### Added
- **Enhanced SEgene_peakprep with differential analysis**:
  - Added differential analysis functionality for ChIP-seq versus Input control data within the bigWig processing workflow.
  - This feature utilizes `bamCompare` from the `deeptools` suite to identify significantly enriched regions by comparing ChIP signals against Input signals.

### Changed
- **Updated SEgene_peakprep code and documentation**:
  - Updated relevant scripts within the `SEgene_peakprep` directory (including `bigwig_peakprep_bamdiff.py`, `bigwig_peakprep_bamdiff_utils.py`, `bigwig_peakprep.py`, etc.) to integrate the new differential analysis capability.
  - Correspondingly revised associated documentation (`README.md`, `README_ja.md`, `bigwig_README.md`, `bigwig_README_ja.md`) to accurately reflect the updated features and their usage instructions.

## [1.7.0] - 2025-04-14

### Added
- **New BigWig method in SEgene_peakprep**:
  - Implemented a new alternative approach for peak quantification using BigWig files
  - Added `deeptools` integration for BAM to BigWig conversion and signal quantification
  - Added `bamCoverage` functionality with multiple normalization options (RPGC, CPM, BPM, RPKM)
  - Added `multiBigwigSummary` processing for extracting signal values from specified regions
  - Implemented automatic log2 transformation of signal values
  - Created comprehensive utilities for processing and handling BigWig files and genomic coordinates
- **Enhanced processing workflow**:
  - Added three-stage workflow with separate scripts for different processing steps
  - Added wrapper script for seamless execution of the entire pipeline
  - Added option to run individual processing steps separately

### Changed
- **Expanded SEgene_peakprep capabilities**:
  - Updated documentation to reflect both CPM and BigWig implementation methods
  - Expanded command-line interface to support both processing approaches
  - Standardized output formats across both methods for downstream compatibility

### Fixed
- **Corrected terminology in SEgene_peakprep**:
  - Fixed references to "merge_SV.tsv" format from v1.6.1, which should correctly be "merge_SE.tsv" (Super-Enhancer) format
  - Updated all related documentation and code comments to use consistent terminology
  - Fixed `--is_mergesv_format` flag name to `--is_mergese_format` for specifying merge_SE.tsv format annotations

## [1.6.1] - 2025-04-11

### Added
- **SEgene_peakprep improvements**:
  - Added direct support for BED and merge_SE.tsv input formats in annotation files
  - Automatic conversion of BED and merge_SE formats to SAF format internally
  - Added `--is_mergese_format` flag for specifying merge_SE.tsv format annotations
  
  [Note: This version incorrectly referred to "merge_SV.tsv" instead of "merge_SE.tsv" in some places. This has been corrected in version 1.7.0.]

### Changed
- **SEgene_peakprep workflow optimization**:
  - Streamlined command-line interface with more intuitive parameter names
  - Adjusted log level output for better readability during execution
  - Improved temporary file handling
- **Documentation updates**:
  - Comprehensive revision of README files (both English and Japanese)
  - Added clearer examples for different annotation file formats
  - Updated command examples to reflect the new parameter structure
  - Improved troubleshooting section

## [1.6.0] - 2025-04-06

### Added
- **New component: SEgene_peakprep**:
  - New data preprocessing pipeline for ChIP-seq data normalization
  - Integrated processing of peak region information from multiple BAM files into normalized data tables
  - Calculation of total mapped reads using `samtools flagstat`
  - Read counting in genomic regions using `featureCounts`
  - CPM calculation and log transformation functionality
  - Integration with existing workflow components
- **Documentation**:
  - Documentation for SEgene_peakprep
  - Revised workflow description positioning SEgene_peakprep as the initial data preparation step

### Changed
- Reorganized project structure to include SEgene_peakprep as the first step in the workflow
- Updated SEgene workflow description:
  - SEgene_peakprep → peak_to_gene_links → SE_to_gene_links → SEgene_RegionAnalyzer

## [1.5.0] - 2025-04-05

### Added
- **New component: SEgene_RegionAnalyzer**:
  - Added a new analytical tool for evaluating super-enhancer activity in genomic regions of interest
  - Integrated with public databases (currently SEdb 2.0)
  - Performs enrichment analysis of tissue-specific super-enhancer associations
  - Supports both SEgene output TSV files and standard BED format inputs
  - Generates comprehensive TSV reports and visual HTML outputs
- **Documentation**:
  - Added detailed documentation for SEgene_RegionAnalyzer
  - Updated main README to include the new component in the project structure

### Changed
- Reorganized the project structure to include SEgene_RegionAnalyzer as an optional component
- Updated the workflow description to reflect the extended analytical capabilities

## [1.4.0] - 2025-03-19

### Added
- **Enhanced network visualization capabilities**:
  - Added detailed graph drawing methods (`draw_network_detailed`, `draw_subnetwork_detailed`, `draw_two_layer_subnetwork_detailed`)
  - Added support for multiple output formats (SVG, PNG, EPS, PDF) across all visualization functions
  - Added customization options for node styling, edge styling, and subgraph boundaries
- **Improved data export functionality**:
  - Added ability to save graph data alongside visualizations
  - Added logging capabilities for reproducibility
  - Added support for exporting data in various formats (CSV, TSV, JSON)
- **Enhanced visualization parameters**:
  - Added resolution (DPI) control for raster format outputs
  - Added figure size and styling customization options
  - Added custom title support for all network visualization methods (`draw_network`, `draw_subnetwork`, `draw_two_layer_subnetwork`)

### Changed
- Updated existing visualization methods to support expanded format options
- Extended graph rendering interfaces with additional customization parameters
- Improved file naming and organization for saved outputs

### Notes
- These visualization enhancements are designed to maintain compatibility with existing network analysis functions

## [1.3.0] - 2025-03-15

### Added
- **Enhanced visualization options** for ROSE summary plots:
  - Added customizable parameters for DPI, figure size, and output format
  - Improved figure saving capabilities with format selection

### Fixed
- Fixed threshold comparison in `p2gl_path_to_filter_df` function to correct filtering behavior

### Notes
- Added pronunciation guide for SEgene ("S-E-gene") in documentation

## [1.2.0] - 2025-03-11

### Added
- **New core analysis capabilities** for SE rank distribution:
  - Added functionality to search for gene-linked super-enhancers across multiple samples
  - Implemented super-enhancer ranking analysis within datasets 
  - Developed statistical tools for percentile calculation and distribution analysis
  - Created visualization methods for SE rank distributions across samples
- **New tutorial notebooks** to demonstrate the new analysis features:
  - `tutorial_book_SE_rank_disribution.ipynb` (English)
  - `tutorial_book_SE_rank_disribution_ja.ipynb` (Japanese)
  - These notebooks guide users through analyzing the ranking status (position and percentile) of super-enhancers corresponding to specific genes

### Notes
- Documentation in README files has been updated to include information about the new notebooks

## [1.1.0] - 2025-02-07

### Added
- **New CLI tools (`cli_tools/`)** for streamlined featureCounts processing, CPM calculation, and SE–gene correlation analysis
  - Scripts such as `generate_gtf.py`, `generate_file_list.sh`, `run_featurecounts_array.sh`, etc.
  - `correlation_analysis.py` computes Pearson correlation between SE (CPM) and genes (TPM)
- **Documentation updates**:
  - Added **English (`README.md`)** and **Japanese (`README_ja.md`)** documentation under `cli_tools/`
  - Updated the main repository documentation (`README.md`) and `SE_to_gene_links/` docs to explain how to integrate the new CLI tools

### Changed
- Refined repository structure by separating stand-alone CLI scripts into `cli_tools/`, independent of `SEgene_package`

### Notes
- The newly added documents provide usage instructions and workflow details in both English and Japanese  
- Enables advanced or extended analyses on top of the existing `SE_to_gene_links` functionalities

## [1.0.0] - 2025-01-15

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

[1.12.0]: https://github.com/hamamoto-lab/SEgene/compare/v1.11.0...v1.12.0
[1.11.0]: https://github.com/hamamoto-lab/SEgene/compare/v1.10.0...v1.11.0
[1.10.0]: https://github.com/hamamoto-lab/SEgene/releases/tag/v1.10.0
[1.9.0]: https://github.com/hamamoto-lab/SEgene/releases/tag/v1.9.0
[1.8.0]: https://github.com/hamamoto-lab/SEgene/releases/tag/v1.8.0
[1.7.0]: https://github.com/hamamoto-lab/SEgene/releases/tag/v1.7.0
[1.6.1]: https://github.com/hamamoto-lab/SEgene/releases/tag/v1.6.1
[1.6.0]: https://github.com/hamamoto-lab/SEgene/releases/tag/v1.6.0
[1.5.0]: https://github.com/hamamoto-lab/SEgene/releases/tag/v1.5.0
[1.4.0]: https://github.com/hamamoto-lab/SEgene/releases/tag/v1.4.0
[1.3.0]: https://github.com/hamamoto-lab/SEgene/releases/tag/v1.3.0
[1.2.0]: https://github.com/hamamoto-lab/SEgene/releases/tag/v1.2.0
[1.1.0]: https://github.com/hamamoto-lab/SEgene/releases/tag/v1.1.0
[1.0.0]: https://github.com/hamamoto-lab/SEgene/releases/tag/v1.0.0