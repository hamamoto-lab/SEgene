# Data Directory

This directory is for storing SEdb data files and other datasets.

## Required Files for SEgene_analyzer

Download the following files from [SEdb 2.0](http://www.licpathway.net/sedb/download.php) and place them in the `SEdb/` subdirectory:

### SEdb 2.0 Files
- `human_sample_information_sedb2.txt` - Sample metadata
- `SE_package_hg38.bed` - Super-enhancer definitions for human (hg38)

### Directory Structure
```
data/
└── SEdb/
    ├── human_sample_information_sedb2.txt
    └── SE_package_hg38.bed
```

## Download Instructions

1. Visit: [SEdb 2.0 Download Page](http://www.licpathway.net/sedb/download.php)
2. Download "Sample information" → `human_sample_information_sedb2.txt`
3. Download "Package of SE, SE element and TE" → Human (hg38) → `SE_package_hg38.bed`
4. Place files in the `SEdb/` subdirectory

## Note

- These data files are **not included in the git repository** due to size and licensing considerations
- The `data/` directory is listed in `.gitignore`
- Each user should download their own copy of the data files