
# P2GL: Peak-to-Gene Links Tool

*(For the Japanese version of this README, please see [README_ja.md](https://github.com/hamamoto-lab/SEgene/blob/main/peak_to_gene_links/README_ja.md).)*

This program is a tool for constructing links from peaks to genes using ChIP-seq peak data and RNA-seq expression data. It operates in **R 4.2.2** and **Julia 1.8.3** environments. A Docker container is also provided for simplified setup.

---

## Requirements

### R Libraries
The following R packages are required:
- CRAN packages:
  - `BiocManager`
  - `data.table`
  - `openxlsx`
  - `optparse`
  - `pbmcapply`
  - `stringr`
- Bioconductor packages:
  - `GenomicRanges`
  - `rhdf5`

Installation commands:
```r
# CRAN packages
install.packages(c("BiocManager", "data.table", "openxlsx", "optparse", "pbmcapply", "stringr"))

# Bioconductor packages
BiocManager::install(c("GenomicRanges", "rhdf5"))
```

### Julia Libraries
The following Julia packages are required:
- `ArgParse`
- `HDF5`
- `RData`
- `StatsBase`

Installation commands:
```julia
using Pkg
Pkg.add(["ArgParse", "HDF5", "RData", "StatsBase"])
```

---

## Quick Start with Docker

This program runs inside a Docker container. Follow the steps below to run the test.

### 1. Download the Test Data
Please download the following test files from Figshare: [Figshare Link](#)  
*Note: The link will be made available once our SEgene paper has been published.*
- `GSE156614_rna_tumor.csv`
- `GSE156614_ChIP_tumor.tsv`

Create the following local directories:
- **Data directory**: `/path/to/data`
- **Output directory**: `/path/to/output`

Place the downloaded files into `/path/to/data`.

### 2. Clone the Repository
Clone the repository and navigate to the appropriate directory:
```bash
git clone https://github.com/your-repo-name/P2GL.git
cd P2GL/SEgene_test/peak_to_gene_links
```

### 3. Build the Docker Image
Build the Docker image using the provided Dockerfile:
```bash
docker build -t p2gl:latest -f docker/Dockerfile .
```

### 4. Run the Program inside Docker
Start the Docker container and execute the program:

1. Enter the Docker container:
```bash
docker run -it --rm   -v /path/to/data:/data   -v /path/to/output:/opt/P2GL/output   p2gl:latest
```

2. Inside the container, execute the P2GL program:
```bash
Rscript peakToGeneLinks.R   -E /data/GSE156614_rna_tumor.csv   -P /data/GSE156614_ChIP_tumor.tsv   -O sample_2000000   --cores 1   --txWidth 2000000
```

### 5. Check the Output
Upon successful execution, the following files will be generated in `/path/to/output`:
- `sample_2000000.tsv`
- `sample_2000000.RData`

### 6. (Optional) Generate Visualization PDFs
To generate visualization PDFs based on the results, execute the following:
```bash
Rscript peakToGeneLinksPlot.R   -I /opt/P2GL/output/sample_2000000   -O sample_2000000_plot.pdf   -x 1:1000
```

This will produce `sample_2000000_plot.pdf` in the local `/path/to/output` folder.

### 7. Exit Docker
After completing the tasks, exit the Docker container:
```bash
exit
```

---

## Notes
- Ensure Docker is installed and properly configured on your system.
- Modify `/path/to/data` and `/path/to/output` to match your local environment.
- The `sample_2000000.tsv` file generated will be similar to the one on [Figshare](#). However, random seeds can affect the output slightly. To fix the output in your environment, use the `--seed` option. Note that the seed's effect may vary across different systems.
*Note: The link will be made available once our SEgene paper has been published.*

---

## Citation

Please refer to the [CITATION](https://github.com/hamamoto-lab/SEgene/blob/main/CITATION) of SEgene.


## License

Please refer to the [LICENSE](https://github.com/hamamoto-lab/SEgene/blob/main/LICENSE) file of SEgene.
