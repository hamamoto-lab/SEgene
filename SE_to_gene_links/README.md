# SEgene_package

*(For the Japanese version of this README, please see [README_ja.md](https://github.com/hamamoto-lab/SEgene/blob/main/SE_to_gene_links/README_ja.md).)*

SEgene_package is a Python package for analyzing the relationship between super-enhancers and genes. This repository provides the necessary setup instructions for both Docker and non-Docker users.

## Table of Contents

- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
  - [For Docker Users](#for-docker-users)
  - [For Non-Docker Users](#for-non-docker-users)
- [Usage](#usage)
- [License](#license)
- [Citation](#Citation)

## Features

- Analyze the relationship between super-enhancers and genes
- Visualize data using graph theory
- Interactive analysis with Jupyter Notebook

## Requirements

This system can be run in a Docker environment or locally using miniforge3.
For Docker usage, the system can run with Docker alone, but Docker Compose V2 is recommended for easier operation.

To check if Docker Compose V2 is installed, run:

```bash
docker compose version
```

You should see a message indicating version 2.

**Note:** The command is `docker compose`, not `docker-compose`.

For installation on Ubuntu, run:

```bash
sudo apt-get update
sudo apt-get install docker-compose-plugin
```

Refer to the [official documentation](https://docs.docker.com/compose/install/linux/) for more details.

For local installation, install miniforge3 on your system.

## Installation Guide

**Common Setup**

Prepare the following directories in your environment:

- A directory for storing analysis data: `/path/to/data`
- A directory for saving intermediate data and analysis results: `/path/to/work`

If running with Docker, these directories must be bound to Docker beforehand.

### For Docker Users

**Steps for starting with Docker Compose**

1. **Clone the repository**

    ```bash
    git clone https://github.com/hamamoto-lab/SEgene.git
    ```

2. **Navigate to the directory**

    ```bash
    cd SEgene/SE_to_gene_links
    ```

3. **Edit docker-compose.yml**

    Open `docker-compose.yml` in an editor and update the following lines:

    ```yaml
      # Replace "/path/to/data" with your "data" folder
      - /path/to/data:/home/jovyan/data
      # Replace "/path/to/work" with your "work" folder
      - /path/to/work:/home/jovyan/work
    ```

4. **Build the Docker image**

    ```bash
    docker-compose build
    ```

5. **Start the container**

    ```bash
    docker-compose up -d
    ```

6. **Access Jupyter Notebook**

    Open your browser and go to:

    ```
    http://localhost:8888
    ```

7. **Use the notebooks**

    Open the notebooks in the `notebooks/` directory to start your analysis.

8. **Stop the container**

    When you're done, stop the container using:

    ```bash
    docker-compose down
    ```

### For Non-Docker Users

1. **Clone the repository**

    ```bash
    git clone https://github.com/hamamoto-lab/SEgene.git
    ```

2. **Navigate to the directory**

    ```bash
    cd SEgene/SE_to_gene_links
    ```

3. **Create a Conda environment**

    Use `environment.yml` to create the environment.

    ```bash
    conda env create -f environment.yml
    ```

4. **Activate the environment**

    ```bash
    conda activate se_gene
    ```

5. **Install the package**

    Install the package.

    ```bash
    pip install .
    ```

6. **Set up IPython startup script (optional)**

    ```bash
    mkdir -p ~/.ipython/profile_default/startup/
    cp docker/startup/00-imports.py ~/.ipython/profile_default/startup/
    ```

7. **Create a data directory (if necessary)**

    ```bash
    mkdir -p ~/data
    ```

8. **Start Jupyter Notebook**

    ```bash
    jupyter notebook
    ```

## Usage

- **Using the notebooks**: 
  - [Main analysis notebook (including network analysis)](https://github.com/hamamoto-lab/SEgene/blob/main/SE_to_gene_links/notebooks/tutorial_book.ipynb): Comprehensive tutorial covering the core analytical functions and network analysis capabilities.
  - [SE rank distribution analysis notebook](https://github.com/hamamoto-lab/SEgene/blob/main/SE_to_gene_links/notebooks/tutorial_book_SE_rank_disribution.ipynb): Added in v1.2.0, this notebook demonstrates how to analyze the ranking status (position and percentile in SE lists) of super-enhancers corresponding to specific genes within SE file datasets.
- **Package features**: SEgene_package includes all tools necessary for analyzing the relationship between super-enhancers and genes.
- **Additional Analysis Tools**: We also provide [command-line tools](https://github.com/hamamoto-lab/SEgene/tree/main/cli_tools) for analyzing correlations between SE regions and gene expression. Please refer to each tool's documentation for details.

## License

This project is licensed under the MIT License. See the [`LICENSE`](https://github.com/hamamoto-lab/SEgene/blob/main/LICENSE) for details.

## Citation

If you use this tool for your research, please refer to the following CITATION file:
[CITATION](https://github.com/hamamoto-lab/SEgene/blob/main/CITATION)
