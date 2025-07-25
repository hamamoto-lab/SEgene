# eRNAbase Metadata

eRNAbase metadata files required for SEgene_analyzer_erna analysis.

**Data Source**: [eRNAbase](https://bio.liclab.net/eRNAbase/index.php) - Enhancer RNA database  
**Data Version**: This metadata corresponds to eRNAbase data as of **July 22, 2025**.

## Files

- **eRNAbase_data.csv**: CSV format
- **eRNAbase_data.parquet**: Parquet format

## Data Statistics

- **Total samples**: 1,012
- **Human samples**: 858, **Mouse samples**: 154
- **Tissue types**: 43, **Cell types**: 119

## Usage

```bash
# Specify metadata file for analysis
erna-analyzer single -m metadata/eRNAbase_data.parquet --chr chr7 --start 1000000 --end 2000000
```

## Data Structure

| Column | Description |
|--------|-------------|
| Sample ID | Sample identifier |
| Species | Species (Homo sapiens / Mus musculus) |
| Tissue Type | Tissue type |
| Cell Type | Cell type |
| Biosample Type | Sample type |
| Experiment Type | Experiment type |