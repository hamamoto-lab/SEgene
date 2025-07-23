#!/usr/bin/env python3
"""
Create sample data for testing SEgene_analyzer_erna

This script generates synthetic eRNAbase-like data for testing and development purposes.
"""

import os
import pandas as pd
import numpy as np
from pathlib import Path
import argparse

def create_sample_bed_files(output_dir, num_samples=5, peaks_per_sample=100):
    """
    Create sample BED files mimicking eRNAbase peak data
    
    Parameters:
    -----------
    output_dir : str
        Directory to save BED files
    num_samples : int
        Number of sample BED files to create
    peaks_per_sample : int
        Number of peaks per BED file
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Define chromosomes and their approximate lengths (simplified)
    chromosomes = {
        'chr1': 249250000, 'chr2': 243200000, 'chr3': 198300000,
        'chr4': 191200000, 'chr5': 180900000, 'chr6': 171100000,
        'chr7': 159100000, 'chr8': 146300000, 'chr9': 141200000,
        'chr10': 135500000, 'chr11': 135000000, 'chr12': 133800000
    }
    
    sample_metadata = []
    
    for i in range(num_samples):
        sample_id = f"sample_{i+1:03d}"
        bed_file = os.path.join(output_dir, f"{sample_id}.bed")
        
        # Generate random peaks
        peaks = []
        for j in range(peaks_per_sample):
            # Randomly select chromosome
            chrom = np.random.choice(list(chromosomes.keys()))
            chrom_length = chromosomes[chrom]
            
            # Generate random start position
            start = np.random.randint(1000000, chrom_length - 2000)  # Leave some buffer
            end = start + np.random.randint(500, 2000)  # Peak length 500-2000bp
            
            # Create peak ID
            peak_id = f"peak_{j+1:04d}"
            
            # Additional columns (simplified eRNAbase format)
            fow_start2 = start - np.random.randint(0, 500)
            fow_end2 = end + np.random.randint(0, 500)
            size = end - start
            
            # Binary features (0 or 1)
            risk_snp = np.random.choice([0, 1], p=[0.9, 0.1])
            common_snp = np.random.choice([0, 1], p=[0.7, 0.3])
            crisps = np.random.choice([0, 1], p=[0.8, 0.2])
            eqtl = np.random.choice([0, 1], p=[0.6, 0.4])
            tfbs = np.random.choice([0, 1], p=[0.5, 0.5])
            enhancer = np.random.choice([0, 1], p=[0.3, 0.7])
            super_enhancer = np.random.choice([0, 1], p=[0.9, 0.1])
            
            peak = [
                chrom, start, end, peak_id, fow_start2, fow_end2, size,
                risk_snp, common_snp, crisps, eqtl, tfbs, enhancer, super_enhancer
            ]
            peaks.append(peak)
        
        # Write BED file
        with open(bed_file, 'w') as f:
            for peak in peaks:
                f.write('\t'.join(map(str, peak)) + '\n')
        
        # Collect sample metadata
        sample_metadata.append({
            'sample_id': sample_id,
            'bed_file': f"{sample_id}.bed"
        })
        
        print(f"Created BED file: {bed_file} ({peaks_per_sample} peaks)")
    
    return sample_metadata

def create_sample_metadata(output_file, sample_metadata, extended=True):
    """
    Create sample metadata file mimicking eRNAbase metadata
    
    Parameters:
    -----------
    output_file : str
        Path to save metadata file
    sample_metadata : list
        List of sample metadata dictionaries
    extended : bool
        Whether to include extended metadata columns
    """
    # Define possible values for metadata fields
    tissues = ['brain', 'liver', 'heart', 'lung', 'kidney', 'muscle', 'skin', 'blood']
    cell_types = [
        'neuron', 'hepatocyte', 'cardiomyocyte', 'pneumocyte', 
        'nephron', 'myocyte', 'keratinocyte', 'lymphocyte'
    ]
    
    # Create enhanced metadata
    for i, sample in enumerate(sample_metadata):
        sample.update({
            'Species': 'Homo sapiens',  # Match actual eRNAbase format
            'Tissue Type': np.random.choice(tissues),  # Match expected column name
            'Cell Type': np.random.choice(cell_types),  # Match expected column name
            'Biosample Type': np.random.choice(['Primary cells', 'Cell line', 'Tissue']),
            'Experiment Type': np.random.choice(['ChIP-seq', 'ATAC-seq', 'DNase-seq']),
            'source': f'lab_{np.random.randint(1, 6)}',
        })
        
        if extended:
            sample.update({
                'gender': np.random.choice(['male', 'female']),
                'age_group': np.random.choice(['young', 'adult', 'elderly']),
                'disease_status': np.random.choice(['healthy', 'diseased'], p=[0.7, 0.3]),
                'batch': f'batch_{np.random.randint(1, 4)}',
                'sequencing_depth': np.random.randint(10000000, 50000000),
                'peak_count': np.random.randint(10000, 50000),
                'quality_score': np.random.uniform(0.7, 1.0)
            })
    
    # Create DataFrame and save
    metadata_df = pd.DataFrame(sample_metadata)
    
    # Save as both CSV and Parquet for testing both formats
    csv_file = output_file.replace('.parquet', '.csv')
    metadata_df.to_csv(csv_file, index=False)
    
    if output_file.endswith('.parquet'):
        metadata_df.to_parquet(output_file, index=False)
        print(f"Created metadata files: {csv_file} and {output_file}")
    else:
        print(f"Created metadata file: {output_file}")
    
    return metadata_df

def create_sample_regions(output_file, num_regions=20):
    """
    Create sample regions TSV file for testing
    
    Parameters:
    -----------
    output_file : str
        Path to save regions file
    num_regions : int
        Number of regions to generate
    """
    regions = []
    
    chromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5']
    
    for i in range(num_regions):
        chrom = np.random.choice(chromosomes)
        start = np.random.randint(1000000, 200000000)
        end = start + np.random.randint(100000, 2000000)  # Region size 100kb-2Mb
        
        # Create region identifier
        region_id = f"{chrom}_{start}_{end}"
        
        regions.append({
            'region_id': region_id,
            'chromosome': chrom,
            'start_position': start,
            'end_position': end,
            'description': f'Test region {i+1}',
            'priority': np.random.choice(['high', 'medium', 'low'])
        })
    
    regions_df = pd.DataFrame(regions)
    regions_df.to_csv(output_file, sep='\t', index=False)
    
    print(f"Created regions file: {output_file} ({num_regions} regions)")
    return regions_df

def main():
    """Main function to create all sample data"""
    parser = argparse.ArgumentParser(
        description='Create sample data for SEgene_analyzer_erna testing',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('--output-dir', default='sample_data',
                       help='Output directory for sample data')
    parser.add_argument('--num-samples', type=int, default=10,
                       help='Number of sample BED files to create')
    parser.add_argument('--peaks-per-sample', type=int, default=100,
                       help='Number of peaks per BED file')
    parser.add_argument('--num-regions', type=int, default=20,
                       help='Number of test regions to create')
    parser.add_argument('--extended-metadata', action='store_true',
                       help='Include extended metadata columns')
    
    args = parser.parse_args()
    
    # Create main output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    print(f"Creating sample data in: {args.output_dir}")
    print(f"Samples: {args.num_samples}, Peaks per sample: {args.peaks_per_sample}")
    print(f"Test regions: {args.num_regions}")
    
    # Create directory structure
    bed_dir = os.path.join(args.output_dir, 'eRNAbase', 'peak')
    regions_dir = os.path.join(args.output_dir, 'regions')
    
    os.makedirs(bed_dir, exist_ok=True)
    os.makedirs(regions_dir, exist_ok=True)
    
    print("\n1. Creating sample BED files...")
    sample_metadata = create_sample_bed_files(
        bed_dir, 
        args.num_samples, 
        args.peaks_per_sample
    )
    
    print("\n2. Creating sample metadata...")
    metadata_file = os.path.join(args.output_dir, 'eRNAbase', 'eRNAbase_data.parquet')
    metadata_df = create_sample_metadata(
        metadata_file, 
        sample_metadata, 
        args.extended_metadata
    )
    
    print("\n3. Creating sample regions...")
    regions_file = os.path.join(regions_dir, 'test_regions.tsv')
    regions_df = create_sample_regions(regions_file, args.num_regions)
    
    print("\nSample data creation completed!")
    print("\nData structure:")
    print(f"  {args.output_dir}/")
    print(f"    ├── eRNAbase/")
    print(f"    │   ├── peak/              # {args.num_samples} BED files")
    print(f"    │   ├── eRNAbase_data.csv  # Metadata (CSV)")
    print(f"    │   └── eRNAbase_data.parquet  # Metadata (Parquet)")
    print(f"    └── regions/")
    print(f"        └── test_regions.tsv   # {args.num_regions} test regions")
    
    print("\nExample usage:")
    print(f"  python examples/production_workflow.py \\")
    print(f"    --bed-dir {bed_dir} \\")
    print(f"    --metadata {metadata_file} \\")
    print(f"    --regions {regions_file} \\")
    print(f"    --output-dir results \\")
    print(f"    --max-rows 5")

if __name__ == "__main__":
    main()