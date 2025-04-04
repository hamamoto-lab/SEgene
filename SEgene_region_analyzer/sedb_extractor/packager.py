"""
SEdb Region Data Packaging Module
Provides functionality to save extracted region data in a standardized format.
"""

import os
import json
import pandas as pd
import logging
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, Optional, List, Union

from sedb_common.logging import get_module_logger
from .extractor import RegionOverlapAnalyzer
#from sedb_extractor.extractor import RegionOverlapAnalyzer



class RegionPackager:
    """Class for packaging extracted region data into a standardized format."""
    
    def __init__(self, logger=None):
        """Initialization."""
        self.logger = logger or get_module_logger(__name__)
    
    def save_region_package(self, extraction_results: Dict[str, Any], output_dir: str,
                          include_distributions: bool = True, compress: bool = True,
                          source_files: Optional[Dict[str, Any]] = None,
                          export_bed: bool = False, bed_format: str = 'classic') -> Dict[str, Any]:
        """
        Save extracted region data as a standardized package.
        
        Parameters
        ----------
        extraction_results : dict
            Results from extract_region().
        output_dir : str
            Output directory path.
        include_distributions : bool, optional
            Whether to include distribution data.
        compress : bool, optional
            Whether to compress the data.
        source_files : dict, optional
            Information about source files (paths and MD5 hashes).
        export_bed : bool, optional
            Whether to also export a BED file.
        bed_format : str, optional
            Output BED format ('classic' or 'extended').
            
        Returns
        -------
        dict
            Information about the saved package.
        """
        # ディレクトリ構造を作成
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        # サブディレクトリを作成
        dist_dir = output_path / "distributions"
        if include_distributions:
            dist_dir.mkdir(exist_ok=True)
        
        # 保存されるファイルのリスト
        files = []
        
        # 1. 領域情報を保存
        region_file = output_path / "region_info.json"
        with open(region_file, 'w') as f:
            json.dump(extraction_results['region'], f, indent=2)
        files.append({"name": "region_info", "path": "region_info.json", "type": "json"})
        
        # 2. 重複するSEデータを保存
        if 'overlapping_se' in extraction_results and len(extraction_results['overlapping_se']) > 0:
            bed_file = output_path / "bed_data.parquet"
            extraction_results['overlapping_se'].to_parquet(
                bed_file, 
                compression='snappy' if compress else None
            )
            files.append({"name": "bed_data", "path": "bed_data.parquet", "type": "parquet"})
            
            # BEDファイルも出力する場合
            if export_bed:

                extractor = RegionOverlapAnalyzer(logger=self.logger)
                
                bed_output_path = output_path / f"overlapping_regions.bed"
                bed_file_path = extractor.export_to_bed(
                    extraction_results['overlapping_se'],
                    bed_output_path,
                    format=bed_format
                )
                
                if bed_file_path:
                    files.append({"name": "bed_file", "path": "overlapping_regions.bed", "type": "bed"})
        
        # 3. サンプル情報を保存
        if 'sample_summary' in extraction_results and len(extraction_results['sample_summary']) > 0:
            sample_file = output_path / "sample_info.parquet"
            extraction_results['sample_summary'].to_parquet(
                sample_file, 
                compression='snappy' if compress else None
            )
            files.append({"name": "sample_info", "path": "sample_info.parquet", "type": "parquet"})
        
        # 4. 分布データを保存（オプション）
        if include_distributions:
            # 組織分布
            if 'tissue_counts' in extraction_results and extraction_results['tissue_counts'] is not None:
                tissue_df = pd.DataFrame({
                    'tissue_type': extraction_results['tissue_counts'].index,
                    'count': extraction_results['tissue_counts'].values
                })
                tissue_file = dist_dir / "tissue_counts.parquet"
                tissue_df.to_parquet(tissue_file, compression='snappy' if compress else None)
                files.append({"name": "tissue_counts", "path": "distributions/tissue_counts.parquet", "type": "parquet"})
            
            # バイオサンプル分布
            if 'biosample_counts' in extraction_results and extraction_results['biosample_counts'] is not None:
                biosample_df = pd.DataFrame({
                    'biosample_type': extraction_results['biosample_counts'].index,
                    'count': extraction_results['biosample_counts'].values
                })
                biosample_file = dist_dir / "biosample_counts.parquet"
                biosample_df.to_parquet(biosample_file, compression='snappy' if compress else None)
                files.append({"name": "biosample_counts", "path": "distributions/biosample_counts.parquet", "type": "parquet"})
            
            # 遺伝子分布
            if 'gene_counts' in extraction_results and extraction_results['gene_counts'] is not None:
                gene_df = pd.DataFrame({
                    'gene': extraction_results['gene_counts'].index,
                    'count': extraction_results['gene_counts'].values
                })
                gene_file = dist_dir / "gene_counts.parquet"
                gene_df.to_parquet(gene_file, compression='snappy' if compress else None)
                files.append({"name": "gene_counts", "path": "distributions/gene_counts.parquet", "type": "parquet"})
        
        # 5. サマリー統計を保存
        summary_file = output_path / "summary.json"
        
        # サマリー統計を拡張
        enhanced_stats = extraction_results['stats'].copy()
        
        # トップ組織を追加
        if 'tissue_counts' in extraction_results and extraction_results['tissue_counts'] is not None:
            top_tissues = []
            for tissue, count in extraction_results['tissue_counts'].nlargest(5).items():
                top_tissues.append({"tissue": tissue, "count": int(count)})
            enhanced_stats['top_tissues'] = top_tissues
        
        # トップバイオサンプルを追加
        if 'biosample_counts' in extraction_results and extraction_results['biosample_counts'] is not None:
            top_biosamples = []
            for biosample, count in extraction_results['biosample_counts'].nlargest(5).items():
                top_biosamples.append({"biosample": biosample, "count": int(count)})
            enhanced_stats['top_biosamples'] = top_biosamples
        
        # トップ遺伝子を追加
        if 'gene_counts' in extraction_results and extraction_results['gene_counts'] is not None:
            top_genes = []
            for gene, count in extraction_results['gene_counts'].nlargest(5).items():
                top_genes.append({"gene": gene, "count": int(count)})
            enhanced_stats['top_genes'] = top_genes
        
        with open(summary_file, 'w') as f:
            json.dump(enhanced_stats, f, indent=2)
        files.append({"name": "summary", "path": "summary.json", "type": "json"})
        
        # 6. メタデータを作成して保存
        metadata = {
            "version": "1.0",
            "created_at": datetime.now().isoformat(),
            "region": extraction_results['region'],
            "extraction_method": extraction_results.get('extraction_method', 'unknown'),
            "stats": extraction_results['stats'],
            "files": files
        }
        
        # ソースファイル情報を追加
        if source_files:
            metadata["source_files"] = source_files
        
        metadata_file = output_path / "metadata.json"
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)
        
        self.logger.info(f"Region data package saved: {output_dir}")
        self.logger.info(f"Number of files saved: {len(files)}")
        
        return {
            "output_dir": str(output_path),
            "metadata": metadata,
            "files": files
        }
