#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
sedb_analyzer.packager - Packaging and Saving of SEdb Analysis Results

This module provides functionalities to save statistical analysis results to the file system and to generate a concise report.
"""

import os
import json
import pandas as pd
import logging
import warnings
import numpy as np 
import sys      
# import pprint   
from pathlib import Path
from typing import Dict, List, Any, Optional, Union, Tuple
from datetime import datetime

# ヘルパーモジュールをインポート
from . import report_utils
from . import html_generator

# Numpy のデータ型を JSON シリアライズ可能にするためのエンコーダー
class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            # Convert NaN or Infinity to None (JSON null)
            if np.isnan(obj) or np.isinf(obj):
                return None
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, pd.Timestamp): # pandas の Timestamp も考慮
            return obj.isoformat()
        # Convert pd.NA to None (JSON null)
        if obj is pd.NA:
            return None
        return super(NpEncoder, self).default(obj)

class AnalysisPackager:
    """
    A class for packaging and saving SEdb analysis results into a specified directory structure.
    """

    def __init__(self, logger=None):
        """
        Initialization
        ----------
        logger : Logger, optional
            The logger instance to be used.
        """
        self.logger = logger or self._setup_default_logger()

    def _setup_default_logger(self):
        """Setup default logger."""
        logger = logging.getLogger(f"{__name__}.{id(self)}")
        if not logger.handlers:
            logger.setLevel(logging.INFO)
            handler = logging.StreamHandler()
            formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
            handler.setFormatter(formatter)
            logger.addHandler(handler)
            logger.propagate = False
        return logger

    def save_analysis_package(self, output_dir: str,
                           region_info: dict,
                           enrichment_results: Optional[pd.DataFrame] = None,
                           sample_analysis: Optional[Dict[str, Any]] = None,
                           gene_analysis: Optional[Dict[str, Any]] = None,
                           figures: Optional[Dict[str, Any]] = None,
                           metadata: Optional[Dict[str, Any]] = None,
                           generate_html: bool = True,
                           image_formats: List[str] = ['png', 'svg'],
                           data_formats_machine: List[str] = ['parquet', 'json'],
                           data_formats_human: List[str] = ['csv', 'txt']) -> Dict[str, Any]:
        """
        Organize analysis results and save them into a directory structure.
        [Parameter descriptions omitted]
        Returns
        -------
        Dict[str, Any]
            A dictionary including paths of saved files and related information.
        """
        # --- Create output directory structure ---
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        machine_dir = output_path / "machine_readable"
        human_dir = output_path / "human_readable"
        machine_data_dir = machine_dir / "data"
        human_tables_dir = human_dir / "tables"
        human_figures_dir = human_dir / "figures"
        machine_dir.mkdir(exist_ok=True)
        human_dir.mkdir(exist_ok=True)
        machine_data_dir.mkdir(exist_ok=True)
        human_tables_dir.mkdir(exist_ok=True)
        human_figures_dir.mkdir(exist_ok=True)

        saved_files = {
            'machine_data': {},
            'human_tables': {},
            'human_figures': {}
        }

        # --- 1. Save region information ---
        if 'json' in data_formats_machine:
            region_file_machine = machine_data_dir / "region_info.json"
            try: # Added try-except for file saving
                 with open(region_file_machine, 'w', encoding='utf-8') as f:
                      json.dump(region_info, f, indent=2)
                 saved_files['machine_data']['region_info_json'] = str(region_file_machine)
            except Exception as e:
                 self.logger.error(f"Failed to save region information JSON ({region_file_machine}): {e}", exc_info=True)

        if 'txt' in data_formats_human:
            region_file_human = human_tables_dir / "region_info.txt"
            try:
                 with open(region_file_human, 'w', encoding='utf-8') as f: # specify encoding
                      f.write(report_utils.dict_to_human_readable(region_info))
                 saved_files['human_tables']['region_info_txt'] = str(region_file_human)
            except Exception as e:
                 self.logger.error(f"Failed to save region information text ({region_file_human}): {e}", exc_info=True)

        # --- 2. Save enrichment results (if available) ---
        if enrichment_results is not None and not enrichment_results.empty:
            if 'parquet' in data_formats_machine:
                enrichment_file_machine = machine_data_dir / "enrichment_results.parquet"
                try:
                    enrichment_results.to_parquet(enrichment_file_machine, index=False)
                    saved_files['machine_data']['enrichment_results_parquet'] = str(enrichment_file_machine)
                except Exception as e:
                    self.logger.warning(f"Failed to save enrichment results in Parquet format: {e}")

            if 'json' in data_formats_machine:
                enrichment_file_machine_json = machine_data_dir / "enrichment_results.json"
                try:
                    # Using custom encoder for safe JSON saving
                    enrichment_results.to_json(enrichment_file_machine_json, orient='records', indent=2, default_handler=str)
                    saved_files['machine_data']['enrichment_results_json'] = str(enrichment_file_machine_json)
                except Exception as e:
                    self.logger.warning(f"Failed to save enrichment results in JSON format: {e}")

            if 'csv' in data_formats_human:
                enrichment_file_human = human_tables_dir / "enrichment_results.csv"
                try:
                    enrichment_results.to_csv(enrichment_file_human, index=False)
                    saved_files['human_tables']['enrichment_results_csv'] = str(enrichment_file_human)
                except Exception as e:
                    self.logger.warning(f"Failed to save enrichment results in CSV format: {e}")

            if 'txt' in data_formats_human:
                enrichment_file_human_txt = human_tables_dir / "enrichment_results.txt"
                try:
                     with open(enrichment_file_human_txt, 'w', encoding='utf-8') as f: # specify encoding
                          f.write(report_utils.dataframe_to_human_readable(enrichment_results))
                     saved_files['human_tables']['enrichment_results_txt'] = str(enrichment_file_human_txt)
                except Exception as e:
                     self.logger.error(f"Failed to save enrichment results text ({enrichment_file_human_txt}): {e}", exc_info=True)


        # --- 3. Save sample analysis results (if available) ---
        if sample_analysis:
            if 'json' in data_formats_machine:
                sample_file_machine = machine_data_dir / "sample_analysis.json"
                try:
                     # Convert Series objects to JSON-safe format
                     json_safe_sample = {}
                     for key, value in sample_analysis.items():
                          if isinstance(value, pd.Series):
                                json_safe_sample[key] = value.to_dict()
                          else:
                                json_safe_sample[key] = value
                     with open(sample_file_machine, 'w', encoding='utf-8') as f: # specify encoding
                          json.dump(json_safe_sample, f, indent=2, cls=NpEncoder)
                     saved_files['machine_data']['sample_analysis_json'] = str(sample_file_machine)
                except Exception as e:
                     self.logger.error(f"Failed to save sample analysis JSON ({sample_file_machine}): {e}", exc_info=True)

            if 'txt' in data_formats_human:
                sample_file_human = human_tables_dir / "sample_analysis.txt"
                try:
                     with open(sample_file_human, 'w', encoding='utf-8') as f: # specify encoding
                          f.write(report_utils.dict_to_human_readable(sample_analysis))
                     saved_files['human_tables']['sample_analysis_txt'] = str(sample_file_human)
                except Exception as e:
                     self.logger.error(f"Failed to save sample analysis text ({sample_file_human}): {e}", exc_info=True)

            # Save tissue distribution
            if 'tissue_counts' in sample_analysis and isinstance(sample_analysis['tissue_counts'], pd.Series):
                tissue_counts = sample_analysis['tissue_counts']
                if not tissue_counts.empty:
                    tissue_df = pd.DataFrame({'tissue_type': tissue_counts.index, 'count': tissue_counts.values})
                    if 'parquet' in data_formats_machine:
                        tissue_file_machine = machine_data_dir / "tissue_distribution.parquet"
                        try:
                             tissue_df.to_parquet(tissue_file_machine, index=False)
                             saved_files['machine_data']['tissue_distribution_parquet'] = str(tissue_file_machine)
                        except Exception as e:
                             self.logger.warning(f"Failed to save tissue distribution in Parquet format: {e}")
                    if 'csv' in data_formats_human:
                        tissue_file_human = human_tables_dir / "tissue_distribution.csv"
                        try:
                             tissue_df.to_csv(tissue_file_human, index=False)
                             saved_files['human_tables']['tissue_distribution_csv'] = str(tissue_file_human)
                        except Exception as e:
                             self.logger.warning(f"Failed to save tissue distribution in CSV format: {e}")


        # --- 4. Save gene analysis results (if available) ---
        if gene_analysis:
            if 'json' in data_formats_machine:
                gene_file_machine = machine_data_dir / "gene_analysis.json"
                try:
                     # Convert Series objects to JSON-safe format
                     json_safe_gene = {}
                     for key, value in gene_analysis.items():
                          if isinstance(value, pd.Series):
                                json_safe_gene[key] = value.to_dict()
                          else:
                                json_safe_gene[key] = value
                     with open(gene_file_machine, 'w', encoding='utf-8') as f: # specify encoding
                          json.dump(json_safe_gene, f, indent=2, cls=NpEncoder)
                     saved_files['machine_data']['gene_analysis_json'] = str(gene_file_machine)
                except Exception as e:
                     self.logger.error(f"Failed to save gene analysis JSON ({gene_file_machine}): {e}", exc_info=True)

            if 'txt' in data_formats_human:
                gene_file_human = human_tables_dir / "gene_analysis.txt"
                try:
                     with open(gene_file_human, 'w', encoding='utf-8') as f: # specify encoding
                          f.write(report_utils.dict_to_human_readable(gene_analysis))
                     saved_files['human_tables']['gene_analysis_txt'] = str(gene_file_human)
                except Exception as e:
                     self.logger.error(f"Failed to save gene analysis text ({gene_file_human}): {e}", exc_info=True)

            # Save gene distribution
            if 'top_genes' in gene_analysis and isinstance(gene_analysis['top_genes'], pd.Series):
                top_genes = gene_analysis['top_genes']
                if not top_genes.empty:
                    gene_df = pd.DataFrame({'gene': top_genes.index, 'count': top_genes.values})
                    if 'parquet' in data_formats_machine:
                        gene_file_machine = machine_data_dir / "gene_distribution.parquet"
                        try:
                             gene_df.to_parquet(gene_file_machine, index=False)
                             saved_files['machine_data']['gene_distribution_parquet'] = str(gene_file_machine)
                        except Exception as e:
                             self.logger.warning(f"Failed to save gene distribution in Parquet format: {e}")
                    if 'csv' in data_formats_human:
                        gene_file_human = human_tables_dir / "gene_distribution.csv"
                        try:
                             gene_df.to_csv(gene_file_human, index=False)
                             saved_files['human_tables']['gene_distribution_csv'] = str(gene_file_human)
                        except Exception as e:
                             self.logger.warning(f"Failed to save gene distribution in CSV format: {e}")

        # --- 5. Save figures ---
        if figures is not None:
            for fig_name, fig in figures.items():
                try:
                     base_path = human_figures_dir / fig_name
                     saved_file_paths = report_utils.save_figure(fig, base_path, formats=image_formats)
                     for file_path in saved_file_paths:
                          fmt = Path(file_path).suffix[1:]
                          saved_files['human_figures'][f"{fig_name}_{fmt}"] = str(file_path)
                          self.logger.debug(f"Figure saved: {file_path}")
                except Exception as e:
                     self.logger.error(f"Error occurred while saving figure ({fig_name}): {e}", exc_info=True)

        # --- 6. Prepare and save summary statistics (revised) ---
        summary_json_file = output_path / "summary.json" # variable name clarified for JSON
        enhanced_stats = {} # initialized as empty
        try:
            # --- Build summary statistics ---
            base_stats = {}
            if sample_analysis:
                # sample_count is provided from analyzer/statistics.py -> passed here
                base_stats['linked_sample_count'] = sample_analysis.get('sample_count', 0) # number of linked samples
                se_stats_val = sample_analysis.get('se_count_stats')
                base_stats['se_count_stats'] = se_stats_val if isinstance(se_stats_val, dict) else {}
                tc_series = sample_analysis.get('tissue_counts')
                base_stats['tissue_type_count'] = len(tc_series) if isinstance(tc_series, pd.Series) else 0
            if gene_analysis:
                base_stats['gene_count'] = gene_analysis.get('gene_count', 0)
                occ_stats_val = gene_analysis.get('occurrence_stats')
                base_stats['occurrence_stats'] = occ_stats_val if isinstance(occ_stats_val, dict) else {}
            if enrichment_results is not None and not enrichment_results.empty:
                if 'Significant' in enrichment_results.columns:
                    base_stats['significant_enrichment_count'] = int(enrichment_results['Significant'].astype(bool).sum())
                if 'Depleted' in enrichment_results.columns:
                    base_stats['significant_depletion_count'] = int(enrichment_results['Depleted'].astype(bool).sum())

            enhanced_stats = base_stats.copy()

            # --- Add top information ---
            if sample_analysis:
                 # Top tissues
                 tc_series = sample_analysis.get('tissue_counts')
                 top_tissues = []
                 if isinstance(tc_series, pd.Series) and not tc_series.empty:
                      for tissue, count in tc_series.nlargest(5).items():
                           top_tissues.append({"tissue": tissue, "count": int(count)})
                 enhanced_stats['top_tissues'] = top_tissues
                 # Top biosamples
                 bc_series = sample_analysis.get('biosample_counts')
                 top_biosamples = []
                 if isinstance(bc_series, pd.Series) and not bc_series.empty:
                      for biosample, count in bc_series.nlargest(5).items():
                           top_biosamples.append({"biosample": biosample, "count": int(count)})
                 enhanced_stats['top_biosamples'] = top_biosamples

            if gene_analysis:
                 # Top genes
                 tg_series = gene_analysis.get('top_genes')
                 top_genes = []
                 if isinstance(tg_series, pd.Series) and not tg_series.empty:
                      for gene, count in tg_series.nlargest(5).items():
                           top_genes.append({"gene": gene, "count": int(count)})
                 enhanced_stats['top_genes'] = top_genes

            # --- Write to JSON file ---
            self.logger.debug(f"Attempting to save summary stats to {summary_json_file}")

            try:
                with open(summary_json_file, 'w', encoding='utf-8') as f:
                    json.dump(enhanced_stats, f, indent=2, cls=NpEncoder)
                self.logger.info(f"Summary statistics (JSON) saved: {summary_json_file}")
                saved_files['summary_json'] = str(summary_json_file)
            except Exception as e:
                self.logger.error(f"Error occurred while saving summary statistics JSON ({summary_json_file}): {e}", exc_info=True)

        except Exception as e:
             self.logger.error(f"Error occurred while building summary statistics data: {e}", exc_info=True)

        # --- 7. Generate summary text ---
        summary_txt_file = human_dir / "summary.txt" # variable name changed
        try:
             summary_text = report_utils.generate_analysis_summary(
                  region_info=region_info,
                  enrichment_results=enrichment_results,
                  sample_analysis=sample_analysis,
                  gene_analysis=gene_analysis
             )
             with open(summary_txt_file, 'w', encoding='utf-8') as f: # specify encoding
                  f.write(summary_text)
             saved_files['summary_txt'] = str(summary_txt_file)
             self.logger.info(f"Summary text generated: {summary_txt_file}")
        except Exception as e:
             self.logger.error(f"Failed to generate or save summary text ({summary_txt_file}): {e}", exc_info=True)

        # --- 8. Generate HTML report (optional) ---
        if generate_html:
             try:
                  html_file = output_path / "index.html"
                  html_content = html_generator.generate_html_report(
                       region_info=region_info,
                       output_base_dir=output_path,
                       enrichment_results=enrichment_results,
                       sample_analysis=sample_analysis,
                       gene_analysis=gene_analysis,
                       figures=figures,
                       saved_files=saved_files
                  )
                  with open(html_file, 'w', encoding='utf-8') as f: # specify encoding
                       f.write(html_content)
                  saved_files['html_report'] = str(html_file)
                  self.logger.info(f"HTML report generated: {html_file}")
             except ImportError as e:
                  self.logger.warning(f"Skipping HTML report generation due to missing required library (Jinja2): {e}")
             except Exception as e:
                  self.logger.error(f"Failed to generate or save HTML report ({html_file}): {e}", exc_info=True)

        # --- 9. Save metadata ---
        metadata_file = machine_dir / "metadata.json"
        if metadata is None:
            metadata = {}

        # Merge basic information and 'enhanced_stats'
        final_metadata = {
            'created_at': datetime.now().isoformat(),
            'region_info': region_info,
            'stats_summary': enhanced_stats,
            'files': saved_files,
            'analysis_parameters': metadata.get('analysis_parameters',{}),
            'input_data': metadata.get('input_data',{}),
            'execution_info': metadata.get('execution_info',{})
        }

        try:
             with open(metadata_file, 'w', encoding='utf-8') as f: # specify encoding
                  json.dump(final_metadata, f, indent=2, cls=NpEncoder)
             saved_files['metadata_json'] = str(metadata_file)
             self.logger.info(f"Metadata JSON saved: {metadata_file}")
        except Exception as e:
             self.logger.error(f"Failed to save metadata JSON ({metadata_file}): {e}", exc_info=True)

        self.logger.info(f"Analysis results package saving process completed: {output_dir}")
        return {
            'output_dir': str(output_path),
            'saved_files': saved_files,
            'metadata': final_metadata
        }
