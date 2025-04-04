"""
sedb_analyzer.main - Main execution script for SEdb statistical analyzer

This module executes the statistical analysis and enrichment testing
on extracted region data from SEdb.
"""

import os
import sys
import logging
import time
import traceback
from pathlib import Path
import pandas as pd
import numpy as np
from typing import Dict, Any, Optional, Union, List, Tuple

from .cli import parse_args, validate_args, setup_analyzer_logger
from .loader import AnalysisDataLoader
from .statistics import StatisticsAnalyzer
from .packager import AnalysisPackager

def main():
    """Main execution function for the SEdb analyzer"""
    # Start timing
    start_time = time.time()
    
    # 1. Parse command line arguments
    args = parse_args()
    
    # 2. Set up logger
    log_level = logging.DEBUG if args.debug else logging.INFO
    logger = setup_analyzer_logger(log_level, args.log_file)
    
    logger.info("Starting SEdb statistical analysis")
    
    # 3. Validate arguments
    if not validate_args(args, logger):
        logger.error("Argument validation failed")
        return 1
    
    try:
        # 4. Load data packages
        loader = AnalysisDataLoader(logger=logger)
        
        # 4.1 Load extracted region data
        logger.info(f"Loading extracted region data from: {args.input}")
        extraction_data = loader.load_extraction_package(args.input)
        
        # 4.2 Load background data
        logger.info(f"Loading background data from: {args.background}")
        background_data = loader.load_preprocessed_package(args.background)
        
        # 5. Prepare data for analysis
        
        # 5.1 Get region information
        region = extraction_data.get('region', {})
        region_name = region.get('name', f"{region.get('chrom', 'unknown')}:{region.get('start', 0)}-{region.get('end', 0)}")
        logger.info(f"Analyzing region: {region_name}")
        
        # 5.2 Get sample and SE data
        overlapping_se = extraction_data.get('overlapping_se', pd.DataFrame())
        sample_info = extraction_data.get('sample_info', pd.DataFrame())
        
        # 5.3 Get tissue distribution data
        region_tissue_counts = None
        if 'distributions' in extraction_data and 'tissue_counts' in extraction_data['distributions']:
            tissue_data = extraction_data['distributions']['tissue_counts']
            
            # Convert to Series if needed
            if isinstance(tissue_data, pd.DataFrame):
                logger.debug(f"Converting tissue data DataFrame to Series, columns: {list(tissue_data.columns)}")
                # Try different column naming patterns
                if 'tissue_type' in tissue_data.columns and 'count' in tissue_data.columns:
                    region_tissue_counts = pd.Series(tissue_data['count'].values, index=tissue_data['tissue_type'])
                elif 'Tissue Type' in tissue_data.columns and 'Count' in tissue_data.columns:
                    region_tissue_counts = pd.Series(tissue_data['Count'].values, index=tissue_data['Tissue Type'])
                elif 'term' in tissue_data.columns and 'count' in tissue_data.columns:
                    region_tissue_counts = pd.Series(tissue_data['count'].values, index=tissue_data['term'])
                else:
                    logger.warning(f"Unrecognized tissue data format with columns: {list(tissue_data.columns)}")
            elif isinstance(tissue_data, pd.Series):
                region_tissue_counts = tissue_data
                logger.debug("Tissue data is already in Series format")
            else:
                logger.warning(f"Unrecognized tissue data type: {type(tissue_data)}")
        else:
            logger.warning("No tissue distribution data found in extraction package")
        
        # 5.4 Get background tissue distribution
        background_tissue_counts = background_data.get('tissue_distribution')
        if background_tissue_counts is None or len(background_tissue_counts) == 0:
            logger.warning("No background tissue distribution data found")
        else:
            logger.info(f"Found background tissue distribution with {len(background_tissue_counts)} tissue types")
        
        # 5.5 Get gene counts data
        gene_counts = None
        if 'distributions' in extraction_data and 'gene_counts' in extraction_data['distributions']:
            gene_data = extraction_data['distributions']['gene_counts']
            
            # Convert to Series if needed
            if isinstance(gene_data, pd.DataFrame):
                logger.debug(f"Converting gene data DataFrame to Series, columns: {list(gene_data.columns)}")
                # Try different column naming patterns
                if 'gene' in gene_data.columns and 'count' in gene_data.columns:
                    gene_counts = pd.Series(gene_data['count'].values, index=gene_data['gene'])
                elif 'Gene' in gene_data.columns and 'Count' in gene_data.columns:
                    gene_counts = pd.Series(gene_data['Count'].values, index=gene_data['Gene'])
                else:
                    logger.warning(f"Unrecognized gene data format with columns: {list(gene_data.columns)}")
            elif isinstance(gene_data, pd.Series):
                gene_counts = gene_data
                logger.debug("Gene data is already in Series format")
            else:
                logger.warning(f"Unrecognized gene data type: {type(gene_data)}")
        else:
            logger.info("No gene distribution data found in extraction package")
        
        # 6. Initialize statistics analyzer
        stat_analyzer = StatisticsAnalyzer(logger=logger)
        
        # 7. Perform statistical analyses
        
        # 7.1 Sample distribution analysis
        sample_analysis = None
        if sample_info is not None and len(sample_info) > 0:
            logger.info(f"Analyzing sample distribution for {len(sample_info)} samples")
            sample_analysis = stat_analyzer.analyze_samples_distribution(sample_info)
            logger.info(f"Sample analysis completed: {sample_analysis.get('sample_count', 0)} samples")
        else:
            logger.warning("No sample information available for analysis")
        
        # 7.2 Enrichment analysis
        enrichment_results = None
        if region_tissue_counts is not None and background_tissue_counts is not None:
            logger.info("Performing tissue enrichment analysis")
            enrichment_results = stat_analyzer.test_tissue_enrichment(
                region_tissue_counts,
                background_tissue_counts,
                method=args.enrichment_method,
                correction=args.correction_method,
                min_region_count=args.min_region_count,
                pvalue_threshold=args.pvalue_threshold,
                fc_threshold=args.fc_threshold
            )
            
            # Count significant results
            if not enrichment_results.empty:
                sig_count = enrichment_results['Significant'].sum() if 'Significant' in enrichment_results.columns else 0
                depleted_count = enrichment_results['Depleted'].sum() if 'Depleted' in enrichment_results.columns else 0
                logger.info(f"Enrichment analysis completed: {len(enrichment_results)} tissues, "
                           f"significant: {sig_count}, depleted: {depleted_count}")
            else:
                logger.warning("Enrichment analysis yielded no results")
        else:
            logger.warning("Missing tissue counts data for enrichment analysis")
        
        # 7.3 Gene distribution analysis
        gene_analysis = None
        if gene_counts is not None and len(gene_counts) > 0:
            logger.info(f"Analyzing gene distribution for {len(gene_counts)} genes")
            gene_analysis = stat_analyzer.analyze_gene_distribution(gene_counts)
            logger.info(f"Gene analysis completed: {gene_analysis.get('gene_count', 0)} genes")
        else:
            logger.info("No gene count data available for analysis")
        
        # 8. Generate figures if requested
        figures = {}
        if args.figures:
            logger.info("Generating figures")
            
            # 8.1 Sample distribution figure
            if sample_analysis:
                logger.debug("Generating sample distribution figure")
                sample_fig = stat_analyzer.plot_samples_distribution(sample_analysis)
                if sample_fig:
                    figures['sample_distribution'] = sample_fig
                    logger.debug("Sample distribution figure generated")
                else:
                    logger.warning("Failed to generate sample distribution figure")
            
            # 8.2 Enrichment analysis figure
            if enrichment_results is not None and not enrichment_results.empty:
                logger.debug("Generating tissue enrichment figure")
                enrichment_fig = stat_analyzer.plot_tissue_enrichment(enrichment_results)
                if enrichment_fig:
                    figures['tissue_enrichment'] = enrichment_fig
                    logger.debug("Tissue enrichment figure generated")
                else:
                    logger.warning("Failed to generate tissue enrichment figure")
            
            # 8.3 Gene distribution figure
            if gene_analysis:
                logger.debug("Generating gene distribution figure")
                gene_fig = stat_analyzer.plot_gene_distribution(gene_analysis)
                if gene_fig:
                    figures['gene_distribution'] = gene_fig
                    logger.debug("Gene distribution figure generated")
                else:
                    logger.warning("Failed to generate gene distribution figure")
            
            logger.info(f"Generated {len(figures)} figures")
        else:
            logger.info("Figure generation disabled")
        
        # 9. Package and save results
        logger.info(f"Packaging analysis results to: {args.output}")
        
        # 9.1 Prepare metadata
        metadata = {
            'analysis_parameters': {
                'enrichment_method': args.enrichment_method,
                'correction_method': args.correction_method,
                'pvalue_threshold': args.pvalue_threshold,
                'fc_threshold': args.fc_threshold,
                'min_region_count': args.min_region_count
            },
            'input_data': {
                'extraction_dir': args.input,
                'background_dir': args.background
            },
            'execution_info': {
                'date': time.strftime('%Y-%m-%d %H:%M:%S'),
                'figures_generated': args.figures,
                'html_report_generated': args.html_report
            }
        }
        
        # 9.2 Initialize packager and save results
        packager = AnalysisPackager(logger=logger)
        
        package_info = packager.save_analysis_package(
            output_dir=args.output,
            region_info=region,
            enrichment_results=enrichment_results,
            sample_analysis=sample_analysis,
            gene_analysis=gene_analysis,
            figures=figures if args.figures else None,
            metadata=metadata,
            generate_html=args.html_report,
            image_formats=args.image_formats,
            data_formats_machine=args.data_formats_machine,
            data_formats_human=args.data_formats_human
        )
        
        # 10. Log summary and completion
        elapsed_time = time.time() - start_time
        logger.info(f"SEdb analysis completed in {elapsed_time:.2f} seconds")
        logger.info(f"Results saved to: {package_info['output_dir']}")
        
        # Log important statistics summary
        if sample_analysis:
            logger.info(f"Analyzed samples: {sample_analysis.get('sample_count', 0)}")
        
        if enrichment_results is not None and not enrichment_results.empty and 'Significant' in enrichment_results.columns:
            sig_results = enrichment_results[enrichment_results['Significant']]
            if not sig_results.empty:
                logger.info(f"Significantly enriched tissues: {len(sig_results)}")
                # Show top 3 results
                top_tissues = sig_results.sort_values('adjusted p-value').head(3)
                for _, row in top_tissues.iterrows():
                    logger.info(f"  - {row['Tissue Type']}: "
                               f"FC={row.get('Fold Change', 0):.2f}, "
                               f"p-adj={row.get('adjusted p-value', 1):.1e}")
        
        return 0
        
    except FileNotFoundError as e:
        logger.error(f"File not found: {e}")
        if args.debug:
            logger.debug(traceback.format_exc())
        return 1
    except ValueError as e:
        logger.error(f"Value error: {e}")
        if args.debug:
            logger.debug(traceback.format_exc())
        return 1
    except ImportError as e:
        logger.error(f"Missing dependency: {e}")
        if args.debug:
            logger.debug(traceback.format_exc())
        return 1
    except Exception as e:
        logger.error(f"Unexpected error during analysis: {e}")
        if args.debug:
            logger.debug(traceback.format_exc())
        return 1

if __name__ == "__main__":
    sys.exit(main())