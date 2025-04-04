#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import subprocess
import os
import re
import json
import logging
import traceback
import argparse
from pathlib import Path
from datetime import datetime
import sys
import pprint  

PREPROCESSOR_MODULE = "sedb_preprocessor.preprocess"
EXTRACTOR_MODULE = "sedb_extractor.main"
ANALYZER_MODULE = "sedb_analyzer.main"

# --- Logger Configuration ---
def setup_batch_logger(log_level=logging.INFO, log_file=None):
    """Configure logger for batch script"""
    logger = logging.getLogger("batch_processor")
    logger.setLevel(log_level)

    if logger.hasHandlers():
        logger.handlers.clear()

    log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    formatter = logging.Formatter(log_format)

    ch = logging.StreamHandler()
    ch.setLevel(log_level) 
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    if log_file:
        try:
            log_file_path = Path(log_file)
            log_file_path.parent.mkdir(parents=True, exist_ok=True) 
            fh = logging.FileHandler(log_file_path, mode='a', encoding='utf-8') 
            fh.setLevel(log_level) 
            fh.setFormatter(formatter)
            logger.addHandler(fh)
        except Exception as e:
            logger.error(f"Failed to create log file handler: {log_file} - {e}")

    return logger

# --- Helper Functions ---
def run_command(command_list, log_prefix=""):
    """Execute a command and return success/failure status and standard error"""
    logger = logging.getLogger("batch_processor") 
    logger.info(f"{log_prefix}Executing command: {' '.join(command_list)}")
    try:
        # Using text=True instead of deprecated universal_newlines=True (Python 3.7+)
        # Not redirecting stderr to stdout, handling them separately
        result = subprocess.run(
            command_list,
            capture_output=True,
            text=True,
            encoding='utf-8', # Specify output encoding
            check=True # Raise CalledProcessError if return code is non-zero
        )
        logger.info(f"{log_prefix}Command completed successfully.")
        # Output stdout/stderr at DEBUG level
        if result.stdout:
            logger.debug(f"{log_prefix}STDOUT:\n{result.stdout.strip()}")
        if result.stderr:
            logger.debug(f"{log_prefix}STDERR:\n{result.stderr.strip()}")
        return True, result.stderr # Return True and stderr on success
    except subprocess.CalledProcessError as e:
        logger.error(f"{log_prefix}Command failed (return code: {e.returncode})")
        # Display error content at ERROR level
        if e.stderr:
            logger.error(f"{log_prefix}STDERR:\n{e.stderr.strip()}")
        if e.stdout: # Error messages may also appear in stdout
            logger.error(f"{log_prefix}STDOUT:\n{e.stdout.strip()}")
        return False, e.stderr # Return False and stderr on failure
    except FileNotFoundError:
        err_msg = f"{log_prefix}Command not found: {command_list[0]}. Please verify your Python environment."
        logger.error(err_msg)
        return False, err_msg
    except Exception as e:
        err_msg = f"{log_prefix}Unexpected error during command execution: {e}"
        logger.error(err_msg)
        logger.error(traceback.format_exc())
        return False, str(e)

def parse_region_id(region_id):
    """Parse ID in chrX_START_END format"""
    match = re.match(r'(chr[\w]+)_(\d+)_(\d+)', str(region_id)) # Consider X, Y, M, convert ID to string
    if match:
        chrom, start, end = match.groups()
        try:
            # Convert start and end to integers
            return chrom, int(start), int(end)
        except ValueError:
            logger.warning(f"Cannot convert coordinates in region ID '{region_id}' to integers.")
            return None, None, None
    logger.warning(f"Region ID '{region_id}' format is not 'chrX_START_END'.")
    return None, None, None





def extract_results(extractor_output_dir, analyzer_output_dir):
    """
    Extracts key results from Extractor and Analyzer output directories.

    Args:
        extractor_output_dir (str or Path): Path to the extractor's output directory.
        analyzer_output_dir (str or Path): Path to the analyzer's output directory.

    Returns:
        dict: A dictionary containing extracted results:
              'se_unique_sample_count': Count of unique samples from SE data (from extractor summary).
              'linked_sample_count': Count of samples successfully linked to metadata (from analyzer summary).
              'significant_tissues': Comma-separated string of significant tissue names.
              'significant_tissues_fc': Comma-separated string of fold changes.
              'significant_tissues_pvalue': Comma-separated string of p-values.
              'significant_tissues_fdr': Comma-separated string of FDR values.
              'error_message': String containing error message if extraction failed, otherwise None.
    """
    logger = logging.getLogger("batch_processor") # Get batch logger
    results = {
        'se_unique_sample_count': 0, # Default 0
        'linked_sample_count': 0,    # Default 0
        'significant_tissues': '',
        'significant_tissues_fc': '',
        'significant_tissues_pvalue': '',
        'significant_tissues_fdr': '',
        'error_message': None
    }
    try:
        # --- Get unique sample count from SE data in Extractor's summary.json ---
        extractor_summary_path = Path(extractor_output_dir) / "summary.json"
        extractor_summary_data = None # Initialize
        if extractor_summary_path.exists():
             logger.debug(f"Reading extractor summary: {extractor_summary_path}")
             with open(extractor_summary_path, 'r', encoding='utf-8') as f:
                 extractor_summary_data = json.load(f)
            #  print(f"\nDEBUG_EXTRACT: Content of {extractor_summary_path}:", file=sys.stderr)
            #  pprint.pprint(extractor_summary_data, stream=sys.stderr, indent=2)
             try:

                  se_count = extractor_summary_data.get('unique_samples', None) 
                  if se_count is None:
                       se_count = extractor_summary_data.get('stats', {}).get('unique_samples', None)
                       if se_count is not None:
                            pass
                            # print(f"DEBUG_EXTRACT: Found se_unique_sample_count under 'stats' key (value={se_count}).", file=sys.stderr)
                       else:
                            # print(f"DEBUG_EXTRACT: Key 'unique_samples' not found at top level or under 'stats' in {extractor_summary_path}.", file=sys.stderr)
                            results['se_unique_sample_count'] = 0 # Set to 0 if not found

                  if se_count is not None: # Convert to int only if value is found
                      results['se_unique_sample_count'] = int(se_count) # Convert to int

                #   print(f"DEBUG_EXTRACT: Extracted se_unique_sample_count = {results['se_unique_sample_count']}", file=sys.stderr) # Verify value
             except (TypeError, ValueError) as e:
                #   print(f"DEBUG_EXTRACT: Error converting se_unique_sample_count (value was {se_count}): {e}", file=sys.stderr)
                  results['se_unique_sample_count'] = 0 # Set to 0 on error
        else:
             logger.warning(f"Extractor summary file not found: {extractor_summary_path}")
             # results['error_message'] = "Extractor Summary JSON not found" # Record error if needed

        # --- Get linked sample count from Analyzer's summary.json ---
        analyzer_summary_path = Path(analyzer_output_dir) / "summary.json"
        analyzer_summary_data = None # Initialize
        if analyzer_summary_path.exists():
            logger.debug(f"Reading analyzer summary: {analyzer_summary_path}")
            with open(analyzer_summary_path, 'r', encoding='utf-8') as f:
                analyzer_summary_data = json.load(f)
            # print(f"\nDEBUG_EXTRACT: Content of {analyzer_summary_path}:", file=sys.stderr)
            # pprint.pprint(analyzer_summary_data, stream=sys.stderr, indent=2)
            try:
                # Look for 'linked_sample_count' or 'sample_count' key
                linked_count = analyzer_summary_data.get('linked_sample_count', None) 
                if linked_count is None:
                    linked_count = analyzer_summary_data.get('sample_count', None) # Try legacy key name
                    if linked_count is not None:
                         pass
                        #  print(f"DEBUG_EXTRACT: Found linked_sample_count using fallback key 'sample_count'.", file=sys.stderr)

                if linked_count is None:
                    #  print(f"DEBUG_EXTRACT: Key 'linked_sample_count' or 'sample_count' not found in {analyzer_summary_path}.", file=sys.stderr)
                     results['linked_sample_count'] = 0
                else:
                     results['linked_sample_count'] = int(linked_count) # Convert to int
                # print(f"DEBUG_EXTRACT: Extracted linked_sample_count = {results['linked_sample_count']}", file=sys.stderr) # Verify value
            except (TypeError, ValueError) as e:
                # print(f"DEBUG_EXTRACT: Error converting linked_sample_count (value was {linked_count}): {e}", file=sys.stderr)
                results['linked_sample_count'] = 0 # Set to 0 on error
        else:
             logger.warning(f"Analyzer summary file not found: {analyzer_summary_path}")
             results['error_message'] = "Analyzer Summary JSON not found"

        # --- Extract enrichment results ---
        enrichment_file_path = None
        enrichment_parquet = Path(analyzer_output_dir) / "machine_readable" / "data" / "enrichment_results.parquet"
        enrichment_csv = Path(analyzer_output_dir) / "human_readable" / "tables" / "enrichment_results.csv"

        if enrichment_parquet.exists():
            enrichment_file_path = enrichment_parquet
        elif enrichment_csv.exists():
            enrichment_file_path = enrichment_csv

        df_enrich = pd.DataFrame() # Initialize
        if enrichment_file_path:
            try:
                if enrichment_file_path.suffix == '.parquet':
                    df_enrich = pd.read_parquet(enrichment_file_path)
                elif enrichment_file_path.suffix == '.csv':
                    df_enrich = pd.read_csv(enrichment_file_path)
                logger.info(f"Loaded enrichment results: {enrichment_file_path}")
            except Exception as e:
                logger.warning(f"Failed to read enrichment results file ({enrichment_file_path}): {e}")
                if not results['error_message']: results['error_message'] = "Enrichment file read error"
        else:
             logger.warning(f"Analyzer enrichment results file not found (Parquet/CSV).")
             # Don't overwrite error message if enrichment exists but summary.json doesn't
             if not results['error_message']: results['error_message'] = "Enrichment results not found"

        # Extract and format significant results
        if not df_enrich.empty and 'Significant' in df_enrich.columns:
            # Handle column name variations
            fdr_col = next((col for col in ['adjusted p-value', 'fdr'] if col in df_enrich.columns), None)
            fc_col = next((col for col in ['Fold Change', 'fold_change'] if col in df_enrich.columns), None)
            pval_col = next((col for col in ['p-value', 'pvalue'] if col in df_enrich.columns), None)
            tissue_col = next((col for col in ['Tissue Type', 'tissue_type', 'tissue'] if col in df_enrich.columns), None)

            if all([fdr_col, fc_col, pval_col, tissue_col]):
                # Consider possibility that Significant column is not boolean
                try:
                    # Convert to boolean type before comparison
                    sig_rows = df_enrich[df_enrich['Significant'].astype(bool) == True].sort_values(fdr_col)
                except Exception as filter_e:
                     logger.error(f"Error filtering significant rows in enrichment results: {filter_e}")
                     sig_rows = pd.DataFrame() # Empty on error

                if not sig_rows.empty:
                    # Extract and format (with error handling)
                    try:
                        results['significant_tissues'] = ','.join(sig_rows[tissue_col].astype(str))
                        results['significant_tissues_fc'] = ','.join(sig_rows[fc_col].map('{:.2f}'.format))
                        results['significant_tissues_pvalue'] = ','.join(sig_rows[pval_col].map('{:.2e}'.format))
                        results['significant_tissues_fdr'] = ','.join(sig_rows[fdr_col].map('{:.2e}'.format))
                    except KeyError as ke:
                         logger.error(f"Missing expected column during significant result formatting: {ke}")
                         if not results['error_message']: results['error_message'] = "Enrichment result formatting error"
                    except Exception as fmt_e:
                         logger.error(f"Error formatting significant results: {fmt_e}")
                         if not results['error_message']: results['error_message'] = "Enrichment result formatting error"
                else:
                    logger.info("No significant enrichment tissues found for this region.")
            else:
                missing_cols_names = [name for col, name in zip([fdr_col, fc_col, pval_col, tissue_col], ['FDR', 'FC', 'PValue', 'Tissue']) if not col]
                logger.warning(f"Enrichment results are missing required columns: {', '.join(missing_cols_names)}")
                if not results['error_message']: results['error_message'] = f"Enrichment column(s) missing ({', '.join(missing_cols_names)})"

        elif not df_enrich.empty:
             logger.warning(f"Enrichment results file ({enrichment_file_path}) is missing the 'Significant' column.")
             if not results['error_message']: results['error_message'] = "Enrichment format error ('Significant' column missing)"

    except Exception as e:
        logger.error(f"Unexpected error during result extraction ({extractor_output_dir}, {analyzer_output_dir}): {e}", exc_info=True)
        results['error_message'] = f"Result extraction failed unexpectedly: {e}"

    # print(f"DEBUG_EXTRACT: Final extracted counts for TSV: SE Unique={results['se_unique_sample_count']}, Linked={results['linked_sample_count']}", file=sys.stderr)
    return results





# --- Main Processing ---
def run_batch_analysis(input_file, input_format, bed_file, metadata_file, main_output_dir, max_rows=None, debug=False, batch_log_file=None):
    """Main function to execute batch analysis"""

    # Configure logger for the batch script itself
    logger = setup_batch_logger(logging.DEBUG if debug else logging.INFO, batch_log_file)

    logger.info("Starting batch analysis process...")
    logger.info(f"Input region file: {input_file} (format: {input_format})")
    logger.info(f"Input BED file: {bed_file}")
    logger.info(f"Input metadata file: {metadata_file}")
    logger.info(f"Main output directory: {main_output_dir}")
    if max_rows:
        logger.info(f"Maximum number of rows to process: {max_rows}")
    if debug:
        logger.info("Debug mode: Enabled")
        logger.info("  -> Detailed logs will be output and debug logs from each tool will be saved to files.")

    # --- Set up output directories ---
    main_output_path = Path(main_output_dir)
    preprocessor_output_dir = main_output_path / "00_preprocessed_data"
    batch_results_dir = main_output_path / "01_batch_region_results"
    final_tsv_dir = main_output_path

    try:
        preprocessor_output_dir.mkdir(parents=True, exist_ok=True)
        batch_results_dir.mkdir(parents=True, exist_ok=True)
        final_tsv_dir.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        logger.error(f"Failed to create output directories: {e}")
        return None

    # --- 1. Execute Preprocessor (once only) ---
    logger.info("--- Step 1: Running Preprocessor ---")
    cmd_preprocess = [
        "python", "-m", PREPROCESSOR_MODULE,
        "--bed-file", bed_file,
        "--metadata", metadata_file,
        "--output", str(preprocessor_output_dir),
        "-f" # Enable standard human chromosome filtering by default
    ]
    if debug:
        cmd_preprocess.append("-d")
        # Optionally save preprocessor log file
        # prep_log = preprocessor_output_dir / "preprocessor.log"
        # cmd_preprocess.extend(["--log-file", str(prep_log)])

    # success, _, stderr_prep = run_command(cmd_preprocess, "[Preprocessor] ")
    # if not success:
    #     logger.error("Preprocessor execution failed. Aborting batch process.")
    #     return None
    # logger.info("--- Preprocessor completed successfully ---")

    # Before correction
    # success, _, stderr_prep = run_command(cmd_preprocess, "[Preprocessor] ")
    # After correction
    success, stderr_prep = run_command(cmd_preprocess, "[Preprocessor] ")
    if not success:
        logger.error(f"Preprocessor execution failed. STDERR: {stderr_prep}") # Add error details to log
        logger.error("Aborting batch process.")
        return None


    # --- 2. Load region list ---
    logger.info("--- Step 2: Loading input region list ---")
    regions_df = None
    if input_format == 'tsv':
        try:
            regions_df = pd.read_csv(input_file, sep='\t')
            # Use first column as 'region_id' for TSV
            regions_df.rename(columns={regions_df.columns[0]: 'region_id'}, inplace=True)
            # Keep a copy of original columns
            input_columns = list(regions_df.columns)
            logger.info(f"Loaded {len(regions_df)} rows from TSV file.")
        except Exception as e:
            logger.error(f"Failed to load region TSV file ({input_file}): {e}")
            return None
    elif input_format == 'bed':
        try:
            # Load BED file (minimum 3 columns: chrom, start, end)
            # Assume no header (header=None)
            df_bed_regions = pd.read_csv(input_file, sep='\t', header=None, usecols=[0, 1, 2], names=['chrom', 'start', 'end'])
            # Generate 'region_id' column (chrX_START_END format)
            df_bed_regions['region_id'] = df_bed_regions.apply(lambda row: f"{row['chrom']}_{row['start']}_{row['end']}", axis=1)
            # Only need region_id
            regions_df = df_bed_regions[['region_id']].copy()
            input_columns = ['region_id'] # Only this column from BED
            logger.info(f"Loaded {len(regions_df)} regions from BED file.")
        except Exception as e:
            logger.error(f"Failed to load region BED file ({input_file}): {e}")
            return None
    else:
        logger.error(f"Unknown input format: {input_format}")
        return None

    if max_rows is not None and max_rows > 0:
        logger.info(f"Limiting processing to maximum {max_rows} rows.")
        regions_df = regions_df.head(max_rows)

    # --- Prepare columns for results ---
    results_columns = [
        'se_unique_sample_count', 'linked_sample_count',
        'significant_tissues', 'significant_tissues_fc',
        'significant_tissues_pvalue', 'significant_tissues_fdr',
        'processing_error'
    ]
    for col in results_columns:
        if col not in regions_df.columns: # Add only if not already existing
             regions_df[col] = None

    # --- 3. Process each region (loop) ---
    logger.info("--- Step 3: Running Extractor and Analyzer for each region ---")
    total_regions = len(regions_df)
    for idx in regions_df.index: # Use index to access with loc/at
        region_id = str(regions_df.loc[idx, 'region_id']) # Get from 'region_id' column
        logger.info(f"--- Processing: Region {idx+1}/{total_regions}: {region_id} ---")

        # Parse coordinates (expecting chrX_START_END format)
        chrom, start, end = parse_region_id(region_id)
        if not all([chrom, start, end]):
            logger.warning(f"Failed to parse region ID '{region_id}'. Skipping.")
            regions_df.loc[idx, 'processing_error'] = "Invalid Region ID Format"
            continue

        # Set up output directories for this region
        # Replace characters that cannot be used in file names/paths (e.g., ':')
        safe_region_id = region_id.replace(":", "_")
        region_results_dir = batch_results_dir / f"region_{idx+1}_{safe_region_id}"
        extractor_output_dir = region_results_dir / "extractor_output"
        analyzer_output_dir = region_results_dir / "analyzer_output"
        # Create directories before running extractor/analyzer (tools should also create with exist_ok=True, but just in case)
        extractor_output_dir.mkdir(parents=True, exist_ok=True)
        analyzer_output_dir.mkdir(parents=True, exist_ok=True)

        # --- 3a. Run Extractor ---
        region_str = f"{chrom}:{start}-{end}" # Format for Extractor
        cmd_extract = [
            "python", "-m", EXTRACTOR_MODULE,
            "--bed-file", bed_file, # Note: Use cleaned BED path here
            "--metadata", metadata_file,
            "--region", region_str,
            "--output", str(extractor_output_dir),
            "--include-distributions", # Always output distribution data
            # "--export-bed", # Control with batch script argument if needed
        ]
        if debug:
             cmd_extract.append("-d")
             # Extractor log file path
             ext_log = extractor_output_dir / f"extractor_{safe_region_id}.log"
             cmd_extract.extend(["--log-file", str(ext_log)])

        logger.info("Running Extractor...")
        # success_ext, _, stderr_ext = run_command(cmd_extract, f"[Extractor {region_id}] ")
        # if not success_ext:
        #     logger.warning(f"Extractor execution failed (region: {region_id}). Skipping analysis for this region.")
        #     regions_df.loc[idx, 'processing_error'] = f"Extractor Failed: {stderr_ext[:200]}" # Record part of error message
        #     continue

        # Before correction
        # success_ext, _, stderr_ext = run_command(cmd_extract, f"[Extractor {region_id}] ")
        # After correction
        success_ext, stderr_ext = run_command(cmd_extract, f"[Extractor {region_id}] ")
        if not success_ext:
            logger.warning(f"Extractor execution failed (region: {region_id}). Skipping analysis for this region.")
            # Record part of error message (modified to convert to string and limit length)
            regions_df.loc[idx, 'processing_error'] = f"Extractor Failed: {str(stderr_ext)[:200]}"
            continue


        # --- 3b. Run Analyzer ---
        cmd_analyze = [
            "python", "-m", ANALYZER_MODULE,
            "--input", str(extractor_output_dir),
            "--background", str(preprocessor_output_dir),
            "--output", str(analyzer_output_dir),
            # "--no-figures", # Control with batch script argument if needed
            # "--no-html-report", # Control with batch script argument if needed
        ]
        if debug:
             cmd_analyze.append("-d")
             # Analyzer log file path
             ana_log = analyzer_output_dir / f"analyzer_{safe_region_id}.log"
             cmd_analyze.extend(["--log-file", str(ana_log)])

        logger.info("Running Analyzer...")
        # success_ana, _, stderr_ana = run_command(cmd_analyze, f"[Analyzer {region_id}] ")
        # if not success_ana:
        #     logger.warning(f"Analyzer execution failed (region: {region_id}).")
        #     regions_df.loc[idx, 'processing_error'] = f"Analyzer Failed: {stderr_ana[:200]}"
        #     continue


        # Before correction
        # success_ana, _, stderr_ana = run_command(cmd_analyze, f"[Analyzer {region_id}] ")
        # After correction
        success_ana, stderr_ana = run_command(cmd_analyze, f"[Analyzer {region_id}] ")
        if not success_ana:
            logger.warning(f"Analyzer execution failed (region: {region_id}).")
            # Record part of error message (modified to convert to string and limit length)
            regions_df.loc[idx, 'processing_error'] = f"Analyzer Failed: {str(stderr_ana)[:200]}"
            continue



        # --- 3c. Extract and store results ---
        logger.info("Extracting analysis results...")
        analysis_summary = extract_results(extractor_output_dir, analyzer_output_dir)

        regions_df.loc[idx, 'se_unique_sample_count'] = analysis_summary['se_unique_sample_count']
        regions_df.loc[idx, 'linked_sample_count'] = analysis_summary['linked_sample_count']
        regions_df.loc[idx, 'significant_tissues'] = analysis_summary['significant_tissues']
        regions_df.loc[idx, 'significant_tissues_fc'] = analysis_summary['significant_tissues_fc']
        regions_df.loc[idx, 'significant_tissues_pvalue'] = analysis_summary['significant_tissues_pvalue']
        regions_df.loc[idx, 'significant_tissues_fdr'] = analysis_summary['significant_tissues_fdr']
        if analysis_summary['error_message']:
             # Append to existing error message rather than overwrite
             current_error = regions_df.loc[idx, 'processing_error']
             new_error = analysis_summary['error_message']
             if pd.isna(current_error) or current_error == "":
                  regions_df.loc[idx, 'processing_error'] = f"Result Extraction Error: {new_error}"
             else: # If there's already another error
                  regions_df.loc[idx, 'processing_error'] = f"{current_error}; Result Extraction Error: {new_error}"

        logger.info(f"--- Region {idx+1}/{total_regions}: {region_id} processing completed ---")

    # --- 4. Save final results TSV ---
    logger.info("--- Step 4: Saving final results to TSV file ---")
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    # Generate output filename based on input filename (without extension)
    input_basename = Path(input_file).stem
    output_tsv_filename = f"{input_basename}_batch_results_{timestamp}.tsv"
    output_tsv_path = final_tsv_dir / output_tsv_filename
    try:
        # Save with original column order + new result columns (use input_columns)
        final_columns = input_columns + [col for col in results_columns if col not in input_columns]
        regions_df.to_csv(output_tsv_path, sep='\t', index=False, columns=final_columns)
        logger.info(f"Batch analysis completed. Results saved to: {output_tsv_path}")
        return str(output_tsv_path)
    except Exception as e:
        logger.error(f"Failed to save final results TSV ({output_tsv_path}): {e}")
        return None

# --- Argument handling when executed as script ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Batch analyze multiple genomic regions using SEdb Extractor and Analyzer, and summarize results in a TSV.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--input_file", required=True,
                        help="Input file defining regions to process (TSV or BED format).")
    parser.add_argument("--input_format", choices=['tsv', 'bed'], default='tsv',
                        help="Input file format ('tsv' or 'bed'). For TSV, the first column must be 'chrX_START_END'; for BED, the first 3 columns must be 'chr', 'start', 'end'.")
    parser.add_argument("--bed_file", required=True,
                        help="Path to original SEdb BED file for analysis (e.g., data/SEdb/SE_package_hg38_cleaned.bed).")
    parser.add_argument("--metadata_file", required=True,
                        help="Path to original metadata file for analysis (e.g., data/SEdb/human_sample_information_sedb2.txt).")
    parser.add_argument("--output_dir", required=True,
                        help="Main directory to save all outputs (intermediate and final results).")
    parser.add_argument("--max_rows", type=int, default=None,
                        help="Maximum number of rows to process from input file (optional, for testing).")
    parser.add_argument("--debug", action="store_true",
                        help="Enable debug mode. Detailed logs will be output and debug logs from each tool will be saved to files.")
    parser.add_argument("--log_file", default=None,
                        help="Path to save log file for this batch script itself (optional).")

    args = parser.parse_args()

    # Run batch processing
    final_tsv = run_batch_analysis(
        input_file=args.input_file,
        input_format=args.input_format,
        bed_file=args.bed_file,
        metadata_file=args.metadata_file,
        main_output_dir=args.output_dir,
        max_rows=args.max_rows,
        debug=args.debug,
        batch_log_file=args.log_file
    )

    if final_tsv:
        logging.info("Batch script completed successfully.")
        sys.exit(0)
    else:
        logging.error("Batch script exited with errors.")
        sys.exit(1)