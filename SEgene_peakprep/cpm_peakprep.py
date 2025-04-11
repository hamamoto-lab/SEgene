# === cpm_peakprep.py ===
import argparse
import logging
import os
import sys
import pandas as pd
import time

# --- cpm_peakprep_utils.py から関数をインポート ---
try:
    # Import only what we need (差分関連削除)
    from cpm_peakprep_utils import (
        get_sample_data,
        run_samtools_flagstat,
        run_featurecounts,
        calculate_logcpm,
        convert_bed_to_saf,
        convert_mergesv_to_saf,
        invert_dictionary
    )
except ImportError as e:
    print(f"Error importing functions from cpm_peakprep_utils.py: {e}", file=sys.stderr)
    print("Please ensure 'cpm_peakprep_utils.py' exists and necessary libraries are installed.", file=sys.stderr)
    sys.exit(1)

# --- 引数パーサーの設定 ---
def parse_arguments():
    """Parses command line arguments."""
    parser = argparse.ArgumentParser(
        description="Simplified ChIP-seq Quantification Pipeline: Runs flagstat, featureCounts, calculates logCPM.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # --- 必須引数 ---
    parser.add_argument("-b", "--bam_dir", required=True,
                        help="Directory containing BAM files.")
    parser.add_argument("-a", "--annotation_file", required=True,
                        help="Path to the annotation file (SAF, BED, or merge_SV format). Specify format with --is_mergesv_format if needed.")
    parser.add_argument("-o", "--output_dir", required=True,
                        help="Directory to save all output files.")

    # --- オプション引数: ファイルフォーマット ---
    parser.add_argument("--is_mergesv_format", action='store_true',
                        help="Specify if the annotation file is in 'merge_SV.tsv' format (chr_start_end in the first column).")

    # --- オプション引数: ファイル/サンプル選択 ---
    parser.add_argument("--filename_pattern", default="*.bam",
                        help="Wildcard pattern for BAM filenames.")
    parser.add_argument("--sample_delimiter", default=None,
                        help="Delimiter string to extract sample names from BAM filenames.")

    # --- オプション引数: ツールパスとスレッド数 ---
    parser.add_argument("--samtools_path", default="samtools", help="Path to samtools.")
    parser.add_argument("--featurecounts_path", default="featureCounts", help="Path to featureCounts.")
    parser.add_argument("-T", "--threads", type=int, default=4, help="Threads for samtools/featureCounts.")
    parser.add_argument("--single_end", action='store_true', help="Run in single-end mode (default: paired-end).")

    # --- オプション引数: featureCounts関連 ---
    parser.add_argument("--fc_basename", default=None,
                        help="Base name for featureCounts output/log files (default: ANNOTATION_fc).")
    parser.add_argument("--fc_options", nargs='+', default=None,
                        help="Additional options for featureCounts (e.g., --fc_options --minOverlap 10).")

    # --- オプション引数: logCPM計算関連 ---
    parser.add_argument("--add_peakid", action='store_true', help="Replace original IDs with sequential PeakIDs.")
    parser.add_argument("--id_column", default="Geneid", help="ID column name in SAF/featureCounts input.")
    parser.add_argument("--output_id_name", default="PeakID", help="ID column name in the final logCPM table.")
    parser.add_argument("--log_base", type=int, default=2, help="Log base (2 or 10). <=0 for no log.")
    parser.add_argument("--pseudocount", type=float, default=1.0, help="Pseudocount for log transform.")

    # --- オプション引数: 出力ファイル名 ---
    parser.add_argument("--logcpm_output_name", default="logCPM.tsv",
                        help="Filename for the logCPM table.")

    # --- オプション引数: ログ設定 ---
    parser.add_argument("--log_level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], help="Logging level.")
    parser.add_argument("--script_log_file", default="pipeline_run.log", help="Filename for script execution log. Set to '' to disable.")

    return parser.parse_args()

# --- メイン実行関数 ---
def main():
    """pipeline main function"""
    start_time = time.time()
    args = parse_arguments()

    # --- ロガー設定 (StreamHandlerを明示的に追加) ---
    log_level_str = args.log_level.upper()
    log_level = getattr(logging, log_level_str, logging.INFO)
    log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    date_format = '%Y-%m-%d %H:%M:%S'

    # ルートロガーを取得し、既存ハンドラをクリア
    root_logger = logging.getLogger()
    # 現在のハンドラをリストのコピーに対してループ処理して安全に削除
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)
        handler.close() # ハンドラを閉じる

    # ルートロガー自体のレベルを設定
    root_logger.setLevel(log_level)

    # 新しいコンソールハンドラ(StreamHandler)を追加し、レベルを設定
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(log_level) # ハンドラにもレベルを設定
    formatter = logging.Formatter(log_format, datefmt=date_format)
    ch.setFormatter(formatter)
    root_logger.addHandler(ch)

    # このスクリプト用の（子）ロガー取得
    logger = logging.getLogger("PipelineMain")

    # 出力ディレクトリ作成
    try:
        os.makedirs(args.output_dir, exist_ok=True)
        logger.info(f"Output directory: {args.output_dir}")
    except OSError as e:
        logger.critical(f"Failed to create output directory {args.output_dir}: {e}")
        sys.exit(1)

    # ファイルハンドラ設定
    if args.script_log_file:
        log_file_path = os.path.join(args.output_dir, args.script_log_file)
        try:
            fh = logging.FileHandler(log_file_path, mode='w')
            fh.setLevel(log_level)
            fh.setFormatter(formatter)
            root_logger.addHandler(fh)
            logger.info(f"Logging script execution to file: {log_file_path}")
        except Exception as e:
            # ファイル書き込みエラーが出ても処理は止めない
            logger.error(f"Failed to configure file logging to {log_file_path}: {e}")

    logger.debug(f"Parsed arguments: {vars(args)}")

    # --- 変数初期化 ---
    logcpm_df = None
    counts_file_path = None
    temp_saf_path = None # 一時SAFファイルのパス
    saf_file_to_use = None # featureCountsに実際に渡すSAFパス

    try: # メイン処理全体を try...finally で囲み、一時ファイルを確実に削除

        # --- アノテーションファイルの準備 ---
        logger.info("--- Preparing annotation file ---")
        annotation_file = args.annotation_file
        if not os.path.exists(annotation_file):
             logger.critical(f"Annotation file not found: {annotation_file}")
             sys.exit(1)

        # ファイル形式を判定し、必要ならSAFに変換
        if args.is_mergesv_format:
            logger.info("Input annotation specified as merge_SV format. Converting to temporary SAF...")
            temp_saf_path = convert_mergesv_to_saf(annotation_file, logger)
            if temp_saf_path is None: sys.exit(1)
            saf_file_to_use = temp_saf_path
        else:
            annot_lower = annotation_file.lower()
            if annot_lower.endswith(('.bed', '.bed.gz')):
                logger.info("Input annotation appears to be a BED file. Converting to temporary SAF...")
                temp_saf_path = convert_bed_to_saf(annotation_file, logger)
                if temp_saf_path is None: sys.exit(1)
                saf_file_to_use = temp_saf_path
            elif annot_lower.endswith('.saf'):
                 logger.info("Input annotation appears to be a SAF file. Using it directly.")
                 saf_file_to_use = annotation_file
            else:
                 logger.warning(f"Unknown extension for annotation file: {annotation_file}. Assuming it is in SAF format.")
                 saf_file_to_use = annotation_file

        # featureCounts出力のベース名を決定
        annot_base = os.path.splitext(os.path.basename(annotation_file))[0]
        logger.info(f"Using annotation for featureCounts: {saf_file_to_use}")
        
        # =============================
        # === BAMファイル処理 ===
        # =============================
        logger.info("--- Processing BAM files ---")
        try:
            # 1. Get sample data
            logger.info("Step 1: Getting sample info for BAM files...")
            sample_dict, bam_list = get_sample_data(
                args.bam_dir, logger, args.sample_delimiter, args.filename_pattern
            )
            if not bam_list: raise ValueError("No BAM files found.")
            bam_to_sample_map = invert_dictionary(sample_dict)
            logger.info(f"Found {len(sample_dict)} samples.")
            
            # 検出されたファイルリストをログ出力
            if bam_list:
                logger.debug(f"BAM files included ({len(bam_list)}):")
                for f_path in bam_list: logger.debug(f"  - {os.path.basename(f_path)}")
            
            # 検出されたファイルリストをファイル保存
            bam_list_path = os.path.join(args.output_dir, "bam_list.txt")
            try:
                with open(bam_list_path, 'w', encoding='utf-8') as f_out:
                    for bam_path in bam_list: f_out.write(bam_path + '\n')
                logger.info(f"Saved list of processed BAM files to: {bam_list_path}")
            except IOError as e:
                logger.warning(f"Could not write BAM list file {bam_list_path}: {e}")

            # 2. Run samtools flagstat
            logger.info("Step 2: Getting total mapped reads for BAM files...")
            flagstat_dir = os.path.join(args.output_dir, "flagstat")
            total_reads_map = run_samtools_flagstat(
                bam_list, flagstat_dir, logger, args.samtools_path, args.threads
            )
            logger.info(f"Extracted reads for {len(total_reads_map)} files.")
            if len(total_reads_map) != len(bam_list): 
                logger.warning("Could not get reads for all BAMs.")

            # 3. Run featureCounts
            logger.info("Step 3: Running featureCounts for BAM files...")
            fc_base = args.fc_basename if args.fc_basename else f"{annot_base}_fc"
            fc_output = os.path.join(args.output_dir, f"{fc_base}_featureCounts.txt")
            fc_log = os.path.join(args.output_dir, f"{fc_base}_featureCounts.log")
            counts_file_path = run_featurecounts(
                bam_list, saf_file_to_use, logger, fc_output, fc_log,
                args.threads, (not args.single_end), args.featurecounts_path, args.fc_options
            )
            if not counts_file_path: raise RuntimeError("featureCounts failed.")
            logger.info(f"featureCounts output: {counts_file_path}")

            # 4. Calculate logCPM
            logger.info("Step 4: Calculating logCPM...")
            try:
                counts_df = pd.read_csv(counts_file_path, sep='\t', comment='#', header=0)
                if counts_df.columns[0].startswith('#'): 
                    counts_df = pd.read_csv(counts_file_path, sep='\t', comment='#', header=1)
            except Exception as e: 
                raise RuntimeError(f"Failed to read counts file {counts_file_path}") from e

            log_base_val = args.log_base if args.log_base in [2, 10] else None
            logcpm_df = calculate_logcpm(
                counts_df=counts_df, 
                total_reads_dict=total_reads_map, 
                logger=logger,
                bam_path_to_sample_name=bam_to_sample_map, 
                add_peakid_column=args.add_peakid,
                id_column_name=args.id_column, 
                output_id_column_name=args.output_id_name,
                log_transform_base=log_base_val, 
                pseudocount=args.pseudocount
            )
            if logcpm_df is None: raise RuntimeError("logCPM calculation failed.")
            logger.info("logCPM calculation completed.")

            # 5. Save logCPM results
            output_path = os.path.join(args.output_dir, args.logcpm_output_name)
            logcpm_df.to_csv(output_path, sep='\t', index=False, float_format='%.4f')
            logger.info(f"Saved logCPM table to: {output_path}")

        except Exception as e:
            # 処理中のエラーは致命的として終了
            logger.critical(f"Critical error during data processing: {e}", exc_info=True)
            sys.exit(1)

        # --- パイプライン終了 ---
        end_time = time.time()
        logger.info(f"Pipeline finished successfully in {end_time - start_time:.2f} seconds.")

    finally:
        # 一時ファイルを確実に削除 - ログレベル変更：INFO→DEBUG
        if temp_saf_path and os.path.exists(temp_saf_path):
            try:
                os.remove(temp_saf_path)
                logger.debug(f"Removed temporary annotation file: {temp_saf_path}")
            except Exception as e:
                logger.warning(f"Failed to remove temporary annotation file {temp_saf_path}: {e}")

# --- スクリプト実行のエントリポイント ---
if __name__ == "__main__":
    main()