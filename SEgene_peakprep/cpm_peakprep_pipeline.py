import argparse
import logging
import os
import sys
import pandas as pd 

try:
    from cpm_peakprep_utils import (
        get_sample_data,
        run_samtools_flagstat,
        run_featurecounts,
        calculate_logcpm
    )
    # natsort, numpy, pandas は pipeline_utils.py 内で import されている想定
except ImportError as e:
    print(f"Error importing functions from pipeline_utils.py: {e}", file=sys.stderr)
    print("Please ensure 'pipeline_utils.py' exists and necessary libraries (natsort, pandas, numpy) are installed.", file=sys.stderr)
    sys.exit(1)
# --- インポートここまで ---

# --- 引数パーサーの設定 ---
def parse_arguments():
    """コマンドライン引数を解析する関数"""
    parser = argparse.ArgumentParser(
        description="ChIP-seq Quantification Pipeline: Extracts reads, runs featureCounts, calculates logCPM.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter # デフォルト値もヘルプに表示
    )

    # --- 必須引数 ---
    parser.add_argument("-b", "--bam_dir", required=True,
                        help="Directory containing input BAM files.")
    parser.add_argument("-a", "--saf_file", required=True,
                        help="Path to the annotation file in SAF format.")
    parser.add_argument("-o", "--output_dir", required=True,
                        help="Directory to save all output files.")

    # --- オプション引数: ファイル/サンプル選択 ---
    parser.add_argument("--filename_pattern", default="*.bam",
                        help="Wildcard pattern to filter BAM filenames.")
    parser.add_argument("--sample_delimiter", default=None,
                        help="Delimiter string to extract sample names from BAM filenames (part before delimiter is used).")

    # --- オプション引数: ツールパスとスレッド数 ---
    parser.add_argument("--samtools_path", default="samtools",
                        help="Path to samtools executable.")
    parser.add_argument("--featurecounts_path", default="featureCounts",
                        help="Path to featureCounts executable.")
    parser.add_argument("-T", "--threads", type=int, default=4,
                        help="Number of threads for samtools and featureCounts.")
    parser.add_argument("--single_end", action='store_true', # 指定するとTrueになるフラグ
                        help="Specify if reads are single-end (default is paired-end).")

    # --- オプション引数: featureCounts と logCPM計算 ---
    parser.add_argument("--fc_output_name", default=None,
                        help="Basename for featureCounts output file (default: derived from SAF filename). Suffix '_featureCounts.txt' is added.")
    parser.add_argument("--fc_log_name", default=None,
                        help="Basename for featureCounts log file (default: derived from SAF filename). Suffix '_featureCounts.log' is added.")
    parser.add_argument("--fc_options", nargs='+', default=None, # '+' は1つ以上の引数をリストとして受け取る
                        help="Additional options for featureCounts, separated by spaces (e.g., --fc_options -g gene_id --minOverlap 10).")
    parser.add_argument("--add_peakid", action='store_true',
                        help="Replace original GeneIDs with sequential PeakIDs (PeakID_XXXX).")
    parser.add_argument("--id_column", default="Geneid",
                        help="Name of the ID column in featureCounts output (input for logCPM).")
    parser.add_argument("--output_id_name", default="PeakID",
                        help="Name of the ID column in the final logCPM table.")
    parser.add_argument("--log_base", type=int, default=2,
                        help="Base for log transformation (e.g., 2 or 10). Set to 0 or negative for no log transform.")
    parser.add_argument("--pseudocount", type=float, default=1.0,
                        help="Pseudocount for log transformation (log(CPM + pseudocount)).")
    parser.add_argument("--logcpm_output_name", default="logCPM_table.tsv",
                        help="Filename for the final logCPM table.")

    # --- オプション引数: ログ設定 ---
    parser.add_argument("--log_level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help="Set logging level.")
    parser.add_argument("--script_log_file", default="pipeline_run.log",
                        help="Filename for the script execution log within the output directory (default: pipeline_run.log). Set to empty ('') to disable file logging.")

    return parser.parse_args()

# --- メイン実行関数 ---
def main():
    """パイプラインのメイン処理を実行する関数"""
    args = parse_arguments()

    # --- ロガー設定 ---
    log_level = getattr(logging, args.log_level.upper())
    log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    date_format = '%Y-%m-%d %H:%M:%S'

    # 基本的なコンソール出力設定
    logging.basicConfig(level=log_level, format=log_format, datefmt=date_format, stream=sys.stdout)
    # このスクリプト用のロガー取得
    logger = logging.getLogger("PipelineMain")

    # 出力ディレクトリ作成 (存在しない場合)
    try:
        os.makedirs(args.output_dir, exist_ok=True)
        logger.info(f"Output directory set to: {args.output_dir}")
    except OSError as e:
        logger.critical(f"Failed to create output directory {args.output_dir}: {e}")
        sys.exit(1)

    # オプションでファイルにもログ出力
    if args.script_log_file: # 空文字列でない場合
        log_file_path = os.path.join(args.output_dir, args.script_log_file)
        try:
            file_handler = logging.FileHandler(log_file_path, mode='w') # 上書きモード
            file_handler.setLevel(log_level)
            file_handler.setFormatter(logging.Formatter(log_format, datefmt=date_format))
            # Root logger にハンドラを追加 (他のモジュールからのログも拾うため)
            logging.getLogger().addHandler(file_handler)
            logger.info(f"Logging script execution to file: {log_file_path}")
        except Exception as e:
             logger.error(f"Failed to configure file logging to {log_file_path}: {e}")
    # --- ロガー設定ここまで ---

    logger.info("Starting ChIP-seq quantification pipeline.")
    # 引数をログに出力 (デバッグ用)
    logger.debug(f"Parsed arguments: {vars(args)}")

    # === ステップ 1: BAMファイルの情報を取得 ===
    logger.info("--- Step 1: Getting sample data and BAM file list ---")
    sample_dict = {}
    bam_list = []
    try:
        sample_dict, bam_list = get_sample_data(
            bam_folder=args.bam_dir,
            logger=logger, # このスクリプトのロガーを渡す
            sample_name_delimiter=args.sample_delimiter,
            filename_pattern=args.filename_pattern
        )
        if not bam_list:
            logger.critical("No matching BAM files found based on pattern and directory. Cannot proceed.")
            sys.exit(1) # BAMファイルがなければ続行不可
        logger.info(f"Found {len(sample_dict)} samples and {len(bam_list)} matching BAM files.")
    except Exception as e:
        logger.critical(f"Critical error during get_sample_data: {e}", exc_info=True)
        sys.exit(1) # ここでのエラーは致命的

    # === ステップ 2: samtools flagstat を実行 ===
    logger.info("--- Step 2: Running samtools flagstat ---")
    # flagstat結果の出力ディレクトリを作成
    flagstat_dir = os.path.join(args.output_dir, "flagstat_output")
    total_reads_map = {} # 結果を格納する辞書
    try:
        total_reads_map = run_samtools_flagstat(
            bam_files=bam_list,
            output_dir=flagstat_dir, # flagstat専用ディレクトリ
            logger=logger,
            samtools_path=args.samtools_path,
            threads=args.threads
        )
        if not total_reads_map:
            # リード数情報が全く得られなかった場合、警告を出す
            logger.warning("Could not extract mapped reads information from any BAM file via flagstat.")
            # CPM計算はできないが、featureCountsは実行できるかもしれないので続行は可能
        else:
            logger.info(f"Successfully extracted mapped reads info for {len(total_reads_map)} BAM files.")
    except Exception as e:
        # flagstatのエラーは致命的ではないかもしれないが、ログには残す
        logger.error(f"Error occurred during run_samtools_flagstat: {e}", exc_info=True)
        # ここで終了するかどうかは要件による (今回は続行)
        logger.warning("Proceeding without total reads information for potentially all samples.")

    # === ステップ 3: featureCounts を実行 ===
    logger.info("--- Step 3: Running featureCounts ---")
    # 出力ファイル名とログファイル名を決定 (デフォルト名生成ロジック)
    saf_base = os.path.splitext(os.path.basename(args.saf_file))[0]
    fc_output_name = args.fc_output_name if args.fc_output_name else f"{saf_base}_featureCounts.txt"
    fc_log_name = args.fc_log_name if args.fc_log_name else f"{saf_base}_featureCounts.log"
    fc_output_path = os.path.join(args.output_dir, fc_output_name)
    fc_log_path = os.path.join(args.output_dir, fc_log_name)
    # ペアエンドフラグを設定
    is_paired = not args.single_end

    counts_file_path = None # 結果パスを初期化
    try:
        counts_file_path = run_featurecounts(
            bam_files=bam_list,
            saf_file=args.saf_file,
            logger=logger,
            output_file=fc_output_path,
            log_file=fc_log_path,
            threads=args.threads,
            is_paired_end=is_paired,
            featurecounts_path=args.featurecounts_path,
            additional_options=args.fc_options
        )
        if counts_file_path is None:
            # featureCountsが失敗したらlogCPMは計算できないので終了
            logger.critical("featureCounts step failed. Cannot proceed to logCPM calculation.")
            sys.exit(1)
        logger.info(f"featureCounts completed. Output counts file: {counts_file_path}")
    except Exception as e:
        # featureCounts実行中の予期せぬエラー
        logger.critical(f"Critical error during run_featurecounts: {e}", exc_info=True)
        sys.exit(1)

    # === ステップ 4: LogCPM 計算 ===
    logger.info("--- Step 4: Calculating logCPM ---")
    # BAMパス -> サンプル名のマッピングを作成 (get_sample_dataの結果から)
    bam_to_sample_map = {v: k for k, v in sample_dict.items()} if sample_dict else {}
    if not bam_to_sample_map and total_reads_map:
        logger.warning("Could not create mapping from BAM path to sample name. Using BAM filenames as sample names in output.")

    # log_base=0 や負数を None に変換 (対数変換なしを示す)
    log_base = args.log_base if args.log_base in [2, 10] else None

    logcpm_df = None # 結果DataFrameを初期化
    try:
        # total_reads_map が空でも関数は呼べるが、結果はNaNになる
        if not total_reads_map:
             logger.warning("Total reads information is missing, logCPM values will likely be NaN.")

        logcpm_df = calculate_logcpm(
            counts_file=counts_file_path,
            total_reads_dict=total_reads_map if total_reads_map else {}, # 空辞書を渡す
            logger=logger,
            bam_path_to_sample_name=bam_to_sample_map, # マッピングを渡す
            add_peakid_column=args.add_peakid,
            id_column_name=args.id_column,
            output_id_column_name=args.output_id_name,
            log_transform_base=log_base,
            pseudocount=args.pseudocount
        )
        if logcpm_df is None:
            # logCPM計算が失敗した場合
            logger.error("LogCPM calculation step failed.")
            # 続行は可能かもしれないが、通常は結果が必要なので終了
            sys.exit(1)
        logger.info("LogCPM calculation completed.")
        logger.debug(f"LogCPM DataFrame head:\n{logcpm_df.head()}")

    except Exception as e:
        # logCPM計算中の予期せぬエラー
        logger.error(f"Error during calculate_logcpm: {e}", exc_info=True)
        sys.exit(1)

    # === ステップ 5: 結果の保存 ===
    logger.info("--- Step 5: Saving final logCPM table ---")
    # 最終的な出力ファイルパスを決定
    final_output_path = os.path.join(args.output_dir, args.logcpm_output_name)
    try:
        # DataFrame を TSV ファイルとして保存
        logcpm_df.to_csv(final_output_path, sep='\t', index=False, float_format='%.4f')
        logger.info(f"Final logCPM table saved successfully to: {final_output_path}")
    except Exception as e:
        # ファイル保存失敗
        logger.error(f"Failed to save the final logCPM table to {final_output_path}: {e}", exc_info=True)
        sys.exit(1) # 保存失敗は致命的

    # --- パイプライン終了 ---
    logger.info("ChIP-seq quantification pipeline finished successfully.")

# --- スクリプトが直接実行された場合のみ main() を呼び出す ---
if __name__ == "__main__":
    main()