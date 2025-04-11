"""
Utility functions for the ChIP-seq Peak Quantification Pipeline.

This module provides functions to:
1. Find BAM files and extract sample information.
2. Run samtools flagstat to get total mapped reads.
3. Run featureCounts to quantify reads in peaks/regions.
4. Calculate Counts Per Million (CPM) and log-transform them.
"""
import os
import logging
import fnmatch
import natsort
import subprocess
import re
import pandas as pd
import numpy as np
from typing import Dict, Optional, List, Tuple, TypeVar
import tempfile

K = TypeVar('K')
V = TypeVar('V')


def get_sample_data(
    bam_folder: str,
    logger: logging.Logger,
    sample_name_delimiter: Optional[str] = None,
    filename_pattern: str = "*.bam"
) -> Tuple[Dict[str, str], List[str]]:
    """
    Retrieves BAM files from a specified folder that match a filename pattern,
    constructs a dictionary mapping sample names to their full file paths,
    and returns both the dictionary and the list of processed file paths,
    both sorted naturally using the 'natsort' library.

    Sample names are derived from filenames. If 'sample_name_delimiter' is provided
    and found in the filename, the part before the delimiter is used as the sample name.
    Otherwise, the filename without the '.bam' extension is used.

    Args:
        bam_folder (str): Directory containing BAM files.
        logger (logging.Logger): Logger instance for logging messages.
        sample_name_delimiter (Optional[str]): A string that acts as a delimiter
            indicating the end of the sample name in the filename. The part of the
            filename *before* this delimiter will be used as the sample name.
            If None or not found, the filename without '.bam' is used. Defaults to None.
        filename_pattern (str): A wildcard pattern to filter filenames
                                  (e.g., "T*H3K27ac*bam"). Defaults to "*.bam" (all .bam files).

    Returns:
        Tuple[Dict[str, str], List[str]]:
            A tuple containing:
            - dict: A dictionary mapping sample names to BAM file paths, sorted naturally by sample name.
            - list: A list of the full paths of the processed BAM files, sorted naturally by filename.
            Returns ({}, []) if no matching BAM files are found.
    """
    # 処理開始のログ
    logger.info(f"Scanning BAM folder: {bam_folder}, filename pattern: '{filename_pattern}'")
    # 結果を格納する辞書とリストを初期化
    sample_dict = {}
    processed_files_list = []
    bam_file_count = 0 # 見つかったBAMファイルの数をカウント

    try:
        # フォルダ内の全アイテムをリストアップ
        all_files = os.listdir(bam_folder)
        # デバッグログ
        logger.debug(f"Found {len(all_files)} items in the directory.")

        # 各アイテムを処理
        for file in all_files:
            # ファイル名がパターンにマッチし、拡張子が.bamかチェック
            if file.endswith(".bam") and fnmatch.fnmatch(file, filename_pattern):
                # デバッグログ
                logger.debug(f"File '{file}' matches pattern '{filename_pattern}'. Processing...")
                bam_file_count += 1
                # 完全なファイルパスを作成
                file_path = os.path.join(bam_folder, file)
                # 処理対象リストに追加
                processed_files_list.append(file_path)

                # --- サンプル名の抽出処理 ---
                sample_name = ""
                # sample_name_delimiterが指定され、ファイル名に含まれるか？
                if sample_name_delimiter and sample_name_delimiter in file:
                    try:
                        # 区切り文字で分割し、最初の部分をサンプル名とする
                        sample_name = file.split(sample_name_delimiter)[0]
                        # デバッグログ
                        logger.debug(f"Extracted sample name '{sample_name}' from '{file}' using delimiter '{sample_name_delimiter}'.")
                    except Exception as e:
                        # 分割エラー時のフォールバック
                        logger.warning(f"Could not extract sample name from '{file}' using delimiter '{sample_name_delimiter}': {e}. Using filename base as key.")
                        sample_name = os.path.splitext(file)[0] # 拡張子を除いたファイル名
                else:
                    # 区切り文字がない、または見つからない場合
                    sample_name = os.path.splitext(file)[0] # 拡張子を除いたファイル名
                    if sample_name_delimiter:
                         # デバッグログ
                         logger.debug(f"Delimiter '{sample_name_delimiter}' not found in '{file}'. Using filename base '{sample_name}' as sample name.")
                    else:
                         # デバッグログ
                         logger.debug(f"No delimiter specified. Using filename base '{sample_name}' as sample name for '{file}'.")
                # --- サンプル名抽出ここまで ---

                # 辞書にサンプル名とパスを追加（重複チェック含む）
                if sample_name in sample_dict:
                    # 警告ログ
                    logger.warning(f"Duplicate sample name '{sample_name}' detected. Overwriting path from '{sample_dict[sample_name]}' to '{file_path}'. Ensure delimiters uniquely identify samples.")
                sample_dict[sample_name] = file_path

    # --- エラーハンドリング ---
    except FileNotFoundError:
        # エラーログ
        logger.error(f"BAM folder not found: {bam_folder}")
        raise # 例外を呼び出し元に投げる
    except Exception as e:
        # エラーログ
        logger.error(f"An error occurred while processing files in {bam_folder}: {e}", exc_info=True)
        raise # 例外を呼び出し元に投げる
    # --- エラーハンドリングここまで ---

    # マッチするBAMファイルが見つからなかった場合
    if not sample_dict:
        # 警告ログ
        logger.warning(f"No BAM files matching pattern '{filename_pattern}' found in folder: {bam_folder}")
        return {}, [] # 空の辞書とリストを返す

    # --- natsort を使用したソート処理 ---
    # 情報ログ
    logger.info(f"Found {bam_file_count} BAM files matching the pattern. Sorting sample data naturally using natsort.")

    # 辞書をキー（サンプル名）で自然順ソート
    sorted_sample_dict = dict(natsort.natsorted(sample_dict.items()))

    # ファイルパスリストをファイル名部分で自然順ソート
    sorted_processed_files_list = natsort.natsorted(processed_files_list, key=os.path.basename)
    # --- ソート処理ここまで ---

    # デバッグログ
    logger.debug(f"Sorted sample dictionary ({len(sorted_sample_dict)} entries) and file list ({len(sorted_processed_files_list)} entries) created.")

    # ソート済みの辞書とリストをタプルで返す
    return sorted_sample_dict, sorted_processed_files_list


# === featureCounts 実行用ヘルパー関数 ===
def _save_featurecounts_log(log_path: str, stdout_content: Optional[str], stderr_content: Optional[str], logger: logging.Logger):
    """Saves stdout and stderr content from featureCounts to the specified log file."""
    # ログ保存試行のログ
    logger.info(f"Attempting to save featureCounts log to: {log_path}")
    try:
        # 保存先ディレクトリが存在しない場合は作成
        log_dir = os.path.dirname(log_path)
        if log_dir and not os.path.exists(log_dir):
            os.makedirs(log_dir, exist_ok=True)
            # ディレクトリ作成ログ
            logger.info(f"Created directory for log file: {log_dir}")

        # ファイルに書き込み (encoding指定)
        with open(log_path, 'w', encoding='utf-8') as f:
            f.write("--- featureCounts stdout ---\n")
            # None や空文字列の場合に対応
            f.write((stdout_content if stdout_content else "[No stdout captured or content was empty]") + "\n")
            f.write("\n--- featureCounts stderr ---\n")
            # None や空文字列の場合に対応
            f.write((stderr_content if stderr_content else "[No stderr captured or content was empty]") + "\n")
        # ログ保存成功ログ
        logger.info(f"Successfully saved featureCounts log to: {log_path}")
    except Exception as e:
        # ログ保存失敗時のエラーログ
        # ログ保存の失敗は run_featurecounts の成否には影響させない
        logger.error(f"Failed to save featureCounts log to {log_path}: {e}")


def run_featurecounts(
    bam_files: List[str],
    saf_file: str,
    logger: logging.Logger,
    output_file: Optional[str] = None,
    log_file: Optional[str] = None,
    threads: int = 1,
    is_paired_end: bool = True,
    featurecounts_path: Optional[str] = 'featureCounts',
    additional_options: Optional[List[str]] = None
) -> Optional[str]:
    """
    Runs featureCounts using subprocess for the given BAM files using a SAF annotation file.
    Optionally saves execution logs. Defaults to paired-end mode.

    Args:
        bam_files (List[str]): List of full paths to input BAM files.
        saf_file (str): Path to the annotation file in SAF (Simplified Annotation Format).
        logger (logging.Logger): Logger object for logging messages.
        output_file (Optional[str]): Path for the output counts file.
            If None, defaults to '<saf_file_basename>_featureCounts.txt'.
            Defaults to None.
        log_file (Optional[str]): Path to save the stdout and stderr from featureCounts.
            If None, logs are not saved to a separate file. Defaults to None.
        threads (int): Number of threads to use (-T option). Defaults to 1.
        is_paired_end (bool): Specifies if reads are paired-end (-p option is added if True).
            Set to False for single-end reads. Defaults to True (paired-end mode).
        featurecounts_path (Optional[str]): Full path to the featureCounts executable.
            Defaults to 'featureCounts'.
        additional_options (Optional[List[str]]): List of additional arguments
            to pass to featureCounts. Defaults to None.

    Returns:
        Optional[str]: The full path to the generated output counts file if successful,
                       None otherwise.
    """
    # --- 入力値の検証 ---
    # BAMファイルリストが空でないか確認
    if not bam_files:
        logger.error("No BAM files provided to run_featurecounts.")
        return None
    # SAFファイルが存在するか確認
    if not os.path.exists(saf_file):
        logger.error(f"SAF annotation file not found: {saf_file}")
        return None
    # SAFファイルの拡張子チェック (オプション、警告のみ)
    if not saf_file.lower().endswith(('.saf', '.txt')): # 一般的な拡張子か確認
        logger.warning(f"Input annotation file '{os.path.basename(saf_file)}' does not have a typical .saf or .txt extension. Ensure it is a valid SAF format file.")
    # featurecounts_pathのデフォルト設定
    if not featurecounts_path:
        featurecounts_path = 'featureCounts'

    # --- 出力ファイルパスの決定 ---
    effective_output_file = ""
    if output_file is None:
        # デフォルトの出力ファイル名を生成 (SAFファイル名ベース)
        saf_basename = os.path.splitext(os.path.basename(saf_file))[0]
        effective_output_file = f"{saf_basename}_featureCounts.txt"
        logger.info(f"Output file path not specified, defaulting to: {effective_output_file}")
    else:
        effective_output_file = output_file
        # 出力ディレクトリ確認・作成
        output_dir = os.path.dirname(effective_output_file)
        if output_dir and not os.path.exists(output_dir):
            try:
                os.makedirs(output_dir, exist_ok=True)
                logger.info(f"Created output directory: {output_dir}")
            except Exception as e:
                logger.error(f"Failed to create output directory {output_dir}: {e}")
                return None

    # --- featureCounts コマンドリストの構築 ---
    cmd = [featurecounts_path]
    # 基本オプション
    cmd.extend(['-T', str(threads)])
    cmd.extend(['-a', saf_file]) # SAFファイルを指定
    cmd.extend(['-F', 'SAF'])   # 形式をSAFと明示的に指定
    cmd.extend(['-o', effective_output_file])
    # ペアエンドオプション (-p は is_paired_end が True の場合のみ追加)
    if is_paired_end:
        cmd.append('-p')
        # ログでモードを明確化
        logger.info("Paired-end mode enabled (-p option added).")
    else:
        # ログでモードを明確化
        logger.info("Single-end mode enabled (no -p option).")
    # 追加オプション
    if additional_options:
        cmd.extend(additional_options)
    # 入力BAMファイルリスト
    cmd.extend(bam_files)
    logger.debug(f"Constructed featureCounts command list: {cmd}")

    # --- サブプロセスとして featureCounts を実行 ---
    # 実行開始ログ
    logger.info(f"Running featureCounts for {len(bam_files)} BAM files...")
    # 実行するコマンド文字列をログに出力 - INFOからDEBUGに変更
    logger.debug(f"Command: {' '.join(cmd)}")

    stdout_log: Optional[str] = None
    stderr_log: Optional[str] = None

    try:
        # コマンド実行
        result = subprocess.run(cmd, check=True, capture_output=True, text=True, encoding='utf-8')
        # 実行結果のstdout/stderrを取得
        stdout_log = result.stdout
        stderr_log = result.stderr
        # 正常終了ログ
        logger.info("featureCounts finished successfully.")
        # stderrのサマリー情報をログ出力 - INFOからDEBUGに変更
        if stderr_log:
            logger.debug(f"featureCounts Summary (from stderr):\n{stderr_log}")
        # stdoutがあればデバッグログ出力
        if stdout_log:
             logger.debug(f"featureCounts stdout:\n{stdout_log}")

        # ログファイルが指定されていれば保存
        if log_file:
            _save_featurecounts_log(log_file, stdout_log, stderr_log, logger) # ヘルパー関数呼び出し

        # 成功したので出力ファイルの絶対パスを返す
        return os.path.abspath(effective_output_file)

    # --- エラーハンドリング ---
    except FileNotFoundError:
        # コマンドが見つからない場合
        logger.error(f"Error: '{featurecounts_path}' command not found.")
        logger.error("Please ensure featureCounts is installed and in your PATH, or provide the correct full path.")
        return None
    except subprocess.CalledProcessError as e:
        # コマンドがエラー終了した場合
        stdout_log = e.stdout
        stderr_log = e.stderr
        logger.error(f"featureCounts failed with exit code {e.returncode}")
        if stdout_log: logger.error(f"Captured stdout on failure:\n{stdout_log}")
        if stderr_log: logger.error(f"Captured stderr on failure:\n{stderr_log}")
        # エラー時もログ保存試行
        if log_file:
            _save_featurecounts_log(log_file, stdout_log, stderr_log, logger)
        # 不完全ファイルの削除試行
        if os.path.exists(effective_output_file):
            logger.warning(f"Attempting to remove potentially incomplete output file: {effective_output_file}")
            try: os.remove(effective_output_file)
            except OSError as rm_err: logger.error(f"Failed to remove incomplete output file {effective_output_file}: {rm_err}")
        return None
    except Exception as e:
        # その他の予期せぬエラー
        logger.error(f"An unexpected error occurred while running featureCounts: {e}", exc_info=True)
        # ログ保存試行
        if log_file:
             _save_featurecounts_log(log_file, stdout_log, stderr_log, logger)
        return None
    # --- エラーハンドリングここまで ---


def run_samtools_flagstat(
    bam_files: List[str],
    output_dir: str,
    logger: logging.Logger,
    samtools_path: str = 'samtools',
    threads: int = 1
) -> Dict[str, int]:
    """
    Runs 'samtools flagstat' on a list of BAM files, saves the output,
    and extracts the number of mapped reads.

    Args:
        bam_files (List[str]): List of full paths to input BAM files.
        output_dir (str): Directory to save the flagstat output files.
        logger (logging.Logger): Logger object for logging messages.
        samtools_path (str): Path to the samtools executable. Defaults to 'samtools'.
        threads (int): Number of threads for samtools (-@ option). Defaults to 1.

    Returns:
        Dict[str, int]: A dictionary mapping input BAM file paths to the
                        number of mapped reads extracted from flagstat output.
                        Returns an empty dictionary if no files are processed or
                        no mapped reads are found. Files with errors during
                        processing or parsing will be excluded.
    """
    # 出力ディレクトリが存在しない場合は作成
    try:
        os.makedirs(output_dir, exist_ok=True)
        logger.info(f"Ensured flagstat output directory exists: {output_dir}")
    except OSError as e:
        logger.error(f"Failed to create output directory {output_dir}: {e}")
        return {} # ディレクトリ作成失敗時は空辞書を返す

    # mappedリード数を格納するための辞書を初期化
    mapped_reads_dict: Dict[str, int] = {}
    # mapped reads数を抽出するための正規表現パターンをコンパイル
    # 例: "39502262 + 0 mapped (100.00% : N/A)" の先頭の数値を取得
    mapped_pattern = re.compile(r'^(\d+) \+ \d+ mapped \(') # QC-passed reads の行にマッチ

    # 入力されたBAMファイルリストをループ処理
    total_files = len(bam_files)
    processed_count = 0
    error_count = 0
    parse_error_count = 0

    logger.info(f"Starting samtools flagstat processing for {total_files} BAM files.")

    for i, bam_path in enumerate(bam_files):
        # 処理中のファイル情報をログ出力
        base_name = os.path.basename(bam_path)
        logger.info(f"Processing file {i+1}/{total_files}: {base_name}")

        # flagstatの出力ファイル名を決定 (元のBAM名 + .flagstat.txt)
        flagstat_filename = os.path.splitext(base_name)[0] + ".flagstat.txt"
        output_flagstat_path = os.path.join(output_dir, flagstat_filename)

        # samtools flagstat コマンドリストを構築
        cmd = [samtools_path, 'flagstat', '-@', str(threads), bam_path]
        logger.debug(f"Running command: {' '.join(cmd)}")

        try:
            # コマンド実行、出力をキャプチャ
            result = subprocess.run(cmd, check=True, capture_output=True, text=True, encoding='utf-8')
            flagstat_output = result.stdout # 標準出力を取得

            # --- flagstat 出力をファイルに保存 ---
            try:
                with open(output_flagstat_path, 'w', encoding='utf-8') as f:
                    f.write(flagstat_output)
                logger.info(f"Saved flagstat output to: {output_flagstat_path}")
            except Exception as e:
                # ファイル保存失敗は警告に留め、抽出処理は続行
                logger.warning(f"Failed to save flagstat output for {base_name} to {output_flagstat_path}: {e}")
            # --- ファイル保存ここまで ---

            # --- mapped reads数の抽出 ---
            mapped_reads = -1 # 初期値（見つからなかったことを示す）
            found_mapped_line = False
            for line in flagstat_output.splitlines():
                match = mapped_pattern.match(line.strip()) # 行頭の空白を削除してマッチ試行
                if match:
                    found_mapped_line = True
                    try:
                        # マッチした最初のグループ（数字）を整数に変換
                        mapped_reads = int(match.group(1))
                        logger.info(f"Extracted mapped reads for {base_name}: {mapped_reads}")
                        # 成功したら辞書に追加 (キーはBAMフルパス)
                        mapped_reads_dict[bam_path] = mapped_reads
                        processed_count += 1
                        break # mapped行が見つかればこのファイルの処理は終了
                    except ValueError:
                        logger.error(f"Could not convert extracted mapped reads '{match.group(1)}' to integer for {base_name}.")
                        parse_error_count += 1
                        break # 変換エラーでも終了
                    except Exception as e:
                         logger.error(f"Unexpected error during mapped reads conversion for {base_name}: {e}")
                         parse_error_count += 1
                         break

            # mapped行が見つからなかった場合の警告
            if not found_mapped_line:
                 logger.warning(f"Could not find the 'mapped' line in flagstat output for {base_name}.")
                 parse_error_count += 1
            # --- mapped reads抽出ここまで ---

        # --- samtools実行自体のエラーハンドリング ---
        except FileNotFoundError:
            logger.error(f"Error: '{samtools_path}' command not found. Skipping {base_name}.")
            error_count += 1
            # samtoolsがない場合、これ以上処理できない可能性があるので続けるか判断が必要
            # ここでは個々のファイルをスキップ
        except subprocess.CalledProcessError as e:
            logger.error(f"samtools flagstat failed for {base_name} with exit code {e.returncode}.")
            # エラー内容をログに出力 (stderr)
            if e.stderr: logger.error(f"Stderr:\n{e.stderr}")
            error_count += 1
        except Exception as e:
            logger.error(f"An unexpected error occurred while processing {base_name}: {e}", exc_info=True)
            error_count += 1
        # --- エラーハンドリングここまで ---

    # 処理終了のサマリーログ
    logger.info(f"Flagstat processing finished.")
    logger.info(f"Successfully extracted mapped reads for: {processed_count}/{total_files} files.")
    if error_count > 0: logger.warning(f"Errors during samtools execution: {error_count} files.")
    if parse_error_count > 0: logger.warning(f"Errors during parsing flagstat output: {parse_error_count} files.")

    # 結果の辞書 (BAMパス -> mappedリード数) を返す
    return mapped_reads_dict


def calculate_logcpm(
    counts_df: pd.DataFrame,
    total_reads_dict: Dict[str, int],
    logger: logging.Logger,
    bam_path_to_sample_name: Optional[Dict[str, str]] = None,
    add_peakid_column: bool = False,
    id_column_name: str = 'Geneid',
    output_id_column_name: str = 'PeakID',
    log_transform_base: Optional[int] = 2,
    pseudocount: float = 1.0
) -> Optional[pd.DataFrame]:
    """
    Calculates log-transformed CPM from a counts DataFrame and total mapped reads.
    Formats the output table with a specified ID column name and sample names.
    The 'Start' coordinate in the output DataFrame is converted to 0-based (BED-like).

    Args:
        counts_df (pd.DataFrame): DataFrame containing counts data, typically from
                                  featureCounts output read into pandas. Expected columns include
                                  'id_column_name', 'Chr', 'Start' (1-based from featureCounts), 'End', etc.,
                                  and count columns named with full BAM paths (usually starting
                                  from the 7th column).
        total_reads_dict (Dict[str, int]): Dictionary mapping BAM file full paths
                                           to total mapped reads. Keys must correspond to
                                           count column names in counts_df.
        logger (logging.Logger): Logger object.
        bam_path_to_sample_name (Optional[Dict[str, str]]): Dictionary mapping BAM file full
                                           paths to desired final sample names (e.g., 'T1').
                                           If provided, count columns will be renamed.
        add_peakid_column (bool): If True, replace original IDs with sequential PeakIDs.
                                  The ID column name will be 'output_id_column_name'.
        id_column_name (str): Name of the feature identifier column in the input DataFrame.
        output_id_column_name (str): Desired name for the ID column in the output DataFrame.
        log_transform_base (Optional[int]): Base for log transformation (2, 10, or None).
                                            If None, returns CPM values without log transformation.
        pseudocount (float): Pseudocount added before log transformation.

    Returns:
        Optional[pd.DataFrame]: A DataFrame with ID (as 'output_id_column_name'), Chr, Start (0-based), End,
                                and log-transformed CPM values. Returns None on error.
    """
    # 処理の種類を決定（ログ用）
    transform_type = f"log{log_transform_base}-CPM" if log_transform_base else "CPM"
    logger.info(f"Calculating {transform_type} from provided DataFrame.")

    try:
        # --- 入力DataFrameの検証 ---
        # DataFrameが有効か、空でないかチェック
        if not isinstance(counts_df, pd.DataFrame) or counts_df.empty:
            logger.error("Input counts_df is not a valid or non-empty DataFrame.")
            return None

        logger.info(f"Processing DataFrame with {len(counts_df)} features and {len(counts_df.columns)} total columns.")

        # 必要なカラム名の検証
        if id_column_name not in counts_df.columns:
             logger.error(f"Specified ID column '{id_column_name}' not found in input DataFrame columns: {counts_df.columns.tolist()}")
             return None
        # 出力に必要なアノテーション列
        required_annot_cols = [id_column_name, 'Chr', 'Start', 'End']
        missing_cols = [col for col in required_annot_cols if col not in counts_df.columns]
        if missing_cols:
            logger.error(f"Required annotation columns are missing from the input DataFrame: {missing_cols}")
            return None

        # カウントカラムを特定 (featureCounts標準出力の7列目以降と仮定)
        if len(counts_df.columns) < 7:
            logger.error("Input DataFrame does not seem to have the expected featureCounts format (needs at least 7 columns: ID, Chr, Start, End, Strand, Length, Counts...).")
            return None
        count_columns = counts_df.columns[6:].tolist() # BAMパス名のカラムリスト

        logger.info(f"Identified {len(count_columns)} potential count columns.")
        # --- 検証ここまで ---

        # --- 結果用DataFrameの準備 ---
        # 出力に必要なアノテーション列のみをコピー
        result_df = counts_df[required_annot_cols].copy()
        # --- 準備ここまで ---

        # --- IDカラム処理 ---
        # 最終的なIDカラム名は output_id_column_name に統一する
        if add_peakid_column:
            # 連番ID (例: PeakID_000001) を生成してID列を置き換える
            logger.info(f"Replacing '{id_column_name}' with sequential IDs named '{output_id_column_name}'.")
            num_rows = len(result_df)
            if num_rows > 0:
                num_digits = len(str(num_rows)) # 桁数を決定
                # 新しいIDリストを生成
                new_ids = [f"{output_id_column_name}_{i:0{num_digits}d}" for i in range(1, num_rows + 1)]
                # 新しいID列を追加
                result_df[output_id_column_name] = new_ids
                # 元のID列が新しいID列名と異なる場合は削除
                if id_column_name != output_id_column_name and id_column_name in result_df.columns:
                    result_df = result_df.drop(columns=[id_column_name])
            else: # DataFrameが空の場合
                 logger.warning("DataFrame is empty, cannot add sequential PeakID column.")
                 # 空でもカラムは作成しておく
                 if id_column_name != output_id_column_name and id_column_name in result_df.columns:
                     result_df = result_df.drop(columns=[id_column_name])
                 if output_id_column_name not in result_df.columns:
                      result_df[output_id_column_name] = pd.NA
        elif id_column_name != output_id_column_name:
            # add_peakid_column=False でも、カラム名を指定の出力名に変更する
            logger.info(f"Renaming ID column '{id_column_name}' to '{output_id_column_name}'.")
            result_df = result_df.rename(columns={id_column_name: output_id_column_name})
        final_id_col = output_id_column_name # 最終的なID列名
        # --- IDカラム処理ここまで ---

        # --- CPM および LogCPM 計算 ---
        processed_samples_count = 0
        skipped_samples = []
        calculated_columns_map: Dict[str, str] = {} # {BAMパス: 一時カラム名}

        logger.info(f"Calculating {transform_type} values...")
        # カウントカラム（BAMファイルパス）をループ
        for bam_path_col in count_columns:
            # DataFrameに実際にそのカラムが存在するか念のため確認
            if bam_path_col not in counts_df.columns:
                 logger.warning(f"Count column '{bam_path_col}' (expected from header) not found in DataFrame. Skipping.")
                 skipped_samples.append(bam_path_col)
                 continue
            # このBAMパスに対応する総リード数が辞書にあるか確認
            if bam_path_col in total_reads_dict:
                total_reads = total_reads_dict[bam_path_col]
                # 一時的な計算後カラム名を生成 (BAMファイル名ベース)
                temp_col_name = os.path.basename(bam_path_col)
                if temp_col_name.endswith('.bam'): temp_col_name = temp_col_name[:-4]

                # 総リード数が0より大きいか確認
                if total_reads > 0:
                    # CPM 計算
                    cpm_values = (counts_df[bam_path_col] / total_reads) * 1_000_000

                    # 対数変換とDataFrameへの追加
                    if log_transform_base == 2:
                        result_df[temp_col_name] = np.log2(cpm_values + pseudocount)
                    elif log_transform_base == 10:
                        result_df[temp_col_name] = np.log10(cpm_values + pseudocount)
                    elif log_transform_base is None: # 対数変換しない場合
                        result_df[temp_col_name] = cpm_values
                    else: # 未サポートの底
                        logger.warning(f"Unsupported log base: {log_transform_base}. Storing raw CPM for {temp_col_name}.")
                        result_df[temp_col_name] = cpm_values

                    calculated_columns_map[bam_path_col] = temp_col_name # 処理したカラムを記録
                    processed_samples_count += 1
                else:
                    # 総リード数が0の場合
                    logger.warning(f"Total mapped reads is 0 for '{bam_path_col}'. Setting column '{temp_col_name}' to NaN.")
                    result_df[temp_col_name] = np.nan # 結果はNaN
                    calculated_columns_map[bam_path_col] = temp_col_name # カラム自体は作成
                    skipped_samples.append(bam_path_col)
            else:
                # 総リード数辞書にBAMパスが見つからない場合
                logger.warning(f"Total mapped reads not found for '{bam_path_col}'. Skipping this sample.")
                skipped_samples.append(bam_path_col)
        # --- 計算ここまで ---

        # --- サンプルカラム名の置換 ---
        final_sample_columns = [] # 最終的なサンプルカラム名のリスト
        if bam_path_to_sample_name:
            # 提供されたマッピングを使ってカラム名を変更する
            logger.info("Renaming sample columns using provided mapping.")
            rename_dict = {} # リネーム用辞書 {一時的な名前: 最終的な名前}
            processed_bam_paths = set() # 処理済みBAMパスを追跡

            # 計算済みの一時カラム名に対してループ
            for bam_path, temp_col_name in calculated_columns_map.items():
                # このBAMパスに対応する最終サンプル名があるか？
                if bam_path in bam_path_to_sample_name:
                    final_sample_name = bam_path_to_sample_name[bam_path]
                    rename_dict[temp_col_name] = final_sample_name
                    final_sample_columns.append(final_sample_name)
                    processed_bam_paths.add(bam_path)
                else:
                    # マッピングがない場合は一時的な名前をそのまま使う
                    logger.warning(f"No sample name mapping found for BAM path '{bam_path}'. Using derived column name '{temp_col_name}'.")
                    final_sample_columns.append(temp_col_name) # リネームしない

            # DataFrameのカラム名をリネーム
            result_df = result_df.rename(columns=rename_dict)
            logger.info(f"Renamed {len(rename_dict)} sample columns based on mapping.")

            # マッピング辞書に含まれていたのに、対応する計算済みカラムがなかった場合（警告）
            unmapped_samples = set(bam_path_to_sample_name.keys()) - processed_bam_paths - set(skipped_samples)
            if unmapped_samples:
                logger.warning(f"Provided mapping contained entries for BAM paths that were not found or skipped: {unmapped_samples}")
        else:
            # マッピングが提供されない場合は、一時的な名前が最終的な名前になる
            final_sample_columns = list(calculated_columns_map.values())
            logger.info("No sample name mapping provided. Using derived column names.")
        # --- カラム名置換ここまで ---

        logger.info(f"Finished calculation for {processed_samples_count} samples.")
        if skipped_samples:
            logger.warning(f"Skipped {len(skipped_samples)} samples due to missing or zero total reads.")

        # --- 最終的なカラム順序の決定と選択 ---
        # IDカラム名を取得 (rename後、または元のまま)
        final_annot_cols = [final_id_col, 'Chr', 'Start', 'End']
        # 存在しないアノテーションカラムを除外 (例: ID列を削除した場合など)
        final_annot_cols = [col for col in final_annot_cols if col in result_df.columns]
        # 最終的なカラムリスト = アノテーション + サンプル (ソートされたサンプル名の順になるはず)
        final_columns_ordered = final_annot_cols + final_sample_columns
        # 存在しないカラムが final_sample_columns に含まれる可能性は低いが、念のためチェック
        final_columns_ordered = [col for col in final_columns_ordered if col in result_df.columns]
        # DataFrameを最終的なカラム順序で選択
        result_df = result_df[final_columns_ordered]
        # --- カラム順序調整ここまで ---

        # ★★★ Start座標を 0-based (BED形式) に戻す ★★★
        if 'Start' in result_df.columns:
            # Start列が数値型であることを確認 (エラーを防ぐため)
            if pd.api.types.is_numeric_dtype(result_df['Start']):
                logger.debug("Converting 'Start' column to 0-based (BED format) by subtracting 1.")
                # .loc を使う (SettingWithCopyWarningを避ける)
                result_df.loc[:, 'Start'] = result_df.loc[:, 'Start'] - 1
            else:
                logger.warning("Could not convert 'Start' column to 0-based because it is not numeric.")
        else:
            logger.warning("Could not find 'Start' column in the result DataFrame to convert to 0-based.")
        # ★★★ 座標変換ここまで ★★★

        # 整形済みのDataFrameを返す
        return result_df

    # --- エラーハンドリング ---
    except ImportError: # pandas, numpy がない場合
        logger.error("Pandas or NumPy library not found. Please install them (`pip install pandas numpy`).")
        return None
    except KeyError as e: # DataFrameカラムや辞書キーが存在しない場合
        logger.error(f"KeyError during processing: {e}. Check column names in DataFrame and keys in total_reads_dict.")
        # エラー発生時にDataFrameや辞書の情報をログに出力するとデバッグに役立つ
        if 'counts_df' in locals() and isinstance(counts_df, pd.DataFrame): logger.error(f"Input DataFrame columns: {counts_df.columns.tolist()}")
        if isinstance(total_reads_dict, dict): logger.error(f"Sample of total_reads_dict keys: {list(total_reads_dict.keys())[:5]}")
        return None
    except Exception as e: # その他の予期せぬエラー
        logger.error(f"An unexpected error occurred during logCPM calculation: {e}", exc_info=True)
        return None
    # --- エラーハンドリングここまで ---


# === invert_dictionary ヘルパー関数の定義 ===
def invert_dictionary(d: Dict[K, V]) -> Dict[V, K]:
    """
    Inverts a dictionary, swapping keys and values.
    Assumes values are unique and hashable to be used as keys.

    Args:
        d (Dict[K, V]): The dictionary to invert.

    Returns:
        Dict[V, K]: A new dictionary with keys and values swapped.
    """
    # 辞書内包表記を使った簡潔版 (重複は後勝ちで上書き)
    return {v: k for k, v in d.items()}


# === BED to SAF 変換関数 ===
def convert_bed_to_saf(bed_path: str, logger: logging.Logger) -> Optional[str]:
    """
    Converts a BED file to a temporary SAF file required by featureCounts.
    Handles standard 3-column BED and attempts to use columns 4 (name) and 6 (strand) if present.
    Generates sequential IDs if name column is missing or unusable. Sets strand to '.' if missing.
    Adjusts start coordinate from 0-based (BED) to 1-based (SAF).

    Args:
        bed_path (str): Path to the input BED file.
        logger (logging.Logger): Logger object.

    Returns:
        Optional[str]: The path to the created temporary SAF file, or None if conversion fails.
                       The caller is responsible for deleting this temporary file later using os.remove().
    """
    # 処理開始ログ
    logger.info(f"Converting BED file to temporary SAF: {os.path.basename(bed_path)}")
    try:
        # --- BEDファイルの読み込み ---
        try:
            # ヘッダーなし、タブ区切りで読み込み
            bed_df = pd.read_csv(bed_path, sep='\t', header=None, dtype=str) # dtype=str で安全に読み込む
        except pd.errors.EmptyDataError:
            logger.error(f"Input BED file is empty or could not be parsed: {bed_path}")
            return None
        except FileNotFoundError:
             logger.error(f"Input BED file not found: {bed_path}")
             return None
        except Exception as e:
             logger.error(f"Failed to read BED file {bed_path}: {e}")
             return None

        # 列数の確認
        num_cols = bed_df.shape[1]
        logger.debug(f"Read BED file with {len(bed_df)} rows and {num_cols} columns.")
        if num_cols < 3:
            logger.error(f"BED file has fewer than 3 required columns (Chr, Start, End): {bed_path}")
            return None
        # --- BED読み込み完了 ---

        # --- SAFデータフレームの準備 ---
        saf_data = {} # SAF列データを格納

        # 必須列: Chr, Start (1-based), End
        try:
            saf_data['Chr'] = bed_df.iloc[:, 0] # 1列目
            # Start: BEDの2列目(0-based) + 1
            saf_data['Start'] = pd.to_numeric(bed_df.iloc[:, 1], errors='coerce') + 1
            # End: BEDの3列目
            saf_data['End'] = pd.to_numeric(bed_df.iloc[:, 2], errors='coerce')
            # 数値変換失敗チェック
            if saf_data['Start'].isnull().any() or saf_data['End'].isnull().any():
                 logger.error("Could not convert Start/End columns to numeric. Check BED file format.")
                 return None
            saf_data['Start'] = saf_data['Start'].astype(int)
            saf_data['End'] = saf_data['End'].astype(int)
        except Exception as e:
            logger.error(f"Error processing required BED columns (Chr, Start, End): {e}. Check format.")
            return None

        # オプション列: GeneID (4列目があれば使用)
        if num_cols >= 4 and not bed_df.iloc[:, 3].isnull().all():
             logger.debug("Using column 4 from BED as GeneID.")
             saf_data['GeneID'] = bed_df.iloc[:, 3]
             if saf_data['GeneID'].duplicated().any():
                 logger.warning(f"Duplicate IDs found in column 4 of BED file {bed_path}. This may affect featureCounts results if IDs are used for merging.")
        else:
             # 連番生成
             logger.debug("Generating sequential GeneIDs (peak_1, peak_2...).")
             num_rows = len(bed_df)
             num_digits = len(str(num_rows)) if num_rows > 0 else 1
             saf_data['GeneID'] = [f"peak_{i:0{num_digits}d}" for i in range(1, num_rows + 1)]

        # オプション列: Strand (6列目があれば使用)
        if num_cols >= 6 and not bed_df.iloc[:, 5].isnull().all():
             logger.debug("Using column 6 from BED as Strand.")
             saf_data['Strand'] = bed_df.iloc[:, 5].copy() # SettingWithCopyWarning回避のためコピー
             # '+' または '-' 以外を '.' に置換
             valid_strands = ['+', '-']
             invalid_strand_mask = ~saf_data['Strand'].isin(valid_strands)
             if invalid_strand_mask.any():
                 logger.warning(f"Found invalid strand values (not '+' or '-') in column 6 of {bed_path}. Replacing with '.'.")
                 saf_data['Strand'].loc[invalid_strand_mask] = '.'
        else:
             # デフォルト '.'
             logger.debug("Setting Strand to '.' (default).")
             saf_data['Strand'] = '.'

        # SAF DataFrameを作成 (列順序指定)
        try:
            saf_df = pd.DataFrame(saf_data)[['GeneID', 'Chr', 'Start', 'End', 'Strand']]
        except KeyError as e:
             logger.error(f"Failed to create SAF DataFrame. Missing column data: {e}")
             return None
        # --- SAFデータフレーム準備完了 ---

        # --- 一時SAFファイルの作成と書き込み ---
        temp_saf_file = None
        try:
            # delete=False でパスを保持、suffixで拡張子指定
            temp_saf_file = tempfile.NamedTemporaryFile(mode='w+', suffix='.saf', delete=False, encoding='utf-8')
            # ヘッダー付き、タブ区切りで書き込み
            saf_df.to_csv(temp_saf_file.name, sep='\t', index=False, header=True)
            temp_saf_path = temp_saf_file.name
            temp_saf_file.close() # ファイルは閉じてもパスは有効
            logger.info(f"Successfully converted BED to temporary SAF file: {temp_saf_path}")
            # 一時ファイルのパスを返す
            return temp_saf_path
        except Exception as e:
            logger.error(f"Failed to create or write temporary SAF file: {e}")
            # 作成途中だったファイルをクリーンアップ
            if temp_saf_file and os.path.exists(temp_saf_file.name):
                try:
                    temp_saf_file.close()
                    os.remove(temp_saf_file.name)
                    logger.debug(f"Cleaned up temporary file on error: {temp_saf_file.name}")
                except Exception as cleanup_e:
                     logger.error(f"Error during temporary file cleanup: {cleanup_e}")
            return None
        # --- 一時ファイル作成完了 ---

    except Exception as e:
        # 予期せぬエラー
        logger.error(f"An unexpected error occurred during BED to SAF conversion for {bed_path}: {e}", exc_info=True)
        return None


def convert_mergesv_to_saf(merge_sv_path: str, logger: logging.Logger) -> Optional[str]:
    """
    Converts a merge_SV.tsv file (expecting 'se_data' column like chr_start_end,
    0-based start) to a temporary SAF file (1-based start).

    Args:
        merge_sv_path (str): Path to the input merge_SV.tsv file.
        logger (logging.Logger): Logger object.

    Returns:
        Optional[str]: Path to the temporary SAF file, or None on error.
                       Caller must delete the temporary file.
    """
    # 処理開始ログ (英語)
    logger.info(f"Converting merge_SV format file to temporary SAF: {os.path.basename(merge_sv_path)}")
    try:
        # --- ファイル読み込みとヘッダー確認 ---
        try:
            # ヘッダーあり、タブ区切りで読み込み
            df = pd.read_csv(merge_sv_path, sep='\t', header=0)
        except pd.errors.EmptyDataError:
            logger.error(f"Input merge_SV file is empty: {merge_sv_path}")
            return None
        except FileNotFoundError:
             logger.error(f"Input merge_SV file not found: {merge_sv_path}")
             return None
        except Exception as e:
             logger.error(f"Failed to read merge_SV file {merge_sv_path}: {e}")
             return None

        # DataFrameが空でないか再確認
        if df.empty:
             logger.error(f"Input merge_SV file is empty after reading: {merge_sv_path}")
             return None

        # 最初の列名が 'se_data' であることを確認 (大文字小文字区別しない)
        se_data_col_name = df.columns[0] # 実際のカラム名を取得
        if se_data_col_name.lower() != 'se_data':
            logger.error(f"Expected first column to be 'se_data' (case-insensitive) in {merge_sv_path}, but found '{se_data_col_name}'.")
            return None
        logger.debug(f"Read merge_SV file with {len(df)} rows.")
        # --- 読み込み完了 ---

        # --- 座標の抽出とSAFデータ作成 ---
        saf_data = {}
        try:
            # 'se_data' 列から chr, start, end を抽出 ( '_' 区切り想定)
            # n=2 で指定し、startとendが数値であることを保証しやすくする
            coords = df[se_data_col_name].str.split('_', expand=True, n=2)
            if coords.shape[1] < 3: # 分割結果が3列未満ならエラー
                 logger.error(f"Could not split '{se_data_col_name}' column into at least 3 parts (chr, start, end) using '_'. Check format.")
                 return None

            saf_data['GeneID'] = df[se_data_col_name] # 元の se_data を ID として使用
            saf_data['Chr'] = coords[0]
            # start は BED基準(0-based)なので +1 して SAF(1-based) にする
            saf_data['Start'] = pd.to_numeric(coords[1], errors='coerce') + 1
            saf_data['End'] = pd.to_numeric(coords[2], errors='coerce')
            saf_data['Strand'] = '.' # Strand 情報はないので '.'

            # 数値変換チェック
            if saf_data['Start'].isnull().any() or saf_data['End'].isnull().any():
                 num_failed = saf_data['Start'].isnull().sum() + saf_data['End'].isnull().sum()
                 failed_indices = df.index[saf_data['Start'].isnull() | saf_data['End'].isnull()].tolist()
                 logger.error(f"Could not convert extracted Start/End to numeric for {num_failed} rows. Check '{se_data_col_name}' column format (expecting chr_start_end). Problematic row indices (0-based): {failed_indices[:10]}...") # 最初の10件表示
                 return None
            saf_data['Start'] = saf_data['Start'].astype(int)
            saf_data['End'] = saf_data['End'].astype(int)

        except Exception as e:
            # 座標抽出や変換中のエラー
            logger.error(f"Error extracting or converting coordinates from '{se_data_col_name}' column: {e}", exc_info=True)
            return None

        # SAF DataFrame 作成 (列順序指定)
        saf_df = pd.DataFrame(saf_data)[['GeneID', 'Chr', 'Start', 'End', 'Strand']]
        # --- SAFデータ作成完了 ---

        # --- 一時ファイル書き込み ---
        temp_saf_file = None
        try:
            # 一時ファイルを作成 (delete=Falseでパスを保持)
            temp_saf_file = tempfile.NamedTemporaryFile(mode='w+', suffix='.saf', delete=False, encoding='utf-8')
            # ヘッダー付き、タブ区切りで書き込み
            saf_df.to_csv(temp_saf_file.name, sep='\t', index=False, header=True)
            temp_saf_path = temp_saf_file.name
            temp_saf_file.close() # ファイルを閉じる
            logger.info(f"Successfully converted merge_SV to temporary SAF file: {temp_saf_path}")
            # 一時ファイルのパスを返す
            return temp_saf_path
        except Exception as e:
            logger.error(f"Failed to write temporary SAF file: {e}")
            # エラー発生時、もしファイルが作られていたら削除試行
            if temp_saf_file and os.path.exists(temp_saf_file.name):
                 try:
                     temp_saf_file.close()
                     os.remove(temp_saf_file.name)
                 except Exception as cleanup_e: logger.error(f"Cleanup error: {cleanup_e}")
            return None
        # --- 一時ファイル書き込み完了 ---

    except Exception as e:
        # その他の予期せぬエラー
        logger.error(f"An unexpected error occurred during merge_SV to SAF conversion: {e}", exc_info=True)
        return None