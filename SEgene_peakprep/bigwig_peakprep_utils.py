# === bigwig_peakprep_utils.py ===

import os
import logging
import pandas as pd
import numpy as np
from typing import Dict, Optional, List, Tuple, TypeVar
import pyranges as pr 

import os
import sys
import time
import subprocess
import re
import tempfile


def _save_bamcoverage_log(log_path: str, stdout_content: Optional[str], stderr_content: Optional[str], logger: logging.Logger):
    """Saves stdout and stderr content from bamCoverage to the specified log file."""
    # ログ保存試行のログ
    logger.info(f"Attempting to save bamCoverage log to: {log_path}")
    try:
        # 保存先ディレクトリが存在しない場合は作成
        log_dir = os.path.dirname(log_path)
        if log_dir and not os.path.exists(log_dir):
            os.makedirs(log_dir, exist_ok=True)
            # ディレクトリ作成ログ
            logger.info(f"Created directory for log file: {log_dir}")

        # ファイルに書き込み (encoding指定)
        with open(log_path, 'w', encoding='utf-8') as f:
            f.write("--- bamCoverage stdout ---\n")
            # None や空文字列の場合に対応
            f.write((stdout_content if stdout_content else "[No stdout captured or content was empty]") + "\n")
            f.write("\n--- bamCoverage stderr ---\n")
            # None や空文字列の場合に対応
            f.write((stderr_content if stderr_content else "[No stderr captured or content was empty]") + "\n")
        # ログ保存成功ログ
        logger.info(f"Successfully saved bamCoverage log to: {log_path}")
    except Exception as e:
        # ログ保存失敗時のエラーログ
        # ログ保存の失敗は run_bamcoverage の成否には影響させない
        logger.error(f"Failed to save bamCoverage log to {log_path}: {e}")


def _get_bamcoverage_version(bamcoverage_path: str, logger: logging.Logger) -> Optional[str]:
    """Gets the version of bamCoverage by running 'bamCoverage --version'."""
    try:
        # バージョン確認コマンド実行
        version_cmd = [bamcoverage_path, "--version"]
        result = subprocess.run(version_cmd, check=True, capture_output=True, text=True, encoding='utf-8')
        version_info = result.stdout.strip() if result.stdout else "Unknown version"
        logger.debug(f"bamCoverage version: {version_info}")
        return version_info
    except Exception as e:
        logger.warning(f"Failed to get bamCoverage version: {e}")
        return None


def run_bamcoverage(
    bam_file: str,
    output_dir: str,
    logger: logging.Logger,
    output_filename: Optional[str] = None,
    bin_size: int = 50,
    normalize_using: str = "RPGC",
    effective_genome_size: int = 2913022398,
    bamcoverage_path: str = 'bamCoverage',
    threads: int = 1,
    log_file: Optional[str] = None,
    additional_options: Optional[List[str]] = None
) -> Optional[str]:
    """
    Runs bamCoverage from deeptools to generate normalized coverage tracks.
    
    Args:
        bam_file (str): Path to the input BAM file.
        output_dir (str): Directory to save the output bigWig file.
        logger (logging.Logger): Logger object for logging messages.
        output_filename (Optional[str]): Name for the output file. 
            If None, defaults to '{base_name}.{normalize_using}.bigwig'.
        bin_size (int): Size of the bins in bases. Defaults to 50.
        normalize_using (str): Normalization method to use. Defaults to "RPGC".
        effective_genome_size (int): The effective genome size used for normalization. 
            Defaults to 2913022398 (human genome).
        bamcoverage_path (str): Path to the bamCoverage executable. Defaults to 'bamCoverage'.
        threads (int): Number of threads to use. Defaults to 1.
        log_file (Optional[str]): Path to save the stdout and stderr from bamCoverage.
            If None, logs are not saved to a separate file. Defaults to None.
        additional_options (Optional[List[str]]): List of additional arguments
            to pass to bamCoverage. Defaults to None.
            
    Returns:
        Optional[str]: The full path to the generated bigWig file if successful,
                       None otherwise.
    """
    # --- 入力値の検証 ---
    if not bam_file:
        logger.error("No BAM file provided to run_bamcoverage.")
        return None
    
    if not os.path.exists(bam_file):
        logger.error(f"BAM file not found: {bam_file}")
        return None
        
    # --- bamCoverage のバージョン確認 ---
    _get_bamcoverage_version(bamcoverage_path, logger)
        
    # --- 出力ディレクトリの確認 ---
    try:
        os.makedirs(output_dir, exist_ok=True)
        logger.info(f"Ensured output directory exists: {output_dir}")
    except OSError as e:
        logger.error(f"Failed to create output directory {output_dir}: {e}")
        return None
        
    # --- 出力ファイル名の決定 ---
    if output_filename is None:
        base_name = os.path.splitext(os.path.basename(bam_file))[0]
        output_filename = f"{base_name}.{normalize_using}.bigwig"
        
    output_file_path = os.path.join(output_dir, output_filename)
    
    # --- bamCoverage コマンドリストの構築 ---
    cmd = [bamcoverage_path]
    
    # 基本オプション
    cmd.extend(['--bam', bam_file])
    cmd.extend(['--outFileName', output_file_path])
    cmd.extend(['--binSize', str(bin_size)])
    cmd.extend(['--numberOfProcessors', str(threads)])
    
    # 正規化オプション
    if normalize_using.upper() == "RPGC":
        cmd.extend(['--normalizeUsing', 'RPGC'])
        cmd.extend(['--effectiveGenomeSize', str(effective_genome_size)])
    else:
        cmd.extend(['--normalizeUsing', normalize_using])
    
    # 追加オプション
    if additional_options:
        cmd.extend(additional_options)
    
    # --- サブプロセスとして bamCoverage を実行 ---
    # 実行開始ログ
    logger.info(f"Running bamCoverage for {os.path.basename(bam_file)}...")
    # 実行するコマンド文字列をログに出力
    logger.debug(f"Command: {' '.join(cmd)}")
    
    stdout_log: Optional[str] = None
    stderr_log: Optional[str] = None
    
    # 実行時間計測の準備
    import time
    start_time = time.time()
    
    try:
        # コマンド実行
        result = subprocess.run(cmd, check=True, capture_output=True, text=True, encoding='utf-8')
        
        # 実行時間の計算と記録
        end_time = time.time()
        execution_time = end_time - start_time
        logger.debug(f"bamCoverage execution took {execution_time:.2f} seconds.")
        
        # 実行結果のstdout/stderrを取得
        stdout_log = result.stdout
        stderr_log = result.stderr
        
        # 正常終了ログ
        logger.info("bamCoverage finished successfully.")
        
        # stdoutがあればデバッグログ出力
        if stdout_log:
            logger.debug(f"bamCoverage stdout:\n{stdout_log}")
        
        # stderrがあればデバッグログ出力
        if stderr_log:
            logger.debug(f"bamCoverage stderr:\n{stderr_log}")
            
        # ログファイルが指定されていれば保存
        if log_file:
            _save_bamcoverage_log(log_file, stdout_log, stderr_log, logger)
            
        # 成功したので出力ファイルの絶対パスを返す
        return os.path.abspath(output_file_path)
        
    # --- エラーハンドリング ---
    except FileNotFoundError:
        # コマンドが見つからない場合
        logger.error(f"Error: '{bamcoverage_path}' command not found.")
        logger.error("Please ensure bamCoverage is installed and in your PATH, or provide the correct full path.")
        return None
    except subprocess.CalledProcessError as e:
        # コマンドがエラー終了した場合
        # 実行時間の計算と記録
        end_time = time.time()
        execution_time = end_time - start_time
        logger.debug(f"Failed bamCoverage execution took {execution_time:.2f} seconds.")
        
        stdout_log = e.stdout
        stderr_log = e.stderr
        logger.error(f"bamCoverage failed with exit code {e.returncode}")
        if stdout_log: logger.error(f"Captured stdout on failure:\n{stdout_log}")
        if stderr_log: logger.error(f"Captured stderr on failure:\n{stderr_log}")
        
        # エラー時もログ保存試行
        if log_file:
            _save_bamcoverage_log(log_file, stdout_log, stderr_log, logger)
            
        # 不完全ファイルの削除試行
        if os.path.exists(output_file_path):
            logger.warning(f"Attempting to remove potentially incomplete output file: {output_file_path}")
            try: 
                os.remove(output_file_path)
            except OSError as rm_err: 
                logger.error(f"Failed to remove incomplete output file {output_file_path}: {rm_err}")
        return None
    except Exception as e:
        # その他の予期せぬエラー
        # 実行時間の計算と記録
        end_time = time.time()
        execution_time = end_time - start_time
        logger.debug(f"Failed bamCoverage execution took {execution_time:.2f} seconds.")
        
        logger.error(f"An unexpected error occurred while running bamCoverage: {e}", exc_info=True)
        
        # ログ保存試行
        if log_file:
            _save_bamcoverage_log(log_file, stdout_log, stderr_log, logger)
        return None


def run_multibigwigsummary(
    bigwig_files: List[str],
    bed_file: str,
    output_dir: str,
    logger: logging.Logger,
    output_filename: Optional[str] = None,
    labels: Optional[List[str]] = None,
    multibigwigsummary_path: str = 'multiBigwigSummary',
    threads: int = 1,
    log_file: Optional[str] = None,
    additional_options: Optional[List[str]] = None
) -> Optional[Tuple[str, str]]:
    """
    Runs multiBigwigSummary from deeptools to generate a matrix of values.
    
    Args:
        bigwig_files (List[str]): List of paths to the input bigWig files.
        bed_file (str): Path to the BED file defining regions for analysis.
        output_dir (str): Directory to save the output files.
        logger (logging.Logger): Logger object for logging messages.
        output_filename (Optional[str]): Base name for the output files (without extension).
            If None, defaults to "bigwig_summary".
        labels (Optional[List[str]]): Labels for the bigWig files.
            If None, uses --smartLabels option to automatically determine labels.
            If provided, these exact labels will be used with --labels option.
        multibigwigsummary_path (str): Path to the multiBigwigSummary executable.
            Defaults to 'multiBigwigSummary'.
        threads (int): Number of threads to use. Defaults to 1.
        log_file (Optional[str]): Path to save the stdout and stderr.
            If None, logs are not saved to a separate file. Defaults to None.
        additional_options (Optional[List[str]]): List of additional arguments
            to pass to multiBigwigSummary. Defaults to None.
            
    Returns:
        Optional[Tuple[str, str]]: A tuple containing paths to the generated npz and raw counts
                                  files if successful, None otherwise.
    """
    # --- 入力値の検証 ---
    if not bigwig_files:
        logger.error("No bigWig files provided to run_multibigwigsummary.")
        return None
    
    for bw_file in bigwig_files:
        if not os.path.exists(bw_file):
            logger.error(f"bigWig file not found: {bw_file}")
            return None
    
    if not bed_file:
        logger.error("No BED file provided to run_multibigwigsummary.")
        return None
        
    if not os.path.exists(bed_file):
        logger.error(f"BED file not found: {bed_file}")
        return None
    
    # --- multiBigwigSummary のバージョン確認 ---
    _get_multibigwigsummary_version(multibigwigsummary_path, logger)
        
    # --- 出力ディレクトリの確認 ---
    try:
        os.makedirs(output_dir, exist_ok=True)
        logger.info(f"Ensured output directory exists: {output_dir}")
    except OSError as e:
        logger.error(f"Failed to create output directory {output_dir}: {e}")
        return None
        
    # --- 出力ファイル名の決定 ---
    if output_filename is None:
        output_filename = "bigwig_summary"
        
    output_npz_path = os.path.join(output_dir, f"{output_filename}.npz")
    output_raw_path = os.path.join(output_dir, f"{output_filename}.tab")
    
    # --- multiBigwigSummary コマンドリストの構築 ---
    cmd = [multibigwigsummary_path, 'BED-file']
    
    # 基本オプション
    cmd.extend(['--BED', bed_file])
    cmd.extend(['--bwfiles'] + bigwig_files)
    
    # ラベル処理の修正：labelsがNoneの場合は--smartLabelsを使用
    if labels is None:
        # ラベル指定なしの場合は--smartLabelsを使用
        logger.info("No labels specified. Using --smartLabels to automatically determine labels from filenames.")
        cmd.append('--smartLabels')
    else:
        # ラベルが指定された場合は従来通り--labelsオプションを使用
        if len(labels) != len(bigwig_files):
            logger.warning(f"Number of labels ({len(labels)}) does not match number of bigWig files ({len(bigwig_files)}). Using --smartLabels instead.")
            cmd.append('--smartLabels')
        else:
            logger.info(f"Using specified labels: {labels}")
            cmd.extend(['--labels'] + labels)
            
    cmd.extend(['--numberOfProcessors', str(threads)])
    cmd.extend(['--outFileName', output_npz_path])
    cmd.extend(['--outRawCounts', output_raw_path])
    cmd.append('-v')  # verbose モードを有効化
    
    # 追加オプション
    if additional_options:
        cmd.extend(additional_options)
    
    # --- サブプロセスとして multiBigwigSummary を実行 ---
    # 実行開始ログ
    logger.info(f"Running multiBigwigSummary for {len(bigwig_files)} bigWig files...")
    # 実行するコマンド文字列をログに出力
    logger.debug(f"Command: {' '.join(cmd)}")
    
    stdout_log: Optional[str] = None
    stderr_log: Optional[str] = None
    
    # 実行時間計測の準備
    import time
    start_time = time.time()
    
    try:
        # コマンド実行
        result = subprocess.run(cmd, check=True, capture_output=True, text=True, encoding='utf-8')
        
        # 実行時間の計算と記録
        end_time = time.time()
        execution_time = end_time - start_time
        logger.debug(f"multiBigwigSummary execution took {execution_time:.2f} seconds.")
        
        # 実行結果のstdout/stderrを取得
        stdout_log = result.stdout
        stderr_log = result.stderr
        
        # 正常終了ログ
        logger.info("multiBigwigSummary finished successfully.")
        
        # stdoutがあればデバッグログ出力
        if stdout_log:
            logger.debug(f"multiBigwigSummary stdout:\n{stdout_log}")
        
        # stderrがあればデバッグログ出力
        if stderr_log:
            logger.debug(f"multiBigwigSummary stderr:\n{stderr_log}")
            
        # ログファイルが指定されていれば保存
        if log_file:
            _save_multibigwigsummary_log(log_file, stdout_log, stderr_log, logger)
            
        # 成功したら出力ファイルのパスをタプルで返す
        return (os.path.abspath(output_npz_path), os.path.abspath(output_raw_path))
        
    # --- エラーハンドリング ---
    except FileNotFoundError:
        # コマンドが見つからない場合
        logger.error(f"Error: '{multibigwigsummary_path}' command not found.")
        logger.error("Please ensure multiBigwigSummary is installed and in your PATH, or provide the correct full path.")
        return None
    except subprocess.CalledProcessError as e:
        # コマンドがエラー終了した場合
        # 実行時間の計算と記録
        end_time = time.time()
        execution_time = end_time - start_time
        logger.debug(f"Failed multiBigwigSummary execution took {execution_time:.2f} seconds.")
        
        stdout_log = e.stdout
        stderr_log = e.stderr
        logger.error(f"multiBigwigSummary failed with exit code {e.returncode}")
        if stdout_log: logger.error(f"Captured stdout on failure:\n{stdout_log}")
        if stderr_log: logger.error(f"Captured stderr on failure:\n{stderr_log}")
        
        # エラー時もログ保存試行
        if log_file:
            _save_multibigwigsummary_log(log_file, stdout_log, stderr_log, logger)
            
        # 不完全ファイルの削除試行
        for output_file in [output_npz_path, output_raw_path]:
            if os.path.exists(output_file):
                logger.warning(f"Attempting to remove potentially incomplete output file: {output_file}")
                try: 
                    os.remove(output_file)
                except OSError as rm_err: 
                    logger.error(f"Failed to remove incomplete output file {output_file}: {rm_err}")
        return None
    except Exception as e:
        # その他の予期せぬエラー
        # 実行時間の計算と記録
        end_time = time.time()
        execution_time = end_time - start_time
        logger.debug(f"Failed multiBigwigSummary execution took {execution_time:.2f} seconds.")
        
        logger.error(f"An unexpected error occurred while running multiBigwigSummary: {e}", exc_info=True)
        
        # ログ保存試行
        if log_file:
            _save_multibigwigsummary_log(log_file, stdout_log, stderr_log, logger)
        return None


def _save_multibigwigsummary_log(log_path: str, stdout_content: Optional[str], stderr_content: Optional[str], logger: logging.Logger):
    """Saves stdout and stderr content from multiBigwigSummary to the specified log file."""
    # ログ保存試行のログ
    logger.info(f"Attempting to save multiBigwigSummary log to: {log_path}")
    try:
        # 保存先ディレクトリが存在しない場合は作成
        log_dir = os.path.dirname(log_path)
        if log_dir and not os.path.exists(log_dir):
            os.makedirs(log_dir, exist_ok=True)
            # ディレクトリ作成ログ
            logger.info(f"Created directory for log file: {log_dir}")

        # ファイルに書き込み (encoding指定)
        with open(log_path, 'w', encoding='utf-8') as f:
            f.write("--- multiBigwigSummary stdout ---\n")
            # None や空文字列の場合に対応
            f.write((stdout_content if stdout_content else "[No stdout captured or content was empty]") + "\n")
            f.write("\n--- multiBigwigSummary stderr ---\n")
            # None や空文字列の場合に対応
            f.write((stderr_content if stderr_content else "[No stderr captured or content was empty]") + "\n")
        # ログ保存成功ログ
        logger.info(f"Successfully saved multiBigwigSummary log to: {log_path}")
    except Exception as e:
        # ログ保存失敗時のエラーログ
        logger.error(f"Failed to save multiBigwigSummary log to {log_path}: {e}")


def _get_multibigwigsummary_version(multibigwigsummary_path: str, logger: logging.Logger) -> Optional[str]:
    """Gets the version of multiBigwigSummary by running 'multiBigwigSummary --version'."""
    try:
        # バージョン確認コマンド実行
        version_cmd = [multibigwigsummary_path, "--version"]
        result = subprocess.run(version_cmd, check=True, capture_output=True, text=True, encoding='utf-8')
        version_info = result.stdout.strip() if result.stdout else "Unknown version"
        logger.debug(f"multiBigwigSummary version: {version_info}")
        return version_info
    except Exception as e:
        logger.warning(f"Failed to get multiBigwigSummary version: {e}")
        return None




##### multibigwig_tab操作関連 #####

def read_multibigwig_counts(
    tab_file: str,
    logger: logging.Logger
) -> Optional[pd.DataFrame]:
    """
    Reads a multiBigwigSummary tab file and returns a pandas DataFrame.
    Handles special format with '#' prefixed headers and quoted column names.
    
    Args:
        tab_file (str): Path to the tab-delimited raw counts file from multiBigwigSummary.
        logger (logging.Logger): Logger object for logging messages.
            
    Returns:
        Optional[pd.DataFrame]: A DataFrame containing the raw counts,
                               or None if an error occurs.
    """
    logger.info(f"Reading multiBigwigSummary tab file: {tab_file}")
    
    try:
        # ファイルの存在確認
        if not os.path.exists(tab_file):
            logger.error(f"Tab file not found: {tab_file}")
            return None
        
        # まず最初の数行を読んで形式を確認
        try:
            with open(tab_file, 'r') as f:
                first_lines = [next(f) for _ in range(5) if f]
            logger.debug(f"File format check - first few lines:\n{''.join(first_lines[:3])}")
        except Exception as e:
            logger.warning(f"Could not preview file format: {e}")
            
        # tabファイルを読み込み - 特殊な形式に対応
        try:
            # まず標準的な方法で読み込みを試みる
            df = pd.read_csv(tab_file, sep='\t', comment=None)
            
            # 特殊な形式の場合（#で始まるヘッダーと引用符で囲まれた列名）を検出して処理
            if df.columns[0].startswith('#'):
                logger.debug("Detected header starting with '#', processing special format")
                
                # ヘッダー行を直接読み込み、処理
                with open(tab_file, 'r') as f:
                    header_line = f.readline().strip()
                
                # '#'を削除し、引用符も削除
                if header_line.startswith('#'):
                    header_line = header_line[1:]
                    # 引用符を取り除く
                    column_names = [col.strip("'\"") for col in header_line.split('\t')]
                    logger.debug(f"Processed column names: {column_names}")
                    
                    # ファイルを再度読み込み、処理したヘッダーを使用
                    df = pd.read_csv(tab_file, sep='\t', comment=None, names=column_names, skiprows=1)
            
            logger.info(f"Successfully read tab file with {len(df)} rows and {len(df.columns)} columns.")
            logger.debug(f"Column names after processing: {df.columns.tolist()}")
            
        except pd.errors.EmptyDataError:
            logger.error(f"Tab file is empty or contains no data: {tab_file}")
            return None
        except Exception as e:
            logger.error(f"Failed to read tab file {tab_file}: {e}")
            return None
            
        # カラム名チェック
        expected_columns = ['chr', 'start', 'end']
        missing_columns = [col for col in expected_columns if col not in df.columns]
        if missing_columns:
            logger.error(f"Required columns missing from tab file: {missing_columns}")
            logger.error(f"Available columns: {df.columns.tolist()}")
            return None
            
        # データサンプルのログ出力（デバッグ用）
        if not df.empty:
            logger.debug(f"First few rows of data:\n{df.head().to_string()}")
            
        # 重複行のチェック
        duplicate_regions = df.duplicated(subset=['chr', 'start', 'end'], keep=False)
        if duplicate_regions.any():
            dup_count = duplicate_regions.sum()
            logger.warning(f"Found {dup_count} duplicate genomic regions in the tab file.")
        
        return df
        
    except Exception as e:
        logger.error(f"An unexpected error occurred while reading tab file: {e}", exc_info=True)
        return None

# ====== 1. DataFrame を PyRanges に変換する関数 ======

def dataframe_to_pyranges(
    df: pd.DataFrame,
    logger: logging.Logger,
    chr_col: str = 'chr',
    start_col: str = 'start',
    end_col: str = 'end'
) -> Optional[pr.PyRanges]:
    """
    Convert DataFrame to PyRanges object
    
    Args:
        df (pd.DataFrame): DataFrame containing genomic region information
        logger (logging.Logger): Logger object
        chr_col (str): Chromosome column name. Default: 'chr'
        start_col (str): Start position column name. Default: 'start'
        end_col (str): End position column name. Default: 'end'
    
    Returns:
        Optional[pr.PyRanges]: Converted PyRanges object, or None if error occurs
    """
    logger.info("Converting DataFrame to PyRanges object")
    
    try:
        # 必要なカラムが存在するか確認
        for col in [chr_col, start_col, end_col]:
            if col not in df.columns:
                logger.error(f"Required column '{col}' not found in DataFrame.")
                return None
        
        # 数値型の確認と変換
        for col in [start_col, end_col]:
            if not pd.api.types.is_numeric_dtype(df[col]):
                logger.warning(f"Column '{col}' is not numeric. Attempting to convert.")
                try:
                    df = df.copy()  # コピーを作成して元のデータフレームを変更しない
                    df[col] = pd.to_numeric(df[col], errors='coerce')
                    
                    # 変換後に NaN がないか確認
                    if df[col].isnull().any():
                        logger.error(f"Column '{col}' contains values that could not be converted to numeric.")
                        return None
                except Exception as e:
                    logger.error(f"Error converting column '{col}' to numeric: {e}")
                    return None
        
        # コピーを作成して元のデータフレームを変更しない
        df_copy = df.copy()
        
        # PyRanges が期待するカラム名にリネーム
        rename_dict = {
            chr_col: 'Chromosome',
            start_col: 'Start',
            end_col: 'End'
        }
        df_copy = df_copy.rename(columns=rename_dict)
        
        # PyRanges オブジェクトを作成
        try:
            gr = pr.PyRanges(df_copy)
            logger.debug(f"Successfully created PyRanges object with {len(gr)} regions")
            return gr
        except Exception as e:
            logger.error(f"Failed to create PyRanges object: {e}")
            return None
            
    except ImportError:
        logger.error("PyRanges library not found. Please install it using 'pip install pyranges'.")
        return None
    except Exception as e:
        logger.error(f"An unexpected error occurred during conversion to PyRanges: {e}", exc_info=True)
        return None


# ====== 2. PyRanges オブジェクトをソートする関数 ======

def sort_pyranges(
    gr: pr.PyRanges,
    logger: logging.Logger
) -> Optional[pr.PyRanges]:
    """
    Sort PyRanges object by chromosome and coordinates
    
    Args:
        gr (pr.PyRanges): PyRanges object to sort
        logger (logging.Logger): Logger object
    
    Returns:
        Optional[pr.PyRanges]: Sorted PyRanges object, or None if error occurs
    """
    logger.info("Sorting PyRanges object")
    
    try:
        # None チェック
        if gr is None:
            logger.error("Input PyRanges object is None")
            return None
            
        # ソート実行
        sorted_gr = gr.sort()
        logger.debug("Successfully sorted PyRanges object")
        return sorted_gr
        
    except Exception as e:
        logger.error(f"Error during PyRanges sorting: {e}", exc_info=True)
        return None


# ====== 3. PyRanges を DataFrame に戻す関数 ======

def pyranges_to_dataframe(
    gr: pr.PyRanges,
    logger: logging.Logger,
    original_chr_col: str = 'chr',
    original_start_col: str = 'start',
    original_end_col: str = 'end'
) -> Optional[pd.DataFrame]:
    """
    Convert PyRanges object to DataFrame
    
    Args:
        gr (pr.PyRanges): PyRanges object to convert
        logger (logging.Logger): Logger object
        original_chr_col (str): Original chromosome column name. Default: 'chr'
        original_start_col (str): Original start position column name. Default: 'start'
        original_end_col (str): Original end position column name. Default: 'end'
    
    Returns:
        Optional[pd.DataFrame]: Converted DataFrame, or None if error occurs
    """
    logger.info("Converting PyRanges object back to DataFrame")
    
    try:
        # None チェック
        if gr is None:
            logger.error("Input PyRanges object is None")
            return None
            
        # PyRanges を DataFrame に変換
        df = gr.df
        
        # 元のカラム名に戻す
        rename_dict = {
            'Chromosome': original_chr_col,
            'Start': original_start_col,
            'End': original_end_col
        }
        df = df.rename(columns=rename_dict)
        
        # インデックスをリセット
        df = df.reset_index(drop=True)
        
        logger.debug(f"Successfully converted PyRanges back to DataFrame with {len(df)} rows")
        return df
        
    except Exception as e:
        logger.error(f"Error converting PyRanges to DataFrame: {e}", exc_info=True)
        return None



##########################




def apply_log2_transform_to_counts(
    df: pd.DataFrame,
    logger: logging.Logger,
    value_columns: Optional[List[str]] = None,
    exclude_patterns: List[str] = [],
    pseudocount: float = 1.0
) -> Optional[pd.DataFrame]:
    """
    Applies log2(x + pseudocount) transformation to value columns and returns a new DataFrame.
    
    Args:
        df (pd.DataFrame): Input DataFrame with genomic regions and values.
        logger (logging.Logger): Logger object for logging messages.
        value_columns (Optional[List[str]]): List of column names to transform.
            If None, all numeric columns except 'chr', 'start', 'end' will be transformed.
        exclude_patterns (List[str]): List of patterns to exclude columns.
            Default is an empty list (no exclusions).
        pseudocount (float): Pseudocount to add before log transformation. Defaults to 1.0.
            
    Returns:
        Optional[pd.DataFrame]: A new DataFrame with log2-transformed values, or None if an error occurs.
    """
    logger.info("Applying log2 transformation to count values...")
    
    try:
        # DataFrameの検証
        if not isinstance(df, pd.DataFrame) or df.empty:
            logger.error("Input is not a valid or non-empty DataFrame.")
            return None
            
        # 必須列の確認
        for col in ['chr', 'start', 'end']:
            if col not in df.columns:
                logger.error(f"Required column '{col}' not found in DataFrame.")
                return None
                
        # 新しいDataFrameを作成（入力をコピー）
        result_df = df.copy()
        
        # 処理対象の列を決定
        if value_columns is None:
            # 'chr', 'start', 'end'以外の数値型の列を抽出
            numeric_cols = result_df.select_dtypes(include=['number']).columns.tolist()
            value_columns = [col for col in numeric_cols if col not in ['chr', 'start', 'end']]
            
            # 除外パターンに一致する列を除外
            if exclude_patterns:
                value_columns = [col for col in value_columns if not any(pattern in col for pattern in exclude_patterns)]
                
        logger.info(f"Columns selected for log2 transformation: {value_columns}")
        
        if not value_columns:
            logger.warning("No columns found for log2 transformation.")
            return result_df
            
        # 変換前の値の範囲を記録 (デバッグ用)
        for col in value_columns:
            min_val = result_df[col].min()
            max_val = result_df[col].max()
            logger.debug(f"Column '{col}' value range before transformation: [{min_val}, {max_val}]")
            
            # 負の値があるか確認
            if min_val < 0:
                logger.warning(f"Column '{col}' contains negative values. Consider how pseudocount affects interpretation.")
        
        # 各列にlog2(x + pseudocount)変換を適用
        transform_count = 0
        for col in value_columns:
            try:
                # log2変換を適用して元の列を上書き
                result_df[col] = np.log2(result_df[col] + pseudocount)
                
                # 変換後の値の範囲を記録
                min_log_val = result_df[col].min()
                max_log_val = result_df[col].max()
                logger.debug(f"Column '{col}' value range after transformation: [{min_log_val}, {max_log_val}]")
                
                transform_count += 1
                
            except Exception as e:
                logger.error(f"Error transforming column '{col}': {e}")
                continue
                
        logger.info(f"Successfully applied log2 transformation to {transform_count} columns.")
        return result_df
        
    except Exception as e:
        logger.error(f"An unexpected error occurred during log2 transformation: {e}", exc_info=True)
        return None

def read_bed_to_dataframe(
    bed_file: str,
    logger: logging.Logger
) -> Optional[pd.DataFrame]:
    """
    Read BED file concisely using PyRanges and convert to pandas DataFrame.
    Assumes PyRanges is installed.
    After loading, explicitly sort by genomic coordinates.
    
    Args:
        bed_file (str): Path to the BED file
        logger (logging.Logger): Logger object
            
    Returns:
        Optional[pd.DataFrame]: DataFrame containing sorted BED regions, or None if error occurs
    """
    logger.info(f"Reading BED file using PyRanges: {bed_file}")
    
    try:
        # BEDファイルの存在確認
        if not os.path.exists(bed_file):
            logger.error(f"BED file not found: {bed_file}")
            return None
        
        # PyRangesでBEDファイルを読み込み
        gr = pr.read_bed(bed_file)
        logger.info(f"Successfully loaded BED file into PyRanges with {len(gr)} regions")
        
        # 明示的にPyRangesオブジェクトをソート
        logger.debug("Explicitly sorting genomic regions in PyRanges object")
        sorted_gr = gr.sort()
        logger.debug("Sorting completed")
        
        # ソート前後の比較情報をログに出力（デバッグ用）
        if len(gr) > 0:
            logger.debug(f"First region before sort - Chromosome: {gr.df['Chromosome'].iloc[0]}, Start: {gr.df['Start'].iloc[0]}")
            logger.debug(f"First region after sort - Chromosome: {sorted_gr.df['Chromosome'].iloc[0]}, Start: {sorted_gr.df['Start'].iloc[0]}")
        
        # ソート済みPyRangesオブジェクトをDataFrameに変換
        bed_df = sorted_gr.df.copy()
        
        # 列名を標準化 (PyRangesの'Chromosome', 'Start', 'End'を'chr', 'start', 'end'に変換)
        rename_dict = {
            'Chromosome': 'chr',
            'Start': 'start',
            'End': 'end'
        }
        
        # その他のBED標準カラムがあれば追加
        standard_bed_columns = {
            'Name': 'name',
            'Score': 'score',
            'Strand': 'strand',
            'ThickStart': 'thickStart',
            'ThickEnd': 'thickEnd',
            'ItemRGB': 'itemRgb',
            'BlockCount': 'blockCount',
            'BlockSizes': 'blockSizes',
            'BlockStarts': 'blockStarts'
        }
        
        # 実際にDataFrameに存在するカラムのみリネーム
        for orig, new in standard_bed_columns.items():
            if orig in bed_df.columns:
                rename_dict[orig] = new
        
        # リネーム実行
        bed_df = bed_df.rename(columns=rename_dict)
        
        # インデックスをリセット
        bed_df = bed_df.reset_index(drop=True)
        
        logger.info(f"Read and sorted BED file with {len(bed_df)} regions and {len(bed_df.columns)} columns")
        return bed_df
            
    except Exception as e:
        logger.error(f"Error reading BED file with PyRanges: {e}", exc_info=True)
        return None
    
def verify_regions_with_bed_df(
    df: pd.DataFrame,
    bed_df: pd.DataFrame,
    logger: logging.Logger,
    chr_col: str = 'chr',
    start_col: str = 'start',
    end_col: str = 'end'
) -> bool:
    """
    Verifies that the genomic regions in the DataFrame match those in the BED DataFrame.
    
    Args:
        df (pd.DataFrame): DataFrame with genomic regions (chr, start, end columns).
        bed_df (pd.DataFrame): BED DataFrame with columns ['chr', 'start', 'end'].
        logger (logging.Logger): Logger object for logging messages.
        chr_col (str): Name of the chromosome column in the DataFrame. Defaults to 'chr'.
        start_col (str): Name of the start position column in the DataFrame. Defaults to 'start'.
        end_col (str): Name of the end position column in the DataFrame. Defaults to 'end'.
            
    Returns:
        bool: True if regions match, False otherwise.
    """
    logger.info("Verifying that DataFrame regions match BED DataFrame")
    
    try:
        # 必要なカラムがDataFrameに存在するか確認
        for col in [chr_col, start_col, end_col]:
            if col not in df.columns:
                logger.error(f"Required column '{col}' not found in DataFrame.")
                return False
        
        # 必要なカラムがBED DataFrameに存在するか確認
        for col in ['chr', 'start', 'end']:
            if col not in bed_df.columns:
                logger.error(f"Required column '{col}' not found in BED DataFrame.")
                return False
                
        # DataFrameの領域情報を抽出
        df_regions = df[[chr_col, start_col, end_col]].copy()
        # 型変換を確保（特にstart/endは整数であるべき）
        df_regions[start_col] = df_regions[start_col].astype(int)
        df_regions[end_col] = df_regions[end_col].astype(int)
        
        # BED DataFrameの必要なカラムだけを取得
        bed_regions = bed_df[['chr', 'start', 'end']].copy()
        bed_regions['start'] = bed_regions['start'].astype(int)
        bed_regions['end'] = bed_regions['end'].astype(int)
            
        # 行数チェック
        if len(df_regions) != len(bed_regions):
            logger.error(f"Region count mismatch: DataFrame has {len(df_regions)} regions, BED DataFrame has {len(bed_regions)} regions.")
            return False
            
        # 両方ともソート（同じ順序で比較するため）
        df_regions = df_regions.sort_values(by=[chr_col, start_col, end_col]).reset_index(drop=True)
        bed_regions = bed_regions.sort_values(by=['chr', 'start', 'end']).reset_index(drop=True)
        
        # 各カラムの一致を確認
        chr_match = df_regions[chr_col].equals(bed_regions['chr'])
        start_match = df_regions[start_col].equals(bed_regions['start'])
        end_match = df_regions[end_col].equals(bed_regions['end'])
        
        if not chr_match:
            logger.error("Chromosome names do not match between DataFrame and BED DataFrame.")
            # 不一致の例を表示
            mismatch_idx = (df_regions[chr_col] != bed_regions['chr']).idxmax()
            logger.error(f"First mismatch at index {mismatch_idx}: DataFrame={df_regions.loc[mismatch_idx, chr_col]}, BED={bed_regions.loc[mismatch_idx, 'chr']}")
            return False
            
        if not start_match:
            logger.error("Start positions do not match between DataFrame and BED DataFrame.")
            # 不一致の例を表示
            mismatch_idx = (df_regions[start_col] != bed_regions['start']).idxmax()
            logger.error(f"First mismatch at index {mismatch_idx}: DataFrame={df_regions.loc[mismatch_idx, start_col]}, BED={bed_regions.loc[mismatch_idx, 'start']}")
            return False
            
        if not end_match:
            logger.error("End positions do not match between DataFrame and BED DataFrame.")
            # 不一致の例を表示
            mismatch_idx = (df_regions[end_col] != bed_regions['end']).idxmax()
            logger.error(f"First mismatch at index {mismatch_idx}: DataFrame={df_regions.loc[mismatch_idx, end_col]}, BED={bed_regions.loc[mismatch_idx, 'end']}")
            return False
            
        # すべてのカラムが一致した場合
        logger.info("DataFrame regions match exactly with BED DataFrame regions.")
        return True
        
    except Exception as e:
        logger.error(f"An unexpected error occurred during region verification: {e}", exc_info=True)
        return False
    






def add_peak_ids_from_bed(
    multibigwig_counts_df: pd.DataFrame,
    bed_df: pd.DataFrame,
    logger: logging.Logger,
    force_chr_start_end_ids: bool = False
) -> Optional[pd.DataFrame]:
    """
    Add PeakID information from bed_df to multibigwig_counts_df.
    After verifying region match with verify_regions_with_bed_df,
    adds the 4th column (name) information from bed_df as PeakID and standardizes column names.
    
    Args:
        multibigwig_counts_df (pd.DataFrame): DataFrame with multiBigwig counts (chr, start, end, sample columns)
        bed_df (pd.DataFrame): DataFrame loaded from BED file (at least containing chr, start, end columns)
        logger (logging.Logger): Logger object
        force_chr_start_end_ids (bool): If True, always use chr_start_end as ID regardless of name column in bed_df
        
    Returns:
        Optional[pd.DataFrame]: DataFrame with added PeakID column and renamed columns. None if error occurs
    """
    logger.info("Adding PeakID information from BED file")
    
    try:
        # DataFrameが空でないか確認
        if multibigwig_counts_df.empty or bed_df.empty:
            logger.error("One or both input DataFrames are empty")
            return None
            
        # 両方のDataFrameに必要なカラムが存在するか確認
        required_columns = ['chr', 'start', 'end']
        for df_name, df in [("multibigwig_counts_df", multibigwig_counts_df), ("bed_df", bed_df)]:
            missing_cols = [col for col in required_columns if col not in df.columns]
            if missing_cols:
                logger.error(f"Required columns missing in {df_name}: {missing_cols}")
                return None
        
        # verify_regions_with_bed_dfを使って領域の一致を確認
        logger.info("Verifying genomic regions match between multibigwig_counts_df and bed_df")
        regions_match = verify_regions_with_bed_df(
            multibigwig_counts_df, bed_df, logger, 
            chr_col='chr', start_col='start', end_col='end'
        )
        
        if not regions_match:
            logger.error("Genomic regions in multibigwig_counts_df do not match those in bed_df")
            return None
        
        # 一致を確認したら新しいデータフレームを作成 (元のデータを変更しないため)
        result_df = multibigwig_counts_df.copy()
        
        # PeakID列の内容を決定
        has_name_col = 'name' in bed_df.columns and not bed_df['name'].isnull().any()
        
        if has_name_col and not force_chr_start_end_ids:
            # bed_dfの4列目 (name) を使用
            logger.info("Using 'name' column from BED file as PeakID")
            peak_ids = bed_df['name'].tolist()
        else:
            # chr_start_end形式のIDを生成
            logger.info("Generating PeakID as 'chr_start_end'")
            if force_chr_start_end_ids:
                logger.debug("Using chr_start_end format as requested by force_chr_start_end_ids=True")
            elif not has_name_col:
                logger.debug("No 'name' column found in BED file, using chr_start_end format")
                
            # chr, start, endを結合してIDを生成
            peak_ids = [f"{c}_{s}_{e}" for c, s, e in zip(
                bed_df['chr'], bed_df['start'], bed_df['end']
            )]
        
        # 列名のリネーム用辞書を作成
        rename_dict = {
            'chr': 'Chr',
            'start': 'Start',
            'end': 'End'
        }
        
        # 列名変更
        result_df = result_df.rename(columns=rename_dict)
        
        # PeakID列をデータフレームの先頭に追加
        result_df.insert(0, 'PeakID', peak_ids)
        
        # 処理結果のサマリー
        logger.info(f"Successfully added PeakID column. Result DataFrame has {len(result_df)} rows and {len(result_df.columns)} columns")
        logger.debug(f"Final column order: {list(result_df.columns)}")
        
        return result_df
        
    except Exception as e:
        logger.error(f"Error while adding PeakID information: {e}", exc_info=True)
        return None
    

def convert_saf_to_bed(saf_path: str, logger: logging.Logger, output_dir: Optional[str] = None) -> Optional[str]:
    """
    Convert SAF format file to BED file
    
    Args:
        saf_path (str): Path to the SAF file
        logger (logging.Logger): Logger object
        output_dir (Optional[str]): Directory to save the converted BED file.
                                    If specified, the file is saved as a log.
        
    Returns:
        Optional[str]: Path to the generated BED file, or None if error occurs
    """
    saf_basename = os.path.basename(saf_path)
    saf_name = os.path.splitext(saf_basename)[0]
    logger.info(f"Converting SAF file to BED: {saf_basename}")
    
    try:
        # SAFファイルを読み込み（ヘッダーあり、タブ区切り）
        df = pd.read_csv(saf_path, sep='\t', header=0)
        
        # 必要なカラムが存在するか確認
        required_cols = ['GeneID', 'Chr', 'Start', 'End', 'Strand']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            logger.error(f"Required columns missing in SAF file: {missing_cols}")
            return None
        
        # BEDデータフレーム用のデータ構築
        bed_data = {
            'Chr': df['Chr'],
            'Start': df['Start'] - 1,  # SAFは1-based、BEDは0-based
            'End': df['End'],
            'Name': df['GeneID'],
            'Score': 0,  # デフォルトスコア
            'Strand': df['Strand'].replace('.', '*')  # 不明なストランドは '*' に変換
        }
        
        bed_df = pd.DataFrame(bed_data)
        
        # ファイルパスを決定
        if output_dir:
            # 指定されたディレクトリにログとして保存
            os.makedirs(output_dir, exist_ok=True)
            conversion_dir = os.path.join(output_dir, "conversion_logs")
            os.makedirs(conversion_dir, exist_ok=True)
            output_bed_path = os.path.join(conversion_dir, f"{saf_name}_converted.bed")
            
            # ファイルに書き込み
            bed_df.to_csv(output_bed_path, sep='\t', index=False, header=False)
            logger.info(f"Saved converted BED file to: {output_bed_path}")
            return output_bed_path
        else:
            # 一時ファイルに書き込み
            temp_bed_file = tempfile.NamedTemporaryFile(mode='w+', suffix='.bed', delete=False, encoding='utf-8')
            bed_df.to_csv(temp_bed_file.name, sep='\t', index=False, header=False)
            temp_bed_path = temp_bed_file.name
            temp_bed_file.close()
            logger.info(f"Created temporary BED file: {temp_bed_path}")
            return temp_bed_path
    
    except Exception as e:
        logger.error(f"Error converting SAF to BED: {e}", exc_info=True)
        return None


def convert_mergese_to_bed(merge_se_path: str, logger: logging.Logger, output_dir: Optional[str] = None) -> Optional[str]:
    """
    Convert merge_SE.tsv format file to BED file
    
    Args:
        merge_se_path (str): Path to the merge_SE.tsv file
        logger (logging.Logger): Logger object
        output_dir (Optional[str]): Directory to save the converted BED file.
                                    If specified, the file is saved as a log.
        
    Returns:
        Optional[str]: Path to the generated BED file, or None if error occurs
    """
    mergese_basename = os.path.basename(merge_se_path)
    mergese_name = os.path.splitext(mergese_basename)[0]
    logger.info(f"Converting merge_SE format file to BED: {mergese_basename}")
    
    try:
        # ファイル読み込み（ヘッダーあり、タブ区切り）
        df = pd.read_csv(merge_se_path, sep='\t', header=0)
        
        # 最初の列名を取得
        se_data_col_name = df.columns[0]
        if se_data_col_name.lower() != 'se_data':
            logger.error(f"Expected first column to be 'se_data' (case-insensitive), but found '{se_data_col_name}'.")
            return None
        
        # se_data列から座標を抽出（chr_start_end形式を想定）
        coords = df[se_data_col_name].str.split('_', expand=True, n=2)
        if coords.shape[1] < 3:
            logger.error(f"Could not split '{se_data_col_name}' column into at least 3 parts (chr, start, end).")
            return None
        
        # BEDデータフレーム用のデータ構築
        bed_data = {
            'Chr': coords[0],
            'Start': pd.to_numeric(coords[1], errors='coerce'),  # すでに0-based
            'End': pd.to_numeric(coords[2], errors='coerce'),
            'Name': df[se_data_col_name],  # 元のIDをそのまま使用
            'Score': 0,  # デフォルトスコア
            'Strand': '*'  # ストランド情報がないので '*'
        }
        
        # 数値変換チェック
        if bed_data['Start'].isnull().any() or bed_data['End'].isnull().any():
            num_failed = bed_data['Start'].isnull().sum() + bed_data['End'].isnull().sum()
            logger.error(f"Could not convert Start/End to numeric values for {num_failed} entries.")
            return None
        
        bed_df = pd.DataFrame(bed_data)
        
        # ファイルパスを決定
        if output_dir:
            # 指定されたディレクトリにログとして保存
            os.makedirs(output_dir, exist_ok=True)
            conversion_dir = os.path.join(output_dir, "conversion_logs")
            os.makedirs(conversion_dir, exist_ok=True)
            output_bed_path = os.path.join(conversion_dir, f"{mergese_name}_converted.bed")
            
            # ファイルに書き込み
            bed_df.to_csv(output_bed_path, sep='\t', index=False, header=False)
            logger.info(f"Saved converted BED file to: {output_bed_path}")
            return output_bed_path
        else:
            # 一時ファイルに書き込み
            temp_bed_file = tempfile.NamedTemporaryFile(mode='w+', suffix='.bed', delete=False, encoding='utf-8')
            bed_df.to_csv(temp_bed_file.name, sep='\t', index=False, header=False)
            temp_bed_path = temp_bed_file.name
            temp_bed_file.close()
            logger.info(f"Created temporary BED file: {temp_bed_path}")
            return temp_bed_path
    
    except Exception as e:
        logger.error(f"Error converting merge_SE to BED: {e}", exc_info=True)
        return None
    

# jupyter code

def get_jupyter_logger(name: str = "JupyterLogger", log_level: str = "INFO") -> logging.Logger:
    """
    Creates a simple logger for use in Jupyter notebooks.
    
    Args:
        name (str): Name for the logger. Defaults to "JupyterLogger".
        log_level (str): Logging level (DEBUG/INFO/WARNING/ERROR/CRITICAL). Defaults to "INFO".
            
    Returns:
        logging.Logger: Configured logger object.
    """
    # ログレベルの文字列を実際のログレベルに変換
    log_level_str = log_level.upper()
    log_level_val = getattr(logging, log_level_str, logging.INFO)
    
    # ロガーの取得と設定
    logger = logging.getLogger(name)
    logger.setLevel(log_level_val)
    
    # 既存のハンドラを全て削除（重複を避けるため）
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)
        handler.close()
    
    # コンソール出力用のハンドラ設定
    console_handler = logging.StreamHandler()
    console_handler.setLevel(log_level_val)
    
    # フォーマッタ設定
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s', 
                                 datefmt='%Y-%m-%d %H:%M:%S')
    console_handler.setFormatter(formatter)
    
    # ハンドラをロガーに追加
    logger.addHandler(console_handler)
    
    return logger
