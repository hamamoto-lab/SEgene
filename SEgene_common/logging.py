# SEgene_common/logging.py
import logging
import os
import sys
import logging.handlers
from logging.handlers import RotatingFileHandler
from typing import Optional

class LoggerManager:
    """
    Unified logger manager class for SEdb tools.
    """

    def __init__(
        self, 
        logger_name: str, 
        enable_console_logging: bool = True, 
        console_verbose: bool = False,
        enable_file_logging: bool = False,
        log_dir: str = './logs',
        log_file_name: Optional[str] = None,
        log_level: int = logging.INFO
    ):
        """
        Initializes the logger manager.
        
        Parameters:
        -----------
        logger_name : str
            Logger name (typically the module name, e.g. __name__).
        enable_console_logging : bool
            Enable console output (default: True).
        console_verbose : bool
            Display detailed log messages on the console (default: False).
            If True, sets the console log level to INFO; if False, to WARNING.
        enable_file_logging : bool
            Enable file output (default: False).
        log_dir : str
            Directory for log files (default: './logs').
        log_file_name : Optional[str]
            Log file name (default: <logger_name>.log).
        log_level : int
            Base log level (default: logging.INFO).
        """

        # # --- デバッグプリント追加 ---
        # import sys # print を stderr に出すため
        # print(f"DEBUG_LM: Initializing LoggerManager for '{logger_name}'", file=sys.stderr)
        # print(f"DEBUG_LM: enable_console_logging={enable_console_logging}, console_verbose={console_verbose}", file=sys.stderr)
        # print(f"DEBUG_LM: enable_file_logging={enable_file_logging}", file=sys.stderr) # ★これが True になっているか？
        # print(f"DEBUG_LM: log_file_name received = '{log_file_name}'", file=sys.stderr) # ★指定したファイル名が渡されているか？
        # print(f"DEBUG_LM: log_level={log_level}", file=sys.stderr)
        # # --- デバッグプリントここまで ---

        # ロガーの取得
        self.logger = logging.getLogger(logger_name)
        self.logger.setLevel(log_level)
        
        # 既存のハンドラを削除（再設定時の重複を避けるため）
        # この部分は選択的に実装可能
        # for handler in self.logger.handlers[:]:
        #     self.logger.removeHandler(handler)

        # フォーマット設定
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

        # コンソールハンドラ
        if enable_console_logging:
            console_level = logging.INFO if console_verbose else logging.WARNING
            
            # 既存のStreamHandlerをチェック
            if not any(isinstance(handler, logging.StreamHandler) and 
                      not isinstance(handler, logging.FileHandler) 
                      for handler in self.logger.handlers):
                console_handler = logging.StreamHandler(sys.stdout)
                console_handler.setLevel(console_level)
                console_handler.setFormatter(formatter)
                self.logger.addHandler(console_handler)

        # ファイルハンドラ
        if enable_file_logging:
            # print(f"DEBUG_LM: File logging IS enabled. Proceeding to add handler.", file=sys.stderr) # ★ここまで到達するか？
            # ログファイルパスの決定ロジック
            log_file_path = None # 初期化
            if log_file_name:
                log_file_path = log_file_name
                # print(f"DEBUG_LM: Using provided log_file_name: {log_file_path}", file=sys.stderr)
                log_file_dir = os.path.dirname(log_file_path)
                if log_file_dir and not os.path.exists(log_file_dir):
                    # print(f"DEBUG_LM: Creating directory {log_file_dir}", file=sys.stderr)
                    os.makedirs(log_file_dir, exist_ok=True)
            elif logger_name:
                if not os.path.exists(log_dir):
                    # print(f"DEBUG_LM: Creating directory {log_dir}", file=sys.stderr)
                    os.makedirs(log_dir, exist_ok=True)
                default_log_filename = f"{logger_name.split('.')[-1]}.log"
                log_file_path = os.path.join(log_dir, default_log_filename)
                # print(f"DEBUG_LM: Using default log path: {log_file_path}", file=sys.stderr)
            else:
                if not os.path.exists(log_dir):
                    # print(f"DEBUG_LM: Creating directory {log_dir}", file=sys.stderr)
                    os.makedirs(log_dir, exist_ok=True)
                log_file_path = os.path.join(log_dir, "default.log")
                # print(f"DEBUG_LM: Using fallback log path: {log_file_path}", file=sys.stderr)

            # log_file_path が決定されたか確認
            if log_file_path:
                # print(f"DEBUG_LM: Determined log_file_path: {os.path.abspath(log_file_path)}", file=sys.stderr) # ★パスを確認
                # 既存ハンドラのチェック
                handler_exists = any(isinstance(handler, logging.handlers.RotatingFileHandler) and
                                hasattr(handler, 'baseFilename') and
                                os.path.abspath(handler.baseFilename) == os.path.abspath(log_file_path)
                                for handler in self.logger.handlers)
                # print(f"DEBUG_LM: File handler with path '{os.path.abspath(log_file_path)}' already exists? {handler_exists}", file=sys.stderr) # ★既存チェックの結果

                if not handler_exists:
                    # print(f"DEBUG_LM: Attempting to add RotatingFileHandler for {log_file_path}", file=sys.stderr) # ★追加処理に入るか？
                    try:
                        file_handler = logging.handlers.RotatingFileHandler(
                            log_file_path,
                            maxBytes=10**6,
                            backupCount=5
                        )
                        file_handler.setLevel(log_level)
                        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
                        file_handler.setFormatter(formatter)
                        self.logger.addHandler(file_handler)
                        # print(f"DEBUG_LM: SUCCESSFULLY added RotatingFileHandler for {log_file_path}", file=sys.stderr) # ★成功したか？
                    except Exception as e:
                        pass
                        # print(f"DEBUG_LM: FAILED to add RotatingFileHandler for {log_file_path}: {e}", file=sys.stderr) # ★失敗したか？
                else:
                    pass
                    # print(f"DEBUG_LM: Skipping add handler because it already exists.", file=sys.stderr)
            else:
                pass
                # print(f"DEBUG_LM: log_file_path could not be determined. Skipping file handler.", file=sys.stderr)
        else:
            # print(f"DEBUG_LM: File logging is NOT enabled.", file=sys.stderr) # ★ここに来ていないか？
            pass
        
        self.logger.propagate = False



    def get_logger(self) -> logging.Logger:
        """Retrieves the configured logger."""
        return self.logger


# モジュールレベルのヘルパー関数
# def get_module_logger(module_name=None, console=True, file_output=False, log_level=logging.INFO):


# def get_module_logger(...): の引数を変更
def get_module_logger(module_name=None, console=True, file_output=False, log_level=logging.INFO, log_file=None): # log_file 引数を追加
    """
    Conveniently retrieves a logger for a module.
    
    Parameters:
    -----------
    module_name : str
        Logger name. If None, the caller's module name is used.
    console : bool
        Enable console output (default: True).
    file_output : bool
        Enable file output (default: False).
    log_level : int
        Logging level (default: logging.INFO).
    
    Returns:
    --------
    logging.Logger: The configured logger instance.
    
    Example:
    ```
    logger = get_module_logger(__name__, file_output=True)
    logger.info("Message")
    ```
    """
    import inspect
    
    # モジュール名が指定されていない場合、呼び出し元から取得
    if module_name is None:
        frame = inspect.stack()[1]
        module = inspect.getmodule(frame[0])
        module_name = module.__name__
    
    # LoggerManagerを使ってロガーを設定して返す
    manager = LoggerManager(
        logger_name=module_name,
        enable_console_logging=console,
        enable_file_logging=file_output,
        log_level=log_level,
        log_file_name=log_file # log_file を log_file_name として渡す
    )
    return manager.get_logger()
