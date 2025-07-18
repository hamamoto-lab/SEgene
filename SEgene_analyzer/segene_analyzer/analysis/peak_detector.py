# sedb_tools/analysis/peak_detector.py
import numpy as np
from scipy import signal
from scipy.ndimage import gaussian_filter1d
from collections import defaultdict
import logging

class SEdbPeakDetector:
    """
    SEdbデータからのピーク検出を専門とするクラス
    
    このクラスは、スーパーエンハンサー読み取りデータからピークを検出し、
    ピーク位置、サンプル分布、読み取り密度などの統計を計算します。
    """
    
    def __init__(self, logger=None):
        """
        初期化
        
        Parameters:
        -----------
        logger : logging.Logger, optional
            ロガーインスタンス。Noneの場合はデフォルトロガーを作成
        """
        self.logger = logger or self._setup_default_logger()
    
    def _setup_default_logger(self):
        """デフォルトロガーを設定"""
        logger = logging.getLogger(f"{__name__}.{id(self)}")
        if not logger.handlers:
            logger.setLevel(logging.INFO)
            handler = logging.StreamHandler()
            formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
            handler.setFormatter(formatter)
            logger.addHandler(handler)
            logger.propagate = False
        return logger
    
    def detect_peaks(self, reads, start, end, bin_size=50, peak_threshold=0.5, 
                    min_peak_distance=200, peak_width_factor=0.5, smooth_window=3):
        """
        読み取りデータからピークを検出
        
        Parameters:
        -----------
        reads : list
            読み取りデータを表す辞書のリスト。各辞書には 'start', 'end', 'sample' キーが必要
        start : int
            解析領域の開始位置
        end : int
            解析領域の終了位置
        bin_size : int, optional
            ビンサイズ（デフォルト: 50）
        peak_threshold : float, optional
            ピーク検出の閾値（0-1、デフォルト: 0.5）
        min_peak_distance : int, optional
            ピーク間の最小距離（bp、デフォルト: 200）
        peak_width_factor : float, optional
            ピーク幅を決定する要素（デフォルト: 0.5）
        smooth_window : int, optional
            平滑化ウィンドウサイズ（デフォルト: 3）
            
        Returns:
        --------
        dict
            検出されたピーク情報を含む辞書
        dict
            ビンカウント情報
        list
            ビンサンプル情報
        """
        # 領域範囲の計算
        region_range = end - start
        num_bins = int(region_range / bin_size) + 1
        bins = np.linspace(start, end, num_bins)
        bin_counts = np.zeros(num_bins - 1)
        
        # サンプルごとのビン追跡
        sample_bins = defaultdict(set)
        bin_samples = [set() for _ in range(num_bins - 1)]
        
        # 各ビンの読み取り数をカウント
        for read in reads:
            read_start_bin = max(0, int((read['start'] - start) / bin_size))
            read_end_bin = min(num_bins - 2, int((read['end'] - start) / bin_size))
            
            sample_name = read['sample']
            
            for bin_idx in range(read_start_bin, read_end_bin + 1):
                bin_counts[bin_idx] += 1
                bin_samples[bin_idx].add(sample_name)
                sample_bins[sample_name].add(bin_idx)
        
        # 元のビンカウントを保存
        original_bin_counts = bin_counts.copy()
        
        # ビンカウントの平滑化
        if smooth_window > 1:
            smoothed_counts = gaussian_filter1d(bin_counts, sigma=smooth_window)
        else:
            smoothed_counts = bin_counts.copy()
        
        # 検出のための正規化
        if smoothed_counts.max() > 0:
            normalized_counts = smoothed_counts / smoothed_counts.max()
        else:
            normalized_counts = smoothed_counts.copy()
        
        # ピーク検出
        actual_threshold = peak_threshold * normalized_counts.max()
        peaks, peak_properties = signal.find_peaks(
            normalized_counts, 
            height=actual_threshold,
            distance=int(min_peak_distance / bin_size)
        )
        
        # ピーク位置のゲノム座標への変換
        peak_positions = start + peaks * bin_size
        peak_heights = peak_properties['peak_heights']
        
        # 閾値に基づくピーク幅の計算
        peak_widths = signal.peak_widths(
            normalized_counts, peaks, rel_height=peak_width_factor
        )
        
        # 幅インデックスのゲノム座標への変換
        peak_left_idx = peak_widths[2].astype(int)
        peak_right_idx = peak_widths[3].astype(int)
        peak_left_pos = start + peak_left_idx * bin_size
        peak_right_pos = start + peak_right_idx * bin_size
        
        # 総サンプル数
        total_samples = len(set(read['sample'] for read in reads))
        
        # ピーク情報の生成
        peak_data = {}
        for i, peak_pos in enumerate(peak_positions):
            peak_id = f"peak_{i+1}"
            
            # ピーク領域内のビンインデックスを取得
            peak_bin_indices = list(range(
                max(0, peak_left_idx[i]),
                min(num_bins - 1, peak_right_idx[i] + 1)
            ))
            
            # ピーク領域内の読み取り統計を計算
            if peak_bin_indices:
                peak_region_counts = original_bin_counts[peak_bin_indices]
                peak_total_reads = int(peak_region_counts.sum())
                peak_max_reads = int(peak_region_counts.max())
                peak_avg_reads = float(peak_region_counts.mean())
            else:
                peak_total_reads = 0
                peak_max_reads = 0
                peak_avg_reads = 0
            
            # ピーク領域内のサンプル統計
            peak_samples = set()
            for bin_idx in peak_bin_indices:
                if 0 <= bin_idx < len(bin_samples):
                    peak_samples.update(bin_samples[bin_idx])
            
            peak_sample_count = len(peak_samples)
            peak_sample_percentage = (peak_sample_count / total_samples * 100) if total_samples > 0 else 0
            
            # ピーク幅
            peak_width = int(peak_right_pos[i] - peak_left_pos[i])
            
            # 読み取り密度
            density = peak_total_reads / peak_width if peak_width > 0 else 0
            
            # 強化されたピークデータを保存
            peak_data[peak_id] = {
                'position': int(peak_pos),
                'start': int(peak_left_pos[i]),
                'end': int(peak_right_pos[i]),
                'width': peak_width,
                'reads': {
                    'total': peak_total_reads,
                    'max': peak_max_reads,
                    'average': round(peak_avg_reads, 2)
                },
                'samples': {
                    'count': peak_sample_count,
                    'percentage': round(peak_sample_percentage, 2),
                    'list': sorted(list(peak_samples))
                },
                'normalized_height': float(peak_heights[i]),
                'density': round(density, 3)
            }
                
        self.logger.info(f"Detected {len(peak_data)} peaks in region {start}-{end}")
        return peak_data, original_bin_counts, bin_samples
    
    def get_peak_display_data(self, peak_data):
        """
        元の表示フォーマットと互換性のあるピークデータを作成
        
        Parameters:
        -----------
        peak_data : dict
            detect_peaks()から返されるピークデータ
            
        Returns:
        --------
        dict
            元のコードと互換性のある表示用データ
        """
        peak_display_data = {}
        for peak_id, peak in peak_data.items():
            peak_display_data[peak_id] = {
                'position': peak['position'],
                'height': peak['normalized_height'],
                'left': peak['start'],
                'right': peak['end'],
                'width': peak['width']
            }
        return peak_display_data