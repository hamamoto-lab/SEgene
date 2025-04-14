# SEgene_peakprep - CPM方式の詳細ガイド

このドキュメントは、SEgene_peakprepの**CPM方式**に特化した詳細な使用方法と実行例を提供します。基本的な概要や共通要件については、[メインREADME](./README_ja.md)を参照してください。

## 全引数リスト

コマンドライン実行時に使用できる全ての引数とその詳細です：

| 引数名 | 短縮形 | デフォルト値 | 説明 |
|--------|--------|------------|------|
| **必須引数** ||||
| `--bam_dir` | `-b` | - | BAMファイルが格納されているディレクトリのパス |
| `--annotation_file` | `-a` | - | アノテーションファイルのパス（BED, SAF, または merge_SE.tsv 形式） |
| `--output_dir` | `-o` | - | 結果を保存するディレクトリ |
| **アノテーション関連オプション** ||||
| `--is_mergese_format` | - | `False` | アノテーションファイルがmerge_SE.tsv形式であることを指定するフラグ |
| **ファイル/サンプル選択オプション** ||||
| `--filename_pattern` | - | `*.bam` | BAMファイル名をフィルタリングするためのワイルドカードパターン |
| `--sample_delimiter` | - | `None` | BAMファイル名からサンプル名を抽出するための区切り文字列 |
| **ツールパスとスレッド数** ||||
| `--samtools_path` | - | `samtools` | samtoolsの実行ファイルへのパス |
| `--featurecounts_path` | - | `featureCounts` | featureCountsの実行ファイルへのパス |
| `--threads` | `-T` | `4` | samtoolsとfeatureCountsで使用するスレッド数 |
| `--single_end` | - | `False` | シングルエンドのリードを指定するフラグ（指定しない場合はペアエンド） |
| **featureCounts関連オプション** ||||
| `--fc_basename` | - | `None` | featureCounts出力/ログファイルのベース名 |
| `--fc_options` | - | `None` | featureCountsに渡す追加オプション（例：`--fc_options --minOverlap 10`） |
| **logCPM計算関連オプション** ||||
| `--add_peakid` | - | `False` | オリジナルのGeneIDを連番のPeakIDに置き換え |
| `--id_column` | - | `Geneid` | SAF/featureCounts出力のID列名 |
| `--output_id_name` | - | `PeakID` | 最終的なlogCPMテーブルのID列名 |
| `--log_base` | - | `2` | 対数変換の底（2または10）。0以下で変換なし |
| `--pseudocount` | - | `1.0` | 対数変換前に追加する擬似カウント |
| **出力ファイル名オプション** ||||
| `--logcpm_output_name` | - | `logCPM.tsv` | logCPMテーブルのファイル名 |
| **ログ関連オプション** ||||
| `--log_level` | - | `INFO` | ロギングレベル（DEBUG/INFO/WARNING/ERROR/CRITICAL） |
| `--script_log_file` | - | `pipeline_run.log` | スクリプト実行ログのファイル名（空文字で無効化） |

## 高度なユースケース例

### featureCountsへの追加オプション指定

特定のオーバーラップ基準やカウント方法を指定する例：

```bash
python cpm_peakprep.py \
    --bam_dir /path/to/bam_files \
    --annotation_file peaks.saf \
    --output_dir results \
    --fc_options --minOverlap 10 --fracOverlap 0.2 --ignoreDup \
    --threads 8
```

### カスタムID生成と対数変換設定

PeakIDを連番生成し、カスタム対数変換を行う例：

```bash
python cpm_peakprep.py \
    --bam_dir /path/to/bam_files \
    --annotation_file peaks.bed \
    --output_dir results \
    --add_peakid \
    --output_id_name "CustomPeakID" \
    --log_base 10 \
    --pseudocount 0.1
```

## Jupyter Notebookでの実行例

### シンプルな使用法

最もシンプルな形では、メインのスクリプトをインポートして使用できます：

```python
import sys
sys.path.append('/path/to/SEgene/SEgene_peakprep')

# モジュールとして実行
import cpm_peakprep

# パラメータ設定
args = [
    "--bam_dir", "/path/to/bam_files",
    "--annotation_file", "peaks.bed",
    "--output_dir", "results",
    "--sample_delimiter", ".sorted.bam",
    "--threads", "8"
]

# main関数を直接呼び出し
sys.argv = ["cpm_peakprep.py"] + args
cpm_peakprep.main()
```

### 詳細なstep by step実行例

以下は、各ステップを個別に制御して実行する詳細な例です：

```python
import os
import sys
import logging
import pandas as pd

# モジュールがインポートできるようにパスを設定
sys.path.append('/path/to/SEgene/SEgene_peakprep')

# 必要な関数を個別にインポート
from cpm_peakprep_utils import (
    get_sample_data,
    run_samtools_flagstat,
    run_featurecounts,
    calculate_logcpm,
    convert_bed_to_saf
)

# パラメータ設定
bam_dir = "/path/to/bam_files"
annotation_file = "peaks.bed"
output_dir = "results"
filename_pattern = "Sample*.sorted.bam"
sample_delimiter = ".sorted.bam"
threads = 10

# 出力ディレクトリを作成
os.makedirs(output_dir, exist_ok=True)

# ロガーの設定
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger("ChIPseqPipeline")

# BEDファイルをSAFに変換（必要な場合）
saf_file = annotation_file
if annotation_file.lower().endswith('.bed'):
    saf_file = convert_bed_to_saf(annotation_file, logger)
    if saf_file is None:
        print("BEDファイルの変換に失敗しました。")
        sys.exit(1)

# ステップ1: BAMファイル情報の取得
sample_dict, bam_list = get_sample_data(
    bam_folder=bam_dir,
    logger=logger,
    sample_name_delimiter=sample_delimiter,
    filename_pattern=filename_pattern
)

# ステップ2: 総リード数の取得
flagstat_dir = os.path.join(output_dir, "flagstat")
total_reads_map = run_samtools_flagstat(
    bam_files=bam_list,
    output_dir=flagstat_dir,
    logger=logger,
    threads=threads
)

# ステップ3: リードカウント（featureCounts）
fc_output = os.path.join(output_dir, f"{os.path.splitext(os.path.basename(saf_file))[0]}_featureCounts.txt")
fc_log = os.path.join(output_dir, f"{os.path.splitext(os.path.basename(saf_file))[0]}_featureCounts.log")
counts_file_path = run_featurecounts(
    bam_files=bam_list,
    saf_file=saf_file,
    logger=logger,
    output_file=fc_output,
    log_file=fc_log,
    threads=threads,
    is_paired_end=True
)

# カウントファイルの読み込みとステップ4: LogCPM計算
counts_df = pd.read_csv(counts_file_path, sep='\t', comment='#', header=0)
if counts_df.columns[0].startswith('#'):
    counts_df = pd.read_csv(counts_file_path, sep='\t', comment='#', header=1)

bam_to_sample_map = {v: k for k, v in sample_dict.items()}
logcpm_df = calculate_logcpm(
    counts_df=counts_df,
    total_reads_dict=total_reads_map,
    logger=logger,
    bam_path_to_sample_name=bam_to_sample_map,
    log_transform_base=2,
    pseudocount=1.0
)

# logCPM結果の保存
logcpm_output_file = os.path.join(output_dir, "logCPM.tsv")
logcpm_df.to_csv(logcpm_output_file, sep='\t', index=False, float_format='%.4f')
print(f"logCPMテーブルが保存されました: {logcpm_output_file}")

# 一時ファイルのクリーンアップ
if 'saf_file' in locals() and saf_file != annotation_file and os.path.exists(saf_file):
    try:
        os.remove(saf_file)
        print(f"一時SAFファイルを削除しました: {saf_file}")
    except Exception as e:
        print(f"一時ファイルの削除に失敗しました: {e}")
```

## 中間ファイルと出力形式

### flagstat出力

`flagstat/` ディレクトリに保存される各BAMファイルのsamtools flagstat結果から、マップされたリード数が抽出され、CPM計算に使用されます。

### featureCounts出力

featureCountsの出力例：
```
# Program:featureCounts v2.0.6; Command:...
Geneid  Chr     Start   End     Strand  Length  bam_files/sample1.bam  bam_files/sample2.bam
peak1   chr1    1000    2000    .       1001    150     220
peak2   chr1    3000    4000    .       1001    342     289
...
```

### logCPMテーブル

最終的なlogCPMテーブル例：
```
PeakID  Chr     Start   End     Sample1  Sample2  Sample3
peak1   chr1    999     2000    2.45     3.12     1.87
peak2   chr1    2999    4000    4.21     3.98     4.56
...
```

> **注意**: すべての出力ファイルの `Start` 列は0-based（BED形式基準）座標です。入力がSAF形式（1-based）の場合でも、出力は0-basedに変換されます。

## 引用

このツールを研究に使用する場合は、以下を引用してください：
**(論文は現在準備中です)**

## ライセンス

このプログラムはMITライセンスで公開されています。詳細については、[LICENSE](./LICENSE)をご覧ください。