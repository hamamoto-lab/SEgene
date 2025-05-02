# cpm_peakprep - CPM方式の詳細ガイド

このドキュメントは、SEgene_peakprepの**CPM方式**に特化した詳細な使用方法と実行例を提供します。基本的な概要や共通要件については、[メインREADME](./README_ja.md)を参照してください。

## 全引数リスト

コマンドライン実行時に使用できる全ての引数とその詳細です：

| 引数名 | 説明 | デフォルト値 |
|--------|------|-------------|
| **必須引数** |||
| `--bam_dir` / `-b` | BAMファイルが格納されているディレクトリのパス | - |
| `--annotation_file` / `-a` | アノテーションファイルのパス（BED, SAF, または merge_SE.tsv 形式） | - |
| `--output_dir` / `-o` | 結果を保存するディレクトリ | - |
| **アノテーション関連オプション** |||
| `--is_mergese_format` | アノテーションファイルがmerge_SE.tsv形式であることを指定するフラグ | `False` |
| **ファイル/サンプル選択オプション** |||
| `--filename_pattern` | BAMファイル名をフィルタリングするためのワイルドカードパターン | `*.bam` |
| `--sample_delimiter` | BAMファイル名からサンプル名を抽出するための区切り文字列 | `None` |
| **ツールパスとスレッド数** |||
| `--samtools_path` | samtoolsの実行ファイルへのパス | `samtools` |
| `--featurecounts_path` | featureCountsの実行ファイルへのパス | `featureCounts` |
| `--threads` / `-T` | samtoolsとfeatureCountsで使用するスレッド数 | `4` |
| `--single_end` | シングルエンドのリードを指定するフラグ（指定しない場合はペアエンド） | `False` |
| **featureCounts関連オプション** |||
| `--fc_basename` | featureCounts出力/ログファイルのベース名 | `None` |
| `--fc_options` | featureCountsに渡す追加オプション（例：`--fc_options --minOverlap 10`） | `None` |
| **Standard log2-CPM計算関連オプション** |||
| `--add_peakid` | オリジナルのGeneIDを連番のPeakIDに置き換え | `False` |
| `--id_column` | SAF/featureCounts出力のID列名 | `Geneid` |
| `--output_id_name` | 最終的なlogCPMテーブルのID列名 | `PeakID` |
| `--log_base` | 対数変換の底（2または10）。0以下で変換なし | `2` |
| `--pseudocount` | 対数変換前に追加する擬似カウント | `1.0` |
| **出力ファイル名オプション** |||
| `--logcpm_output_name` | log2-CPMテーブルのファイル名 | `logCPM.tsv` |
| **ログ関連オプション** |||
| `--log_level` | ロギングレベル（DEBUG/INFO/WARNING/ERROR/CRITICAL） | `INFO` |
| `--script_log_file` | スクリプト実行ログのファイル名（空文字で無効化） | `pipeline_run.log` |

## 依存パッケージとバージョン要件

### Python
- **Python**: ≥3.10
- **pandas**: ≥1.5
- **numpy**: ≥1.23
- **natsort**: ≥8.3

### 標準CPM計算のみで必要
- **samtools**: ≥1.13
- **featureCounts** (Subread): ≥2.0.0

### edgeR正規化CPM（calcnorm）で追加で必要
- **R言語**: ≥4.2.0
- **edgeR**: ≥3.40.0
- **rpy2**: ≥3.5.0

## edgeR正規化CPM（calcnorm）モード

`--calcnorm`フラグを付けると、featureCountsから得られた生カウントデータに対して、edgeRの`calcNormFactors()`関数を適用し、正規化されたライブラリサイズに基づくCPM値を計算します。この機能により、サンプル間の組成バイアスを補正した正規化が可能になります。

デフォルトでは**Standard log2-CPM**（`logCPM.tsv`）と**calcnorm CPM**（`calcnorm.tsv`）の両方を出力します。それぞれの出力は`--logcpm_output_name`と`--calcnorm_output_name`で制御できます。

| 引数名 | 説明 | デフォルト値 |
|--------|------|-------------|
| `--calcnorm` | edgeR正規化CPMを有効化 | `False` |
| `--calcnorm_method` | 正規化手法（upperquartile, TMM, RLE, none） | `upperquartile` |
| `--min_cpm` | フィルタリングのための最小CPM閾値 | `1.0` |
| `--min_samples` | CPM閾値を超えるべき最小サンプル数（0=フィルタリングなし） | `0` |
| `--calcnorm_output_name` | calcnorm CPM出力ファイル名（空文字列で保存しない） | `calcnorm.tsv` |
| `--remove_extensions` | 一般的なBAMファイル拡張子（.bam, .sorted.bamなど）を除去 | `False` |
| `--pattern_rules` | サンプル名クリーニングのための追加パターンルールを含むJSONファイル | `None` |
| `--full_metadata_output` | 標準出力に加えて全featureCountsメタデータ列を含む出力を保存 | `None` |

### 正規化手法の選択

`--calcnorm_method`オプションでは、edgeRパッケージで実装されている以下の正規化手法を選択できます：

- **upperquartile**（デフォルト）: 上位四分位数正規化。ゼロでないカウント値の75パーセンタイルに基づく正規化で、ChIP-seqデータのような疎なデータセットに適しています。
- **TMM**: Trimmed Mean of M-values。サンプル間のログ比のトリム平均に基づく正規化で、RNA-seqデータ解析でよく使用されます。
- **RLE**: Relative Log Expression。遺伝子ごとの幾何平均を基準とした相対的な発現量に基づく正規化です。
- **none**: 正規化を適用せず、ライブラリサイズの違いのみを考慮したCPMを計算します。

各正規化手法の詳細な理論的背景については、edgeRのドキュメントを参照してください。calcnorm CPMの詳細な使用方法や背景については、[edgeR正規化CPMガイド](./cpm_calcnorm_README_ja.md)に記載されています。

## CPM値計算の違い

### Standard log2-CPM と edgeR正規化CPMの計算の違い

1. **Standard log2-CPM**:
   - samtoolsのflagstatで取得した総マップリード数に基づいて計算
   - 計算式: CPM = (count × 10⁶) / total_mapped_reads
   - その後、log₂(CPM + pseudocount)で対数変換

2. **edgeR正規化CPM（calcnorm CPM）**:
   - featureCountsで計数された該当領域のリードに基づく計算
   - edgeRの正規化係数を用いてサンプル間のライブラリサイズを調整
   - サンプル組成の違いを考慮した正規化が可能

この違いにより、特に組成の異なるサンプル間を比較する場合などに結果が異なる場合があります。

## 基本的な使用例

### シンプルな実行例

```bash
python cpm_peakprep.py \
    --bam_dir /path/to/bam_files \
    --annotation_file peaks.bed \
    --output_dir results \
    --filename_pattern "Sample*.sorted.bam" \
    --sample_delimiter=".sorted.bam" \
    --threads 10
```

### featureCountsへの追加オプション指定

```bash
python cpm_peakprep.py \
    --bam_dir /path/to/bam_files \
    --annotation_file peaks.saf \
    --output_dir results \
    --fc_options --minOverlap 10 --fracOverlap 0.2 --ignoreDup \
    --threads 8
```

### カスタムID生成と対数変換設定

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

### edgeR正規化CPMを使用する例

```bash
python cpm_peakprep.py \
    --bam_dir /path/to/bam_files \
    --annotation_file peaks.bed \
    --output_dir results \
    --calcnorm \
    --calcnorm_method TMM \
    --threads 10
```

### calcnorm CPMのみを保存する例

```bash
python cpm_peakprep.py \
    --bam_dir /path/to/bam_files \
    --annotation_file peaks.bed \
    --output_dir results \
    --calcnorm \
    --logcpm_output_name "" \
    --threads 10
```

### 発現フィルタリングを適用した例

```bash
python cpm_peakprep.py \
    --bam_dir /path/to/bam_files \
    --annotation_file peaks.bed \
    --output_dir results \
    --calcnorm \
    --min_cpm 1.0 \
    --min_samples 3 \
    --threads 10
```

この例では、少なくとも3つのサンプルでCPM値が1.0以上のピークのみが保持されます。

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

### log2-CPMテーブル

最終的なlog2-CPMテーブル例：
```
PeakID  Chr     Start   End     Sample1  Sample2  Sample3
peak1   chr1    999     2000    2.45     3.12     1.87
peak2   chr1    2999    4000    4.21     3.98     4.56
...
```

> **注意**: すべての出力ファイルの `Start` 列は0-based（BED形式基準）座標です。入力がSAF形式（1-based）の場合でも、出力は0-basedに変換されます。

## 関連ドキュメント

- [メインREADME](./README_ja.md): ツール全体の概要
- [edgeR正規化CPMガイド](./cpm_calcnorm_README_ja.md): edgeR正規化CPMの詳細情報
- [BigWig方式の詳細ガイド](./bigwig_README_ja.md): BigWig方式の詳細情報

## ライセンス

このプログラムはMITライセンスで公開されています。詳細については、[LICENSE](https://github.com/hamamoto-lab/SEgene/blob/main/LICENSE)をご覧ください。