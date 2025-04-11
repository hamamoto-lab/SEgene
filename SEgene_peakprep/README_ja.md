# SEgene_peakprep

*(For the English version of this README, please see [README.md](./README.md).)*

**SEgene_peakprep** は、シーケンスデータ（BAMファイル）から指定されたゲノム領域における読み取りカウントを正規化し、CPM（Counts Per Million）値を計算するためのPythonパイプラインです。このツールは、特定の参照ゲノムには依存せず、例えばhg38などにマッピングされたデータで利用することができます。

## 開発状況

> **注意**: このツールは**現在開発中**であり、プログラムのファイル構成や使用方法はバージョンによって変化する可能性があります。

## 概要

SEgene_peakprepは、以下の機能を提供します：

- 複数のBAMファイルの効率的な処理とサンプル情報の自動抽出
- `samtools flagstat`を使用した総マップ済みリード数の計算
- `featureCounts`を使用した指定ゲノム領域（BED, SAF, または merge_SV.tsv 形式）におけるリードカウント
- CPM値の計算とログ変換（デフォルトはlog2(CPM+1)）
- 結果の整形と保存

**SEgeneプロジェクト内での位置づけ：**

このパイプラインは[SEgeneプロジェクト](https://github.com/hamamoto-lab/SEgene)において、初期データ準備の役割を果たします：

1. **入力**: マッピング済みシーケンスデータ（BAMファイル）とアノテーションファイル（BED/SAF/merge_SV.tsv形式）
2. **出力**: CPM/logCPM値テーブル（TSVファイル）
3. **フロー**: この出力は、`peak_to_gene_links`モジュールの主要な入力となり、RNA-seqのTPMデータと組み合わせて遺伝子制御ネットワーク解析に使用されます

## 要件

- **Python**: Python 3.7以上
- **必須ライブラリ**:
  - pandas
  - numpy
  - natsort
  
  ※必要なライブラリはcondaまたはpipを使用してインストールしてください：
  ```bash
  # condaの場合
  conda install pandas numpy natsort
  
  # pipの場合
  pip install pandas numpy natsort
  ```
  
- **外部ツール**:
  - samtools (推奨: v1.13以上)
  - featureCounts (Subreadパッケージの一部、推奨: v2.0.0以上)
  
  ※これらのツールはPATHを通すか、実行時にオプションでパスを指定してください。
  ```bash
  # condaでのインストール例
  conda install -c bioconda samtools subread
  ```
  
  詳細なインストール方法については、[samtools](http://www.htslib.org/)と[Subread/featureCounts](http://subread.sourceforge.net/)の公式サイトを参照してください。

## インストール

```bash
# リポジトリのクローン
git clone https://github.com/hamamoto-lab/SEgene.git
cd SEgene/SEgene_peakprep
```

## 入力データの準備

### アノテーションファイルの準備
このパイプラインは、以下の形式のアノテーションファイルを入力として受け付けます：

1. **BED形式**： 一般的なBEDフォーマット（0-based開始座標）。パイプラインは内部で自動的にSAF形式に変換します。

2. **SAF形式** (Simplified Annotation Format)： タブ区切りのテキストファイルで、以下の列を含みます：
   - `GeneID` - 各領域の一意な識別子
   - `Chr` - 染色体名
   - `Start` - 開始位置（1-based）
   - `End` - 終了位置
   - `Strand` - ストランド（+、-、または.）

   例：
   ```
   GeneID  Chr     Start   End     Strand
   peak1   chr1    1000    2000    .
   peak2   chr1    3000    4000    .
   ```

3. **merge_SV.tsv形式**： 1列目が `chr_start_end` 形式の領域ID（例：`chr1_1000_2000`）を含むタブ区切りファイル。この形式を使用する場合は `--is_mergesv_format` フラグを指定する必要があります。

アノテーションファイルは、`--annotation_file` 引数で指定します。BEDファイルやmerge_SVファイルを指定した場合、実行時に自動的にSAF形式に変換されます。

### BAMファイルの準備
BAMファイルは標準的なシーケンスマッピングパイプラインを使用して生成してください。BAMファイルはインデックス付け（.baiファイル）されていることを推奨します。

## 使用方法

パイプラインは、コマンドラインまたはJupyterノートブックから実行できます。

### コマンドライン実行

**BEDファイルを使用した実行コマンド（基本）：**

```bash
python cpm_peakprep.py \
    --bam_dir /path/to/bam_files \
    --annotation_file peaks.bed \
    --output_dir results \
    --filename_pattern "Sample*.sorted.bam" \
    --sample_delimiter=".sorted.bam" \
    --threads 10
```

**SAFファイルを使用した実行コマンド：**

```bash
python cpm_peakprep.py \
    --bam_dir /path/to/bam_files \
    --annotation_file peaks.saf \
    --output_dir results \
    --filename_pattern "Sample*.sorted.bam" \
    --sample_delimiter=".sorted.bam" \
    --threads 10
```

**merge_SV.tsv形式のファイルを使用した実行コマンド：**

```bash
python cpm_peakprep.py \
    --bam_dir /path/to/bam_files \
    --annotation_file merged_se.tsv \
    --is_mergesv_format \
    --output_dir results \
    --filename_pattern "Sample*.sorted.bam" \
    --sample_delimiter=".sorted.bam" \
    --threads 10
```

> **注意**: 値にハイフンを含む場合は、イコール記号（`=`）を使用して指定する必要があります（例：`--sample_delimiter="-H3K27ac.sorted.bam"`）。

#### 引数の説明

| 引数名 | 短縮形 | デフォルト値 | 説明 |
|--------|--------|------------|------|
| **必須引数** ||||
| `--bam_dir` | `-b` | - | BAMファイルが格納されているディレクトリのパス |
| `--annotation_file` | `-a` | - | アノテーションファイルのパス（BED, SAF, または merge_SV.tsv 形式） |
| `--output_dir` | `-o` | - | 結果を保存するディレクトリ |
| **アノテーション関連オプション** ||||
| `--is_mergesv_format` | - | `False` | アノテーションファイルがmerge_SV.tsv形式であることを指定するフラグ |
| **ファイル/サンプル選択オプション** ||||
| `--filename_pattern` | - | `*.bam` | BAMファイル名をフィルタリングするためのワイルドカードパターン |
| `--sample_delimiter` | - | `None` | BAMファイル名からサンプル名を抽出するための区切り文字列 |
| **ツールパスとスレッド数** ||||
| `--samtools_path` | - | `samtools` | samtoolsの実行ファイルへのパス |
| `--featurecounts_path` | - | `featureCounts` | featureCountsの実行ファイルへのパス |
| `--threads` | `-T` | `4` | samtoolsとfeatureCountsで使用するスレッド数 |
| `--single_end` | - | `False` | シングルエンドのリードを指定するフラグ（指定しない場合はペアエンドと見なされる） |
| **featureCounts関連オプション** ||||
| `--fc_basename` | - | `None` | featureCounts出力/ログファイルのベース名（デフォルト：`ANNOTATION_fc`） |
| `--fc_options` | - | `None` | featureCountsに渡す追加オプション（例：`--fc_options --minOverlap 10`） |
| **logCPM計算関連オプション** ||||
| `--add_peakid` | - | `False` | オリジナルのGeneIDを連番のPeakIDに置き換えるフラグ |
| `--id_column` | - | `Geneid` | SAF/featureCounts出力のID列名 |
| `--output_id_name` | - | `PeakID` | 最終的なlogCPMテーブルのID列名 |
| `--log_base` | - | `2` | 対数変換の底（2または10）。0以下を設定すると変換なし |
| `--pseudocount` | - | `1.0` | 対数変換前に追加する擬似カウント |
| **出力ファイル名オプション** ||||
| `--logcpm_output_name` | - | `logCPM.tsv` | logCPMテーブルのファイル名 |
| **ログ関連オプション** ||||
| `--log_level` | - | `INFO` | ロギングレベル（DEBUG/INFO/WARNING/ERROR/CRITICAL） |
| `--script_log_file` | - | `pipeline_run.log` | スクリプト実行ログのファイル名（空文字で無効化） |

*(`デフォルト値は変更される可能性があります。最新の情報については `python cpm_peakprep.py --help` で確認してください`)*

### Jupyterノートブックでの実行

以下は、Jupyterノートブックでパイプラインを実行する基本的な例です：

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
annotation_file = "peaks.bed"  # BED, SAF, または merge_SV.tsv ファイル
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

## 出力フォーマット

パイプラインは以下のファイルを生成します：

1. **logCPMテーブル** (デフォルト: `logCPM.tsv`):
   - 各ピーク/領域におけるログ変換済みCPM値
   - 列：PeakID、Chr、Start（0-based、BED形式基準）、End、各サンプルのlogCPM値

2. **中間ファイル**:
   - **flagstat出力** (`flagstat/` ディレクトリ):
     - 各BAMファイルのsamtools flagstat結果（総マップ済みリード数などの情報）
   - **featureCounts出力** (例: `annotation_base_fc_featureCounts.txt`):
     - 各領域における生のリードカウント
   - **実行ログ** (例: `pipeline_run.log`):
     - パイプライン実行の詳細ログ

例（logCPMテーブル）：
```
PeakID  Chr     Start   End     Sample1  Sample2  Sample3
peak1   chr1    999     2000    2.45     3.12     1.87
peak2   chr1    2999    4000    4.21     3.98     4.56
...
```

> **注意**: すべての出力ファイルの `Start` 列は0-based（BED形式基準）座標です。入力がSAF形式（1-based）の場合でも、出力は0-basedに変換されます。

## SEgeneワークフローでの位置づけ

このパイプラインは、SEgeneプロジェクトにおいて以下のようなワークフローの一部として機能します：

1. **シーケンスデータの前処理**：標準的なマッピングパイプラインを使用してBAMファイルを生成
2. **CPM/logCPM値の計算**（本パイプライン）
3. **peak-to-gene linksの構築**（`peak_to_gene_links`モジュール）
   - 入力：本パイプラインで生成したlogCPMファイル
   - 入力：RNA-seq TPMデータ
4. **スーパーエンハンサーの同定と解析**（`SE_to_gene_links`）
5. **オプショナル：公共データに基づく領域評価**（`SEgene_RegionAnalyzer`）

## トラブルシューティング

よくある問題と解決策：

1. **samtools/featureCountsが見つからない**：
   - PATHに追加されているか確認、または完全なパスを`--samtools_path`と`--featurecounts_path`で指定

2. **BAMファイルが見つからない**：
   - `--filename_pattern`引数を確認

3. **サンプル名抽出の問題**：
   - `--sample_delimiter`の値を確認
   - ハイフンを含むデリミタは`--delimiter="値"`の形式で指定

4. **アノテーションファイル形式の問題**：
   - BEDファイルを使用する場合は拡張子が`.bed`であることを確認
   - merge_SV形式を使用する場合は`--is_mergesv_format`フラグを指定
   - SAFファイルのフォーマットが正しいか確認（タブ区切り）

5. **featureCountsがエラーを返す**：
   - BAMファイルがインデックス付けされているか確認
   - ペアエンド/シングルエンドの設定が正しいか確認（`--single_end`フラグ）

## 引用

このツールを研究に使用する場合は、以下を引用してください：
**(論文は現在準備中です)**

## ライセンス

このプログラムはMITライセンスで公開されています。詳細については、[LICENSE](./LICENSE)をご覧ください。