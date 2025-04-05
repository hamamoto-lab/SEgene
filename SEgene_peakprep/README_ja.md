# SEgene_peakprep

*(For the English version of this README, please see [README.md](./README.md).)*

**SEgene_peakprep** は、ChIP-seqデータ（BAMファイル）から指定されたゲノム領域における読み取りカウントを正規化し、CPM（Counts Per Million）値を計算するためのPythonパイプラインです。このツールは **hg38参照ゲノム** にマッピングされたデータを対象に設計されています。

## 開発状況

> **注意**: このツールは**現在開発中**であり、プログラムのファイル構成や使用方法はバージョンによって大きく変化する可能性があります。

## 概要

SEgene_peakprepは、以下の機能を提供します：

- 複数のBAMファイルの効率的な処理とサンプル情報の自動抽出
- `samtools flagstat`を使用した総マップ済みリード数の計算
- `featureCounts`を使用した指定ゲノム領域（SAFフォーマット）におけるリードカウント
- CPM値の計算とオプションのログ変換（log2、log10）
- 結果の整形と保存

**SEgeneプロジェクト内での位置づけ：**

このパイプラインは[SEgeneプロジェクト](https://github.com/hamamoto-lab/SEgene)において、初期データ準備の役割を果たします：

1. **入力**: マッピング済みChIP-seqデータ（BAMファイル）と関心領域（SAFフォーマット）
2. **出力**: CPM/logCPM値テーブル（TSVファイル）
3. **フロー**: この出力は、`peak_to_gene_links`モジュールの主要な入力となり、RNA-seqのTPMデータと組み合わせて遺伝子制御ネットワーク解析に使用されます

## 要件

- **Python**: Python 3.7以上
- **必須ライブラリ**:
  - pandas
  - numpy
  - natsort
  
  ※必要なライブラリはcondaまたはpipを使用してインストールしてください
  
- **外部ツール**:
  - samtools (flagstat機能用) - PATHを通すか、実行時にオプションでパスを指定
  - featureCounts (Subreadパッケージの一部) - PATHを通すか、実行時にオプションでパスを指定

## インストール

本パイプラインを使用するには、必要なPythonライブラリ（pandas, numpy, natsort）をcondaまたはpipを使用してインストールしてください。

このリポジトリをクローンするか、スクリプトをダウンロードします：

```bash
git clone https://github.com/hamamoto-lab/SEgene.git
cd SEgene/SEgene_peakprep
```

## 入力データの準備

### SAFファイルの準備
このパイプラインは入力として、Simplified Annotation Format (SAF) ファイルを使用します。SAFファイルは以下の列を含むタブ区切りのテキストファイルです：

1. `GeneID` - 各領域の一意な識別子
2. `Chr` - 染色体名
3. `Start` - 開始位置（1-based）
4. `End` - 終了位置
5. `Strand` - ストランド（+、-、または.）

例：
```
GeneID  Chr     Start   End     Strand
peak1   chr1    1000    2000    .
peak2   chr1    3000    4000    .
```

BEDファイルからSAFファイルへの変換が必要な場合は、以下のようなスクリプトを使用できます：

```python
import pandas as pd

# BEDファイルを読み込み
bed_df = pd.read_csv("your_file.bed", sep='\t', header=None)

# SAF形式に変換
saf_df = pd.DataFrame({
    'GeneID': [f"peak_{i+1}" for i in range(len(bed_df))],
    'Chr': bed_df[0],  # 染色体列
    'Start': bed_df[1] + 1,  # BEDは0-basedなので1を加える
    'End': bed_df[2],
    'Strand': '.'  # ストランド情報がなければ '.' を使用
})

# SAFファイルとして保存
saf_df.to_csv("output.saf", sep='\t', index=False)
```

### BAMファイルの準備
BAMファイルは標準的なChIP-seqマッピングパイプライン（nf-core/chipseq など）を使用して生成してください。BAMファイルはインデックス付け（.baiファイル）されていることを推奨します。

## 使用方法

パイプラインは、コマンドラインまたはJupyterノートブックから実行できます。

### コマンドライン実行

基本的な実行コマンド：

```bash
python cpm_peakprep_pipeline.py \
    --bam_dir /path/to/bam_files \
    --saf_file peaks.saf \
    --output_dir results \
    --filename_pattern "Sample*.sorted.bam" \
    --sample_delimiter ".sorted.bam" \
    --threads 10
```

#### 引数の説明

| 引数 | 説明 |
|------|------|
| `--bam_dir` | BAMファイルが格納されているディレクトリのパス（**必須**） |
| `--saf_file` | 解析対象のゲノム領域を定義したSAFファイルのパス（**必須**） |
| `--output_dir` | 結果を保存するディレクトリ（**必須**） |
| `--filename_pattern` | BAMファイル名をフィルタリングするためのワイルドカードパターン（例：`*.bam`、`T*H3K27ac*.bam`）（オプション、デフォルト：`*.bam`） |
| `--sample_delimiter` | BAMファイル名からサンプル名を抽出するための区切り文字列（オプション）<br>例：ファイル名が`T9-H3K27ac_R1.mLb.clN.sorted.bam`で、区切り文字列が`-H3K27ac_R1.mLb.clN.sorted.bam`なら、サンプル名は`T9`となる |
| `--samtools_path` | samtoolsの実行ファイルへのパス（オプション、デフォルト：`samtools`） |
| `--featurecounts_path` | featureCountsの実行ファイルへのパス（オプション、デフォルト：`featureCounts`） |
| `--threads` | samtoolsとfeatureCountsで使用するスレッド数（オプション、デフォルト：4） |
| `--single_end` | シングルエンドのリードを指定するフラグ。指定しない場合はペアエンドと見なされる（オプション、デフォルト：ペアエンド） |
| `--fc_output_name` | featureCountsの出力ファイル名（オプション、デフォルト：SAFファイル名から自動生成） |
| `--fc_log_name` | featureCountsのログファイル名（オプション） |
| `--fc_options` | featureCountsに渡す追加オプション（オプション） |
| `--add_peakid` | オリジナルのGeneIDを連番のPeakID（PeakID_XXXX）に置き換えるフラグ（オプション） |
| `--id_column` | featureCounts出力のID列名（オプション、デフォルト：`Geneid`） |
| `--output_id_name` | 最終的なlogCPMテーブルのID列名（オプション、デフォルト：`PeakID`） |
| `--log_base` | 対数変換の底（2または10）。0以下を設定すると変換なし（オプション、デフォルト：2） |
| `--pseudocount` | 対数変換前に追加する擬似カウント（オプション、デフォルト：1.0） |
| `--logcpm_output_name` | 最終的なlogCPMテーブルのファイル名（オプション、デフォルト：`logCPM_table.tsv`） |
| `--log_level` | ロギングレベル（DEBUG/INFO/WARNING/ERROR/CRITICAL）（オプション、デフォルト：INFO） |

### Jupyterノートブックでの実行

以下は、Jupyterノートブックでパイプラインを実行する基本的な例です：

```python
import os
import logging
import pandas as pd
# ユーティリティ関数をインポート
from cpm_peakprep_utils import (
    get_sample_data,
    run_samtools_flagstat,
    run_featurecounts,
    calculate_logcpm
)

# パラメータ設定
bam_dir = "/path/to/bam_files"
saf_file = "peaks.saf"
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

# ステップ1: BAMファイル情報の取得
sample_dict, bam_list = get_sample_data(
    bam_folder=bam_dir,
    logger=logger,
    sample_name_delimiter=sample_delimiter,
    filename_pattern=filename_pattern
)

# ステップ2: 総リード数の取得
flagstat_dir = os.path.join(output_dir, "flagstat_output")
total_reads_map = run_samtools_flagstat(
    bam_files=bam_list,
    output_dir=flagstat_dir,
    logger=logger,
    threads=threads
)

# ステップ3: リードカウント（featureCounts）
fc_output_path = os.path.join(output_dir, f"{os.path.splitext(os.path.basename(saf_file))[0]}_featureCounts.txt")
fc_log_path = os.path.join(output_dir, f"{os.path.splitext(os.path.basename(saf_file))[0]}_featureCounts.log")
counts_file_path = run_featurecounts(
    bam_files=bam_list,
    saf_file=saf_file,
    logger=logger,
    output_file=fc_output_path,
    log_file=fc_log_path,
    threads=threads,
    is_paired_end=True
)

# ステップ4: LogCPM計算
bam_to_sample_map = {v: k for k, v in sample_dict.items()}
logcpm_df = calculate_logcpm(
    counts_file=counts_file_path,
    total_reads_dict=total_reads_map,
    logger=logger,
    bam_path_to_sample_name=bam_to_sample_map,
    log_transform_base=2,
    pseudocount=1.0
)

# 結果の保存
logcpm_output_file = os.path.join(output_dir, "logCPM_table.tsv")
logcpm_df.to_csv(logcpm_output_file, sep='\t', index=False, float_format='%.4f')
print(f"最終的なlogCPMテーブルが保存されました: {logcpm_output_file}")
```

## 出力フォーマット

パイプラインは以下のファイルを生成します：

1. **flagstat出力** (`flagstat_output/` ディレクトリ)：
   - 各BAMファイルのsamtools flagstat結果（総マップ済みリード数などの情報）

2. **featureCounts出力** (例: `peaks_featureCounts.txt`)：
   - 各領域における生のリードカウント
   - 列：GeneID、Chr、Start、End、Strand、Length、各BAMファイルのカウント

3. **logCPMテーブル** (例: `logCPM_table.tsv`)：
   - 最終的なログ変換済みCPM値
   - 列：PeakID、Chr、Start、End、各サンプルのlogCPM値

例（logCPMテーブル）：
```
PeakID  Chr     Start   End     T1      T2      T3
peak1   chr1    1000    2000    2.45    3.12    1.87
peak2   chr1    3000    4000    4.21    3.98    4.56
...
```

## パイプラインのワークフロー

1. **BAMファイル情報の取得**：
   - 指定されたディレクトリから条件に合うBAMファイルを検索
   - ファイル名からサンプル名を抽出

2. **総マップ済みリード数の取得**：
   - `samtools flagstat`を使用して各BAMファイルの総マップ済みリード数を計算
   - 結果はCPM正規化に使用

3. **領域ごとのリード数のカウント**：
   - `featureCounts`を使用して各BAMファイルから指定された領域（SAFファイル）内のリードをカウント
   - ペアエンドまたはシングルエンドリードに対応

4. **CPM計算とログ変換**：
   - 各リードカウントを総マップ済みリード数で正規化（CPM）
   - オプションでログ変換（log2またはlog10）を適用
   - サンプル名をわかりやすい形式に変換

5. **結果の出力**：
   - 最終的なテーブルをTSV形式で保存

## SEgeneワークフローでの位置づけ

このパイプラインは、SEgeneプロジェクトにおいて以下のようなワークフローの一部として機能します：

1. **ChIP-seqデータの前処理**：標準的なマッピングパイプライン（nf-core/chipseq など）を使用してBAMファイルを生成
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
   - `--filename_pattern`引数を確認、パターンが適切か確認

3. **featureCountsがエラーを返す**：
   - BAMファイルがインデックス付けされているか確認
   - ペアエンド/シングルエンドの設定が正しいか確認（`--single_end`フラグ）
   - SAFファイルのフォーマットが正しいか確認

## ライセンス

このプログラムはMITライセンスで公開されています。詳細については、[LICENSE](./LICENSE)をご覧ください。

## 引用

このツールを研究に使用する場合は、以下を引用してください：
**(論文は現在準備中です)**