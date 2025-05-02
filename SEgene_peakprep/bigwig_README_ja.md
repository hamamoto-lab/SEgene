# SEgene_peakprep - BigWig方式の詳細ガイド

このドキュメントでは、**SEgene_peakprep**の**BigWig方式**に関する詳細な使用手順と例を提供します。基本的な概要や共通要件については、[メインREADME](./README_ja.md)を参照してください。

## パイプライン構造と処理フロー

BigWig方式は、処理ニーズに応じて使用できる2つの主要なパイプラインで構成されています：

### 1. 単一サンプル解析パイプライン

**メインスクリプト**: `bigwig_peakprep.py`

このパイプラインは、個々のBAMファイルから正規化されたbigWigファイルを生成し、特定の領域（ピークなど）におけるシグナル値を定量化します。内部的には、以下の手順で処理が行われます：

1. `bigwig_peakprep_bamconvert.py`: BAMファイルを正規化されたbigWigファイルに変換
   - bamCoverageを使用して各サンプルのbigWigファイルを生成
   - 出力: bigWigファイルと変換詳細ファイル (`bam_to_bigwig_details.tsv`)

2. `bigwig_peakprep_summary.py`: 生成されたbigWigファイルから指定領域のシグナル値を集計
   - multiBigwigSummaryを使用して各領域からシグナル値を抽出
   - 出力: 生カウントテーブルとlog2変換テーブル

**使用例**:
```bash
python bigwig_peakprep.py \
    --bam_dir="/path/to/bam_files" \
    --annotation_file="/path/to/peaks.bed" \
    --output_dir="/path/to/output"
```

### 2. 差分解析パイプライン

**メインスクリプト**: `bigwig_peakprep_bamdiff.py`

このパイプラインは、ChIP-seqサンプルとInputコントロール間の比較（主にlog2比）を実行し、バックグラウンド補正されたシグナル値を計算します。内部的には、以下の手順で処理が行われます：

1. `bigwig_peakprep_bamdiff_generatediff.py`: サンプルとコントロールのペアリングを行い、差分比較を実行
   - bamCompareを使用して各サンプル/コントロールペアのlog2比bigWigファイルを生成
   - 出力: 差分bigWigファイルと変換詳細ファイル (`bamdiff_to_bigwig_details.tsv`)

2. `bigwig_peakprep_summary.py`: 差分bigWigファイルからシグナル値を集計
   - 注: 差分bigWigファイルは既にlog2比となっているため、追加のlog2変換は行われません

**使用例**:
```bash
python bigwig_peakprep_bamdiff.py \
    --sample_bam_dir="/path/to/chipseq_bams" \
    --control_bam_dir="/path/to/input_bams" \
    --annotation_file="/path/to/peaks.bed" \
    --output_dir="/path/to/output"
```

### ユーティリティモジュール

- `bigwig_peakprep_utils.py`: 基本関数と共通ユーティリティを提供
  - BAMからbigWigへの変換（bamCoverage実行）
  - PyRangesを使用したBEDファイル処理
  - multiBigwigSummaryの実行と結果処理
  - log2変換と各種データフォーマット変換機能

- `bigwig_peakprep_bamdiff_utils.py`: 差分解析に特化した関数を提供
  - サンプル-コントロールペアリング
  - bamCompare実行用関数
  - YAML変換プランの管理
  - 結果ファイル生成機能

### 使用するスクリプトの選択

- **単一サンプル解析**: `bigwig_peakprep.py` を使用
  - 各BAMファイルを個別に解析し、正規化されたシグナル値を取得したい場合

- **ChIP-seq/Input差分解析**: `bigwig_peakprep_bamdiff.py` を使用
  - バックグラウンド補正されたシグナル値（log2比）を取得したい場合
  - サンプルとコントロールのペアがある場合に最適

## 出力ディレクトリ構造

典型的な出力ディレクトリ構造は以下のようになります：

### 単一サンプル解析（bigwig_peakprep.py）の場合

```
output_dir/
├── bigwig/                              # bigWigファイル保存ディレクトリ
│   ├── sample1.RPGC.bigwig              # 正規化されたbigWigファイル
│   ├── sample2.RPGC.bigwig
│   ├── ...
│   └── sampleN.RPGC.bigwig
├── summary/                             # multiBigwigSummary出力ディレクトリ
│   ├── bigwig_summary.npz               # バイナリデータ（NumPy形式）
│   ├── bigwig_summary.tab               # タブ区切りテキストデータ
│   └── bigwig_summary.log               # 実行ログ
├── conversion_logs/                     # 変換ログディレクトリ（該当する場合）
│   └── peaks_converted.bed              # 変換されたBEDファイル
├── bam_to_bigwig_details.tsv            # BAM→bigWig変換の詳細情報
├── multibigwig_counts.tsv               # 生カウントテーブル
├── multibigwig_counts_log2.tsv          # log2変換カウントテーブル
└── pipeline_run.log                     # パイプライン実行ログ
```

### 差分解析（bigwig_peakprep_bamdiff.py）の場合

```
output_dir/
├── bigwig/                              # bigWigファイル保存ディレクトリ
│   ├── sample1_vs_input1.log2ratio.bigwig  # log2比bigWigファイル
│   ├── sample2_vs_input2.log2ratio.bigwig
│   ├── ...
│   └── sampleN_vs_inputN.log2ratio.bigwig
├── summary/                             # multiBigwigSummary出力ディレクトリ
│   ├── bigwig_summary.npz               # バイナリデータ（NumPy形式）
│   ├── bigwig_summary.tab               # タブ区切りテキストデータ
│   └── bigwig_summary.log               # 実行ログ
├── conversion_logs/                     # 変換ログディレクトリ（該当する場合）
│   └── peaks_converted.bed              # 変換されたBEDファイル
├── conversion_plan.yaml                 # 変換計画YAML
├── conversion_plan.tsv                  # 変換計画TSV形式
├── bamdiff_to_bigwig_details.tsv        # 差分bigWig生成の詳細情報
├── bamdiff_to_bigwig_details.yaml       # 詳細情報YAML形式
├── multibigwig_counts.tsv               # 差分値カウントテーブル
└── pipeline_run.log                     # パイプライン実行ログ
```

## 引数一覧

ここでは、コマンドライン実行に使用できる全ての引数を、2つのメインスクリプトごとに分けて説明します。

### 単一サンプル解析（bigwig_peakprep.py）の引数

| 引数 | 説明 | デフォルト値 |
|--------|------|-------------|
| **必須引数** |||
| `--bam_dir` | BAMファイルが格納されているディレクトリ（BAM→bigWig変換に必要） | - |
| `--annotation_file` | アノテーションファイル（BED、SAF、またはmerge_SE形式） | - |
| `--output_dir` | 結果を保存するディレクトリ | - |
| `--bigwig_details_file` | bam_to_bigwig_details.tsvファイルへのパス、サマリーのみ実行時（`--run_summary_only`）に必要 | None |
| **ツールパス関連** |||
| `--tools_dir` | deeptoolsの実行ファイルを含むディレクトリ | None |
| `--bamcoverage_path` | bamCoverageの実行ファイルへのパス（`--tools_dir`指定時は無視） | "bamCoverage" |
| `--multibigwigsummary_path` | multiBigwigSummaryの実行ファイルへのパス（`--tools_dir`指定時は無視） | "multiBigwigSummary" |
| **ファイル/サンプル選択** |||
| `--filename_pattern` | BAMファイル名のワイルドカードパターン | "*.bam" |
| `--sample_delimiter` | BAMファイル名からサンプル名を抽出するための区切り文字列 | None |
| **bamCoverage関連** |||
| `--bigwig_dir` | bigWigファイルを保存するディレクトリ | output_dir/bigwig |
| `--bin_size` | bamCoverageのビンサイズ（ベース単位） | 50 |
| `--normalize_using` | bamCoverageの正規化方法（RPGC、CPMなど） | "RPGC" |
| `--effective_genome_size` | RPGC正規化の有効ゲノムサイズ | 2913022398 |
| **multiBigwigSummary関連** |||
| `--summary_dir` | multiBigwigSummaryファイルを保存するディレクトリ | output_dir/summary |
| `--summary_basename` | multiBigwigSummary出力ファイルのベース名 | "bigwig_summary" |
| **データテーブル関連** |||
| `--force_chr_start_end_ids` | BEDファイルに名前列があってもchr_start_end形式をPeakIDとして使用 | False |
| `--pseudocount` | log2変換の擬似カウント | 1.0 |
| `--raw_counts_name` | 生カウントテーブルのファイル名 | "multibigwig_counts.tsv" |
| `--log2_counts_name` | log2変換カウントテーブルのファイル名 | "multibigwig_counts_log2.tsv" |
| **実行制御** |||
| `--threads` / `-T` | deeptoolsコマンドのスレッド数 | 4 |
| `--force_overwrite` | 既存のbigWigファイルを強制的に上書き（デフォルトではスキップ） | False |
| `--run_bamconvert_only` | BAM→bigWig変換のみ実行 | False |
| `--run_summary_only` | bigWigサマリー処理のみ実行 | False |

### 差分解析（bigwig_peakprep_bamdiff.py）の引数

| 引数 | 説明 | デフォルト値 |
|--------|------|-------------|
| **必須引数** |||
| `--sample_bam_dir` | サンプルBAMファイルが格納されているディレクトリ | - |
| `--control_bam_dir` | コントロールBAMファイル（Input/IgG）が格納されているディレクトリ | - |
| `--annotation_file` | アノテーションファイル（BED、SAF、またはmerge_SE形式） | - |
| `--output_dir` | 結果を保存するディレクトリ | - |
| **ツールパス関連** |||
| `--tools_dir` | deeptoolsの実行ファイルを含むディレクトリ | None |
| `--bamcompare_path` | bamCompareの実行ファイルへのパス | "bamCompare" |
| **ファイル選択** |||
| `--sample_pattern` | サンプルBAMファイルのグロブパターン | "*.bam" |
| `--control_pattern` | コントロールBAMファイルのグロブパターン | "*Input*.bam" |
| `--sample_delimiter` | サンプル名抽出のための区切り文字 | "_" |
| **差分解析関連** |||
| `--operation` | bamCompareの操作タイプ（log2、ratio、subtract、add、meanなど） | "log2" |
| `--bin_size` | bamCompareのビンサイズ（ベース単位） | 50 |
| `--pseudocount` | ゼロ除算を避けるための擬似カウント値 | 1.0 |
| `--effective_genome_size` | 有効ゲノムサイズ | 2913022398 |
| **実行制御** |||
| `--threads` / `-T` | スレッド数 | 4 |
| `--force_overwrite` | 既存のbigWigファイルを強制上書き | False |
| `--run_diff_only` | 差分bigWig生成のみ実行（サマリー処理をスキップ） | False |

## 単一サンプル解析パイプライン（bigwig_peakprep.py）の詳細な使用例

### 基本的な使用法

```bash
python bigwig_peakprep.py \
    --bam_dir="/path/to/bam_files" \
    --annotation_file="/path/to/peaks.bed" \
    --output_dir="/path/to/output" \
    --tools_dir="/path/to/deeptools/bin" \
    --sample_delimiter="-H3K27ac.sorted.bam" \
    --normalize_using="RPGC" \
    --bin_size=50 \
    --threads=10
```

### 高度なデータテーブル設定

```bash
python bigwig_peakprep.py \
    --bam_dir="/path/to/bam_files" \
    --annotation_file="/path/to/peaks.bed" \
    --output_dir="/path/to/output" \
    --tools_dir="/path/to/deeptools/bin" \
    --force_chr_start_end_ids \
    --pseudocount=0.1 \
    --raw_counts_name="my_raw_counts.tsv" \
    --log2_counts_name="my_log2_counts.tsv" \
    --threads=10
```

### 2段階処理: BAM変換のみ

```bash
python bigwig_peakprep.py \
    --bam_dir="/path/to/bam_files" \
    --annotation_file="/path/to/peaks.bed" \
    --output_dir="/path/to/output" \
    --tools_dir="/path/to/deeptools/bin" \
    --normalize_using="CPM" \
    --threads=10 \
    --run_bamconvert_only
```

### 2段階処理: サマリー処理のみ

```bash
python bigwig_peakprep.py \
    --annotation_file="/path/to/peaks.bed" \
    --output_dir="/path/to/output" \
    --tools_dir="/path/to/deeptools/bin" \
    --bigwig_details_file="/path/to/output/bam_to_bigwig_details.tsv" \
    --threads=10 \
    --run_summary_only
```

### 個別ステップスクリプトを使用した直接実行

個別スクリプトは、より細かな制御が必要な場合や特定の処理ステップのみを実行したい場合に便利です。各スクリプトには独自のオプションがあります。`--help`フラグを付けて各スクリプトを実行すると、詳細なオプションを確認できます。

BAM→bigWig変換（ステップ1）:
```bash
python bigwig_peakprep_bamconvert.py \
    --bam_dir="/path/to/bam_files" \
    --output_dir="/path/to/output" \
    --tools_dir="/path/to/deeptools/bin" \
    --sample_delimiter="-H3K27ac.sorted.bam" \
    --threads=10
```

bigWigサマリー処理（ステップ2）:
```bash
python bigwig_peakprep_summary.py \
    --bigwig_details_file="/path/to/output/bam_to_bigwig_details.tsv" \
    --annotation_file="/path/to/peaks.bed" \
    --output_dir="/path/to/output" \
    --tools_dir="/path/to/deeptools/bin" \
    --threads=10
```

## 差分解析パイプライン（bigwig_peakprep_bamdiff.py）の詳細な使用例

差分解析パイプラインは、ChIP-seqサンプルとInputコントロールを比較し、バックグラウンドシグナルを補正して結合部位の強度を計算します。主にlog2比計算を実行し、より精度の高い結合部位解析を可能にします。

### 基本的な使用法

```bash
python bigwig_peakprep_bamdiff.py \
    --sample_bam_dir="/path/to/chipseq_bams" \
    --control_bam_dir="/path/to/input_bams" \
    --annotation_file="/path/to/peaks.bed" \
    --output_dir="/path/to/output" \
    --tools_dir="/path/to/deeptools/bin" \
    --threads=10
```

### サンプルとコントロールの自動マッチング

`bigwig_peakprep_bamdiff.py`は、ファイル名パターンに基づいてサンプル（ChIP-seq）とコントロール（Input）を自動的にマッチングします。例えば：

- サンプル: `T1-H3K27ac.bam`, `T2-H3K27ac.bam`
- コントロール: `T1-Input.bam`, `T2-Input.bam`

この場合、T1サンプルはT1コントロールと、T2サンプルはT2コントロールと自動的にペアリングされます。また、単一のコントロールを全てのサンプルに適用することも可能です。

### 特定のパターンを使用したサンプルとコントロールの指定

```bash
python bigwig_peakprep_bamdiff.py \
    --sample_bam_dir="/path/to/bams" \
    --control_bam_dir="/path/to/bams" \
    --sample_pattern="*-H3K27ac.bam" \
    --control_pattern="*-Input.bam" \
    --annotation_file="/path/to/peaks.bed" \
    --output_dir="/path/to/output" \
    --threads=10
```

### 差分bigWig生成のみ実行（サマリー処理をスキップ）

```bash
python bigwig_peakprep_bamdiff.py \
    --sample_bam_dir="/path/to/chipseq_bams" \
    --control_bam_dir="/path/to/input_bams" \
    --annotation_file="/path/to/peaks.bed" \
    --output_dir="/path/to/output" \
    --run_diff_only \
    --threads=10
```

### 操作タイプの変更（log2以外）

```bash
python bigwig_peakprep_bamdiff.py \
    --sample_bam_dir="/path/to/chipseq_bams" \
    --control_bam_dir="/path/to/input_bams" \
    --annotation_file="/path/to/peaks.bed" \
    --output_dir="/path/to/output" \
    --operation="subtract" \
    --threads=10
```

## 出力ファイルの詳細

### 単一サンプル解析パイプライン（bigwig_peakprep.py）の出力

1. **bigWigファイル** (`bigwig/`ディレクトリ内):
   - 各BAMファイルから生成された正規化bigWigファイル
   - ファイル名形式: `{sample_name}.{normalization_method}.bigwig`
   - 例: `sample1.RPGC.bigwig`

2. **bam_to_bigwig_details.tsv**:
   - ステップ1（ステップ2の入力）からの詳細出力ファイル
   - コメント行にはパラメータ情報と実行メタデータが含まれる
   - 列: `Sample_name`, `BAM_filename`, `BigWig_filename`, `BAM_fullpath`, `BigWig_fullpath`

3. **multiBigwigSummary出力** (`summary/`ディレクトリ内):
   - `{summary_basename}.npz`: NumPy形式のバイナリデータ
   - `{summary_basename}.tab`: タブ区切りテキストデータ

4. **カウントテーブル**:
   - `multibigwig_counts.tsv`: 生カウントテーブル（PeakID情報付き）
   - `multibigwig_counts_log2.tsv`: log2変換カウントテーブル（PeakID情報付き）

### 差分解析パイプライン（bigwig_peakprep_bamdiff.py）の出力

1. **差分bigWigファイル** (`bigwig/`ディレクトリ内):
   - サンプルとコントロール間のlog2比を表すbigWigファイル
   - ファイル名形式: `{sample_name}_vs_{control_name}.log2ratio.bigwig`

2. **bamdiff_to_bigwig_details.tsv**:
   - 差分解析の詳細情報を含むTSVファイル
   - 列: `Sample_name`, `BigWig_filename`, `BigWig_fullpath`
   - 重要: このファイルにはメタデータ`already_log2_transformed: true`が含まれる（サマリー処理中の追加log2変換を防ぐため）

3. **カウントテーブル** (サマリー処理実行時):
   - `multibigwig_counts.tsv`: 各ピーク領域のlog2比値を含むテーブル
   - 注: 値は既にlog2変換されているため、追加の変換は実行されません（log2_counts_nameオプションで指定されたファイルは生成されない場合があります）

カウントテーブルの例:
```
PeakID  Chr     Start   End     Sample1  Sample2  Sample3
peak1   chr1    999     2000    2.45     3.12     1.87
peak2   chr1    2999    4000    4.21     3.98     4.56
...
```

> **注**: 出力ファイルの`Start`列は0-based（BED形式標準）です。

## 高度な機能: YAML変換プラン

差分解析パイプラインは、YAMLフォーマットでサンプルとコントロールの変換プランを管理します。この機能により、複雑な解析設定を簡単に保存して再利用できます。

### YAML変換プランの構造

変換プランファイル（`conversion_plan.yaml`）には以下の情報が含まれます：

```yaml
meta_info:
  generated_date: "2023-01-01 00:00:00"
  pipeline: "bigwig_peakprep_bamdiff pipeline"
parameters:
  operation: "log2"
  bin_size: 50
  pseudocount: 1.0
  effective_genome_size: 2913022398
  threads: 4
samples:
  T1-H3K27ac:
    bam_file: "T1-H3K27ac.bam"
    bam_path: "/path/to/T1-H3K27ac.bam"
controls:
  T1-Input:
    bam_file: "T1-Input.bam"
    bam_path: "/path/to/T1-Input.bam"
conversion_plan:
  - pair_id: "T1-H3K27ac_vs_T1-Input"
    sample_name: "T1-H3K27ac"
    control_name: "T1-Input"
    sample_bam: "T1-H3K27ac.bam"
    control_bam: "T1-Input.bam"
    output_bigwig: "T1-H3K27ac_vs_T1-Input.log2ratio.bigwig"
    sample_bam_path: "/path/to/T1-H3K27ac.bam"
    control_bam_path: "/path/to/T1-Input.bam"
    output_bigwig_path: "/path/to/output/bigwig/T1-H3K27ac_vs_T1-Input.log2ratio.bigwig"
```

### YAML処理結果

処理結果もYAML形式（`bamdiff_to_bigwig_details.yaml`）で保存され、各ジョブの成功/失敗に関する詳細情報が含まれます：

```yaml
results:
  summary:
    total_jobs: 5
    completed_jobs: 4
    skipped_jobs: 1
    failed_jobs: 0
  job_results:
    - pair_id: "T1-H3K27ac_vs_T1-Input"
      status: "completed"
      start_time: "2023-01-01T00:00:00"
      end_time: "2023-01-01T00:05:00"
      duration_seconds: 300
      output_file: "/path/to/output/bigwig/T1-H3K27ac_vs_T1-Input.log2ratio.bigwig"
  successful_files: ["/path/to/output/bigwig/T1-H3K27ac_vs_T1-Input.log2ratio.bigwig", ...]
  failed_pairs: []
```

これらの変換プランと結果YAMLファイルは、解析設定の再現や後続の解析で特定のサンプルを再処理するのに役立ちます。

## 正規化方法の選択

### 単一サンプル解析（bigwig_peakprep.py）の正規化

`--normalize_using`オプションを使用して、以下の正規化方法から選択できます：

- **RPGC**（デフォルト）: Reads Per Genomic Content。指定されたゲノムサイズに基づいて正規化
- **CPM**: Counts Per Million mapped reads
- **BPM**: Bins Per Million mapped reads
- **RPKM**: Reads Per Kilobase per Million mapped reads
- **None**: 正規化なし

例えば、CPM正規化を使用するには：
```bash
python bigwig_peakprep.py \
    --bam_dir="/path/to/bam_files" \
    --annotation_file="/path/to/peaks.bed" \
    --output_dir="/path/to/output" \
    --normalize_using="CPM" \
    --threads=10
```

### 差分解析（bigwig_peakprep_bamdiff.py）の操作タイプ

`--operation`オプションを使用して、以下の操作タイプから選択できます：

- **log2**（デフォルト）: log2(sample/control)比（擬似カウントは`--pseudocount`で指定）
- **ratio**: sample/control比
- **subtract**: sample - control減算
- **add**: sample + control加算
- **mean**: (sample + control)/2 平均
- **reciprocal_ratio**: control/sample逆比
- **first**: サンプルのみ出力
- **second**: コントロールのみ出力

例えば、単純な減算を実行するには：
```bash
python bigwig_peakprep_bamdiff.py \
    --sample_bam_dir="/path/to/chipseq_bams" \
    --control_bam_dir="/path/to/input_bams" \
    --annotation_file="/path/to/peaks.bed" \
    --output_dir="/path/to/output" \
    --operation="subtract" \
    --threads=10
```

## トラブルシューティング

よくある問題と解決策：

1. **ハイフンを含む引数値のエラー**:
   - 引数に等号形式を使用： `--sample_delimiter="-H3K27ac.sorted.bam"`
   - スペース区切り形式を使用しないこと： ~~`--sample_delimiter "-H3K27ac.sorted.bam"`~~

2. **bamCoverage/multiBigwigSummaryが見つからないエラー**:
   - `--tools_dir`オプションで正しいdeeptoolsのbinディレクトリを指定
   - または、`--bamcoverage_path`と`--multibigwigsummary_path`で完全なパスを指定

3. **サンプル-コントロール自動マッチングの問題**:
   - サンプルとコントロールのファイル名に共通の識別子（T1、T2など）があることを確認
   - 必要に応じて`--sample_pattern`と`--control_pattern`を使用して明示的にパターンを指定

4. **メモリ不足エラー**:
   - 多数のファイルや大きなゲノム領域を処理する場合は、より多くのメモリを持つマシンで実行
   - ファイル数を減らすか、バッチで処理することを検討

5. **処理速度が遅い**:
   - `--threads`オプションでスレッド数を増やす
   - インデックス付きのBAMファイルを使用
   - `--bin_size`値を大きくすると処理が速くなりますが、解像度は低下します

### 差分解析特有の問題

6. **サンプル-コントロールマッチングの失敗**:
   - ログに「Could not create sample-control pairs」警告メッセージが表示される場合、マッチングに失敗しています
   - `--sample_pattern`と`--control_pattern`を使用してファイル名パターンを明示的に指定
   - ファイル名に共通の識別子（T1、T2など）があることを確認
   - 全てのサンプルに単一のコントロールを使用する場合は、コントロールディレクトリに1つのファイルのみが含まれていることを確認

7. **差分bigWig生成後のエラー**:
   - 差分bigWigファイルは存在するがサマリー処理に失敗する場合、バイナリ形式が破損している可能性があります
   - `--force_overwrite`オプションを使用して差分bigWigファイルを再生成
   - 生成に失敗したファイルのログ（`bamCompare.log`）を確認

8. **正しいスクリプトが使用されているかの確認**:
   - サンプルとコントロールの差分解析には、`bigwig_peakprep.py`ではなく`bigwig_peakprep_bamdiff.py`を使用していることを確認
   - 逆に、単一サンプル解析には`bigwig_peakprep.py`を使用

## 引用

この研究ツールを使用する場合は、以下を引用してください：
**(論文は現在準備中です)**

## ライセンス

このプログラムはMITライセンスで公開されています。詳細については、[LICENSE](https://github.com/hamamoto-lab/SEgene/blob/main/LICENSE)をご覧ください。