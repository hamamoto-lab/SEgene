# SEgene_peakprep

*(For the English version of this README, please see [README.md](./README.md).)*

**SEgene_peakprep** は、シーケンスデータ（BAMファイル）から指定されたゲノム領域におけるシグナル値を定量化・正規化し、下流のSEgene解析で使用する形式に整形するためのPythonパイプラインです。このツールを使用して、**(i) BigWig**、**(ii) Standard log2-CPM**、**(iii) edgeR正規化CPM（calcnorm CPM）** の3種類の定量テーブルを生成できます。`--calcnorm`オプションを使用するとedgeR正規化CPMを追加出力します。

## 開発状況

> **注意:** このツールは**現在開発中**であり、プログラムのファイル構成や使用方法はバージョンによって変化する可能性があります。

## 概要

SEgene_peakprepは、[SEgeneプロジェクト](https://github.com/hamamoto-lab/SEgene)において、初期データ準備の役割を果たします：

1. **入力**: マッピング済みシーケンスデータ（BAMファイル）とアノテーションファイル（BEDファイルなど）
2. **出力**: 正規化されたカウント値テーブル
3. **フロー**: この出力は、`peak_to_gene_links`モジュールの主要な入力となり、RNA-seqのTPMデータと組み合わせて遺伝子制御ネットワーク解析に使用されます

このパイプラインは以下の2つの実装方法を提供します：

1. **CPM方式** (`cpm_peakprep.py`): `featureCounts`を使用して直接Counts Per Million（CPM）値を計算
2. **BigWig方式** (`bigwig_peakprep.py`): `deeptools`を使用してBAMからbigWigファイルに変換し、multiBigwigSummaryでピークを定量化

どちらの方法も、BED、SAF、merge_SE形式のアノテーションファイルに対応しています。用途に応じて適切な方法を選択できます  
**注意**: merge_SE形式はSEgeneプロジェクト特有の形式で、パイプライン内の出力および中間ファイルとして使用されます

## CPM方式における2モード

| モード | 切替オプション | 正規化方法 | 既定ファイル |
|--------|---------------|------------|-------------|
| Standard log2-CPM | （指定なし） | ライブラリ総リード数でCPM → log₂変換 | `logCPM.tsv` |
| edgeR正規化CPM<br>(calcnorm CPM) | `--calcnorm` | edgeR `calcNormFactors()`(UQ/TMM/RLE) → CPM | `calcnorm.tsv` |

`--calcnorm`を付けると**両方のテーブルが同時に出力**されます。  
ファイル名は`--logcpm_output_name`、`--calcnorm_output_name`で変更可能です。

## 共通要件・依存関係

以下の共通要件を満たしていることを前提とします。

### Python環境
- Python 3.10 以上

### Pythonライブラリ（pipまたはcondaでインストール）
| ライブラリ | 用途 | 推奨バージョン | 必要な機能 |
|------------|------|---------------|------------|
| pandas | テーブル操作 | ≥1.5 | すべての機能 |
| numpy | 数値計算 | ≥1.23 | すべての機能 |
| natsort | 自然順ソート | ≥8.3 | すべての機能 |
| rpy2 | R言語との連携 | ≥3.5.0 | edgeR正規化CPMモード |
| pyranges | ゲノム領域の操作 | ≥0.14 | BigWig方式のみ |

```bash
# pipでのインストール例
pip install pandas numpy natsort
pip install rpy2  # edgeR正規化CPMモードで必要
pip install pyranges  # BigWig方式で必要
```

```bash
# condaでのインストール例
conda install pandas numpy natsort
conda install -c conda-forge rpy2  # edgeR正規化CPMモードで必要
conda install -c bioconda pyranges  # BigWig方式で必要
```

### 外部ツール（condaやバイナリでインストール）
| ツール名 | 用途 | 推奨バージョン | 必要な機能 |
|----------|------|---------------|------------|
| samtools | BAM操作 | ≥1.13 | CPM方式 |
| featureCounts | CPMカウント | ≥2.0.0 | CPM方式 |
| R | 統計解析環境 | ≥4.2.0 | edgeR正規化CPMモード |
| edgeR | RNA-seq解析用Bioconductorパッケージ | ≥3.40.0 | edgeR正規化CPMモード |
| deeptools | BigWig生成・差分解析 | ≥3.5.0 | BigWig方式 |
| - bamCoverage | BAMからbigWig変換 | (deeptoolsに含む) | BigWig単一サンプル解析 |
| - multiBigwigSummary | bigWigからのシグナル集計 | (deeptoolsに含む) | BigWig単一サンプル解析 |
| - bamCompare | サンプル/コントロール比較 | (deeptoolsに含む) | BigWig差分解析(bamdiff) |

```bash
# condaでのインストール例
# CPM方式の基本ツール
conda install -c bioconda samtools subread

# edgeR正規化CPMモードのためのR/edgeR
conda install -c conda-forge r-base=4.2
conda install -c conda-forge r-essentials
conda install -c bioconda bioconductor-edger=3.40

# BigWig方式
conda install -c bioconda deeptools
```

### 共通の入力ファイル

両方の実装方法で、以下の入力ファイルが必要です：

1. **BAMファイル**: マッピング済みのシーケンスデータファイル（推奨：インデックス付き`.bai`ファイル）
2. **アノテーションファイル**: ゲノム領域を定義するファイルで、以下の形式に対応:
   - **BED形式**: 標準的なBEDフォーマット（0-based開始座標）
     ```
     # 例: peaks.bed（タブ区切り）
     chr1    1000    2000    peak1    0    +
     chr1    3000    4000    peak2    0    +
     chr2    5000    6000    peak3    0    -
     ```
   
   - **SAF形式** (Simplified Annotation Format): タブ区切りのテキストファイル（**1-based開始座標**）
     ```
     # 例: peaks.saf（タブ区切り）
     GeneID  Chr     Start   End     Strand
     peak1   chr1    1001    2000    +
     peak2   chr1    3001    4000    +
     peak3   chr2    5001    6000    -
     ```
     **注意:** SAF形式の座標は1-based（Start=1001は0-basedでは1000）ですが、出力は全て0-basedに統一されます
   
   - **merge_SE.tsv形式**: 1列目が `chr_start_end` 形式の領域ID（例：`chr1_1000_2000`）
     ```
     # 例: merged_se.tsv（タブ区切り）
     se_data            feature
     chr1_1000_2000     feature1
     chr1_3000_4000     feature2
     chr2_5000_6000     feature3
     ```
     **注意:** 
      - merge_SE形式では座標は0-basedと解釈されます
      - merge_SE形式はSEgeneプロジェクト特有の形式で、SEgeneパイプラインの出力および中間ファイルとして使用されます
      - 座標は0-based（chr1_1000_2000のStartは0-basedで1000）と解釈されます
      - このフォーマットを使うと、SEgeneの他のコンポーネントとシームレスに連携できます

## CPM方式（cpm_peakprep）

CPM方式は、`featureCounts`を使用して指定ゲノム領域（BED/SAF/merge_SE.tsv形式）におけるリードを直接カウントし、以下の2つのモードで正規化値を計算します：

1. **Standard log2-CPM**（デフォルト）: 総リード数に基づくシンプルなCPM（Counts Per Million）計算とlog2変換
2. **edgeR正規化CPM**（以下、**calcnorm CPM**と略）: edgeRパッケージのスケーリングファクターでライブラリサイズ補正を行ったCPM

モードは`--calcnorm`オプションで切り替えます：

```bash
python cpm_peakprep.py \
    --bam_dir /path/to/bam_files \
    --annotation_file peaks.bed \
    --output_dir results \
    --calcnorm \
    --calcnorm_method upperquartile \
    --threads 10
```

calcnorm CPMを使用するには、Rおよび必要なパッケージ（edgeR）がインストールされている必要があります。詳細については、[CPM方式詳細ガイド](./cpm_README_ja.md)を参照してください。

### 主要な引数

| 引数名 | 説明 | デフォルト値 |
|--------|------|-------------|
| **必須引数** |||
| `--bam_dir` / `-b` | BAMファイルが格納されているディレクトリ | - |
| `--annotation_file` / `-a` | アノテーションファイル (BED, SAF, merge_SE形式) | - |
| `--output_dir` / `-o` | 結果を保存するディレクトリ | - |
| **アノテーション関連** |||
| `--is_mergese_format` | アノテーションファイルがmerge_SE.tsv形式であることを指定 | False |
| **ファイル/サンプル選択** |||
| `--filename_pattern` | BAMファイル名のワイルドカードパターン | "*.bam" |
| `--sample_delimiter` | サンプル名抽出のための区切り文字列 | None |
| **ツールパス・スレッド数** |||
| `--samtools_path` | samtoolsの実行ファイルへのパス | "samtools" |
| `--featurecounts_path` | featureCountsの実行ファイルへのパス | "featureCounts" |
| `--threads` / `-T` | スレッド数 | 4 |
| `--single_end` | シングルエンドモードで実行（指定しない場合はペアエンド） | False |
| **edgeR正規化CPM（calcnorm）関連** |||
| `--calcnorm` | Standard log2-CPM計算の代わりにedgeR正規化CPMを使用する | False |
| `--calcnorm_method` | 正規化手法（upperquartile, TMM, RLE, none） | "upperquartile" |
| `--min_cpm` | フィルタリングのための最小CPM閾値 | 1.0 |
| `--min_samples` | CPM閾値を超えるべき最小サンプル数（0=フィルタリングなし） | 0 |
| `--calcnorm_output_name` | calcnorm CPM出力ファイル名（空文字列で保存しない） | "calcnorm.tsv" |

全てのオプションを確認するには、以下のコマンドを実行してください：
```bash
python cpm_peakprep.py --help
```

より詳細な設定オプションや出力の説明(jupyter notebook実行例含む)については、[CPM方式詳細ガイド](./cpm_README_ja.md)を参照してください。また、edgeR正規化CPMの詳細については[edgeR正規化CPMガイド](./cpm_calcnorm_README_ja.md)を参照してください。

## BigWig方式（bigwig_peakprep）

BigWig方式は、`deeptools`を使用してBAMファイルを処理し、ゲノム領域内のシグナル値を定量化します。**単一サンプルの解析**と**サンプル-コントロールの差分解析**の2種類のパイプラインを提供します。

### パイプラインの構造

BigWig方式は2つの主要な実行パスを提供しています：

1. **単一サンプル解析パイプライン** (`bigwig_peakprep.py`)
   - BAMファイルから正規化されたbigWigファイルを生成し、指定領域のシグナル値を定量化
   - 内部的に以下のスクリプトを順次実行：
     1. `bigwig_peakprep_bamconvert.py`: BAMファイルからbigWigファイルへの変換
     2. `bigwig_peakprep_summary.py`: bigWigファイルからピーク領域のカウント値テーブルを生成

2. **差分解析パイプライン** (`bigwig_peakprep_bamdiff.py`)
   - ChIP-seqサンプルとInputコントロールのlog2比較解析を実行
   - 内部的に以下のスクリプトを順次実行：
     1. `bigwig_peakprep_bamdiff_generatediff.py`: サンプルとコントロールのlog2比bigWigファイルを生成
     2. `bigwig_peakprep_summary.py`: 生成された差分bigWigファイルから領域カウント値テーブルを生成

**ユーティリティモジュール**:
- `bigwig_peakprep_utils.py`: 基本的な関数と共通ユーティリティを提供
- `bigwig_peakprep_bamdiff_utils.py`: 差分解析に特化した関数を提供

### どのパイプラインを選ぶべきか

- **単一サンプル解析パイプライン**（`bigwig_peakprep.py`）を使用する場合：
  - 各サンプルの独立した解析や比較を行いたい場合
  - 個々のBAMファイルから正規化されたシグナル値を取得したい場合

- **差分解析パイプライン**（`bigwig_peakprep_bamdiff.py`）を使用する場合：
  - ChIP-seqとInputコントロールを比較し、バックグラウンドノイズを除去した解析を行いたい場合
  - サンプルとコントロールのペアがある場合（例：T1-H3K27acとT1-Input）

デフォルトでは、bamCoverageは**RPGC (Reads Per Genomic Content)** 正規化を使用しますが、`--normalize_using`オプションで他の正規化方法（CPM, BPM, RPKM, RPGC, None）に変更可能です。

### 基本的な使用方法

#### 単一サンプル解析の場合（bigwig_peakprep.py）

```bash
python bigwig_peakprep.py \
    --bam_dir="/path/to/bam_files" \
    --annotation_file="/path/to/peaks.bed" \
    --output_dir="/path/to/output" \
    --tools_dir="/path/to/deeptools/bin" \
    --sample_delimiter="-H3K27ac.sorted.bam" \
    --threads=10
```

#### サンプル-コントロール差分解析の場合（bigwig_peakprep_bamdiff.py）

```bash
python bigwig_peakprep_bamdiff.py \
    --sample_bam_dir="/path/to/chipseq_bams" \
    --control_bam_dir="/path/to/input_bams" \
    --annotation_file="/path/to/peaks.bed" \
    --output_dir="/path/to/output" \
    --tools_dir="/path/to/deeptools/bin" \
    --sample_delimiter="-H3K27ac.sorted.bam" \
    --threads=10
```

### 主要な引数

#### 単一サンプル解析（bigwig_peakprep.py）

| 引数名 | 説明 | デフォルト値 |
|--------|------|-------------|
| **必須引数** |||
| `--bam_dir` | BAMファイルが格納されているディレクトリ | - |
| `--annotation_file` | アノテーションファイル（BED, SAF, merge_SE形式） | - |
| `--output_dir` | 結果を保存するディレクトリ | - |
| **ツールパス関連** |||
| `--tools_dir` | deeptoolsの実行ファイルを含むディレクトリ | None |
| `--bamcoverage_path` | bamCoverageの実行ファイルへのパス | "bamCoverage" |
| `--multibigwigsummary_path` | multiBigwigSummaryの実行ファイルへのパス | "multiBigwigSummary" |
| **実行ステップ制御** |||
| `--run_bamconvert_only` | BAM→bigWig変換のみを実行 | False |
| `--run_summary_only` | bigWigサマリー処理のみを実行 | False |
| **正規化関連** |||
| `--normalize_using` | bamCoverageの正規化方法 (RPGC, CPM, BPM, RPKM, None) | "RPGC" |
| `--effective_genome_size` | RPGC正規化用の有効ゲノムサイズ | 2913022398 |

#### サンプル-コントロール差分解析（bigwig_peakprep_bamdiff.py）

| 引数名 | 説明 | デフォルト値 |
|--------|------|-------------|
| **必須引数** |||
| `--sample_bam_dir` | ChIP-seqサンプルBAMファイルが格納されているディレクトリ | - |
| `--control_bam_dir` | Inputコントロール（Input/IgG）BAMファイルが格納されているディレクトリ | - |
| `--annotation_file` | アノテーションファイル（BED, SAF, merge_SE形式） | - |
| `--output_dir` | 結果を保存するディレクトリ | - |
| **ファイル選択** |||
| `--sample_pattern` | サンプルBAMファイルのグロブパターン | "*.bam" |
| `--control_pattern` | コントロールBAMファイルのグロブパターン | "*Input*.bam" |
| **処理関連** |||
| `--operation` | bamCompareの操作方法（log2, ratio, subtract, add, meanなど） | "log2" |
| `--bin_size` | bamCompareのbinサイズ（ベース単位） | 50 |
| `--pseudocount` | ゼロ除算を避けるための擬似カウント値 | 1.0 |
| `--run_diff_only` | 差分bigWig生成のみを実行（サマリー処理をスキップ） | False |

全てのオプションを確認するには、以下のコマンドを実行してください：
```bash
python bigwig_peakprep.py --help
python bigwig_peakprep_bamdiff.py --help
```

より詳細な設定オプションや出力の説明については、[BigWig方式の詳細ドキュメント](./bigwig_README_ja.md)を参照してください。

## 出力フォーマット

> **重要な注意:** 入力アノテーションファイルとしてSAF形式（1-based座標）を使用した場合でも、**全ての出力テーブルの座標（Start列）は0-based（BED形式基準）に統一されます。**

### CPM方式の出力

1. **Standard log2-CPM** (デフォルト: `logCPM.tsv`):
   - 各ピーク/領域におけるログ変換済みCPM値を含むテーブル
   - 列: PeakID、Chr、Start（0-based、BED形式基準）、End、各サンプルのlog2-CPM値

2. **edgeR正規化CPM（calcnorm）** (デフォルト: `calcnorm.tsv`):
   - edgeRの正規化メソッドを適用したCPM値を含むテーブル
   - 列: PeakID、Chr、Start（0-based、BED形式基準）、End、各サンプルの正規化CPM値

3. **中間ファイル**:
   - flagstat出力ファイル（総マップ済みリード数）
   - featureCounts出力ファイル（生のリードカウント）

### BigWig方式の出力

#### 単一サンプル解析の場合（bigwig_peakprep.py）

1. **bigWigファイル** (ステップ1の出力):
   - 各BAMファイルから生成された正規化されたbigWigファイル
   - ファイル名形式: `{サンプル名}.{正規化方法}.bigwig`
   - 例: `sample1.RPGC.bigwig`

2. **bam_to_bigwig_details.tsv**:
   - ステップ1の詳細出力ファイル（ステップ2への入力）
   - コメント行にはパラメータ情報と実行時メタデータが含まれる
   - 列: `Sample_name`, `BAM_filename`, `BigWig_filename`, `BAM_fullpath`, `BigWig_fullpath`

3. **カウントテーブル** (ステップ2の出力):
   - `multibigwig_counts.tsv`: 生のカウントテーブル
   - `multibigwig_counts_log2.tsv`: log2変換済みカウントテーブル

#### サンプル-コントロール差分解析の場合（bigwig_peakprep_bamdiff.py）

1. **差分bigWigファイル**:
   - サンプルとコントロールのlog2比を表すbigWigファイル
   - ファイル名形式: `{sample_name}_vs_{control_name}.log2ratio.bigwig`
   - 例: `T1-H3K27ac_vs_T1-Input.log2ratio.bigwig`

2. **bamdiff_to_bigwig_details.tsv**:
   - 差分解析の詳細情報を含むTSVファイル
   - カラム: `Sample_name`, `BigWig_filename`, `BigWig_fullpath`
   - 重要: このファイルには `already_log2_transformed: true` メタデータが含まれる（サマリー処理での追加log2変換を防ぐため）

3. **変換プラン**:
   - `conversion_plan.yaml`: 差分解析の実行計画
   - `conversion_plan.tsv`: 変換計画のTSV形式（表形式）

4. **カウントテーブル**:
   - `multibigwig_counts.tsv`: 各ピーク領域でのlog2比の値を含むテーブル
   - 注意: すでにlog2変換済みのため、追加のlog2変換は行われない

出力テーブルの形式例:
```
PeakID  Chr     Start   End     Sample1  Sample2  Sample3
peak1   chr1    1000    2000    2.45     3.12     1.87
peak2   chr1    3000    4000    4.21     3.98     4.56
...
```

### 出力ファイル構造の例

#### 単一サンプル解析の出力ディレクトリ構造:
```
output_dir/
├── bigwig/
│   ├── sample1.RPGC.bigwig
│   ├── sample2.RPGC.bigwig
│   └── ...
├── summary/
│   ├── bigwig_summary.npz
│   └── bigwig_summary.tab
├── bam_to_bigwig_details.tsv
├── multibigwig_counts.tsv
├── multibigwig_counts_log2.tsv
└── pipeline_run.log
```

#### 差分解析の出力ディレクトリ構造:
```
output_dir/
├── bigwig/
│   ├── T1-H3K27ac_vs_T1-Input.log2ratio.bigwig
│   ├── T2-H3K27ac_vs_T2-Input.log2ratio.bigwig
│   └── ...
├── summary/
│   ├── bigwig_summary.npz
│   └── bigwig_summary.tab
├── bamdiff_to_bigwig_details.tsv
├── bamdiff_to_bigwig_details.yaml
├── conversion_plan.yaml
├── conversion_plan.tsv
├── multibigwig_counts.tsv
└── pipeline_run.log
```

## トラブルシューティング

よくある問題と解決策：

1. **ツールが見つからないエラー**:
   - samtoolsやfeatureCounts、deeptoolsがインストールされているか確認
   - 完全なパスを対応するオプションで指定
   ```bash
   --tools_dir="/path/to/deeptools/bin"
   --samtools_path="/path/to/samtools"
   --featurecounts_path="/path/to/featureCounts"
   ```

2. **BAMファイルが見つからない**:
   - `--filename_pattern`引数、`--sample_pattern`、`--control_pattern`を確認
   - ファイルの権限と読み取りアクセスを確認

3. **サンプル名抽出の問題**:
   - `--sample_delimiter`の値を適切に設定し、イコール記号（`=`）を使用して指定
   - ハイフンを含む値の場合は、必ずイコール記号を使用する: `--sample_delimiter="-H3K27ac.sorted.bam"`

4. **メモリ不足エラー**:
   - より多くのメモリを持つマシンで実行
   - 処理するファイル数を減らす、またはバッチに分けて処理

5. **処理が遅い場合**:
   - `--threads`オプションでスレッド数を増やす
   - インデックス付きのBAMファイルを使用する
   - `--bin_size`値を大きくすると処理が速くなりますが、解像度は低下します

6. **サンプルとコントロールのマッチング失敗**:
   - ログで警告メッセージ "Could not create sample-control pairs" が表示される場合、マッチングに失敗しています
   - ファイル名パターンを `--sample_pattern` と `--control_pattern` で明示的に指定する
   - ファイル名に共通の識別子（T1, T2など）があることを確認する
   - 単一のコントロールを全サンプルに使用する場合は、コントロールディレクトリに1つのファイルだけが含まれるようにする

## 詳細ドキュメント
* [cpm_README_ja.md](./cpm_README_ja.md) – CPM方式の詳細な使い方
* [cpm_calcnorm_README_ja.md](./cpm_calcnorm_README_ja.md) – edgeR正規化CPMのアルゴリズム解説
* [bigwig_README_ja.md](./bigwig_README_ja.md) – BigWig方式の詳細な使い方

## 引用 / Citation

このツールを研究に使用する場合は、以下を引用してください：

Shinkai, N., Asada, K., Machino, H., Takasawa, K., Takahashi, S., Kouno, N., Komatsu, M., Hamamoto, R., & Kaneko, S. (2025). SEgene identifies links between super enhancers and gene expression across cell types. *npj Systems Biology and Applications*, 11(1), 49. https://doi.org/10.1038/s41540-025-00533-x

## ライセンス

このプログラムはMITライセンスで公開されています。詳細については、[LICENSE](https://github.com/hamamoto-lab/SEgene/blob/main/LICENSE)をご覧ください。