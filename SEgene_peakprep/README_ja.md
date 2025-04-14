# SEgene_peakprep

*(For the English version of this README, please see [README.md](./README.md).)*

**SEgene_peakprep** は、シーケンスデータ（BAMファイル）から指定されたゲノム領域におけるシグナル値を定量化・正規化し、下流のSEgene解析で使用する形式に整形するためのPythonパイプラインです。このツールは、特定の参照ゲノムには依存せず、例えばhg38などにマッピングされたデータで利用することができます。

## 開発状況

> **注意:** このツールは**現在開発中**であり、プログラムのファイル構成や使用方法はバージョンによって変化する可能性があります。

## 概要

SEgene_peakprepは、[SEgeneプロジェクト](https://github.com/hamamoto-lab/SEgene)において、初期データ準備の役割を果たします：

1. **入力**: マッピング済みシーケンスデータ（BAMファイル）とアノテーションファイル（BEDファイルなど）
2. **出力**: 正規化されたカウント値テーブル
3. **フロー**: この出力は、`peak_to_gene_links`モジュールの主要な入力となり、RNA-seqのTPMデータと組み合わせて遺伝子制御ネットワーク解析に使用されます

このパイプラインは以下の2つの実装方法を提供します：

1. **CPM方式**: `featureCounts`を使用して直接Counts Per Million（CPM）値を計算します
2. **BigWig方式**: `deeptools`を使用してBAMからbigWigファイルに変換し、multiBigwigSummaryでピークを定量化します

どちらの方法も、BED、SAF、merge_SE形式のアノテーションファイルに対応しています。用途に応じて適切な方法を選択できます。

## 共通要件

- **Python**: Python 3.10以上
- **基本ライブラリ**:
  - pandas
  - numpy
  - natsort
  
  ```bash
  # condaの場合
  conda install pandas numpy natsort
  
  # pipの場合
  pip install pandas numpy natsort
  ```

### 共通の入力ファイル

両方の実装方法で、以下の入力ファイルが必要です：

1. **BAMファイル**: マッピング済みのシーケンスデータファイル（推奨：インデックス付き.baiファイル）
2. **アノテーションファイル**: ゲノム領域を定義するファイルで、以下の形式に対応:
   - **BED形式**: 標準的なBEDフォーマット（0-based開始座標）
     ```
     # 例: peaks.bed（タブ区切り）
     chr1    999     2000    peak1    0    +
     chr1    2999    4000    peak2    0    +
     chr2    4999    6000    peak3    0    -
     ```
   
   - **SAF形式** (Simplified Annotation Format): タブ区切りのテキストファイル（**1-based開始座標**）
     ```
     # 例: peaks.saf（タブ区切り）
     GeneID  Chr     Start   End     Strand
     peak1   chr1    1000    2000    +
     peak2   chr1    3000    4000    +
     peak3   chr2    5000    6000    -
     ```
     **注意:** SAF形式の座標は1-based（Start=1000は0-basedでは999）ですが、出力は全て0-basedに統一されます
   
   - **merge_SE.tsv形式**: 1列目が `chr_start_end` 形式の領域ID（例：`chr1_1000_2000`）
     ```
     # 例: merged_se.tsv（タブ区切り）
     se_data            feature
     chr1_999_2000      feature1
     chr1_2999_4000     feature2
     chr2_4999_6000     feature3
     ```
     **注意:** merge_SE形式では座標は0-based（chr1_999_2000のStartは0-basedで999）と解釈されます

## CPM方式（cpm_peakprep）

CPM方式は、`featureCounts`を使用して指定ゲノム領域（BED/SAF/merge_SE.tsv形式）におけるリードを直接カウントし、CPM値とその対数変換値を計算します。

### 追加要件

- **外部ツール**:
  - samtools (推奨: v1.13以上、v1.15.1で動作確認済み)
  - featureCounts (Subreadパッケージの一部、推奨: v2.0.0以上、v2.0.6で動作確認済み)
  
  ```bash
  # condaでのインストール例
  conda install -c bioconda samtools subread
  ```

### 詳細ガイド

CPM方式の詳細な設定オプションや使用方法については、以下を参照してください：

### 基本的な使用方法

```bash
python cpm_peakprep.py \
    --bam_dir /path/to/bam_files \
    --annotation_file peaks.bed \
    --output_dir results \
    --filename_pattern "Sample*.sorted.bam" \
    --sample_delimiter=".sorted.bam" \
    --threads 10
```

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

全てのオプションを確認するには、以下のコマンドを実行してください：
```bash
python cpm_peakprep.py --help
```

より詳細な設定オプションや出力の説明(jupyter notebook実行例含む)については、[CPM方式の詳細ドキュメント](./cpm_README_ja.md)を参照してください。

## BigWig方式（bigwig_peakprep）

BigWig方式は、`deeptools`の`bamCoverage`を使用してBAMファイルを正規化されたbigWigファイルに変換し、続いて`multiBigwigSummary`を使用して指定領域のカウント値を抽出します。

### 追加要件

- **追加ライブラリ**:
  - pyranges
  
  ```bash
  # condaの場合
  conda install -c bioconda pyranges
  
  # pipの場合
  pip install pyranges
  ```

- **外部ツール**:
  - deeptools (推奨: v3.5.0以上、v3.5.5で動作確認済み)
  
  ```bash
  # condaでのインストール例
  conda install -c bioconda deeptools
  ```

### パイプラインの構造

BigWig方式は以下の3つのスクリプトから構成されています：

1. **bigwig_peakprep.py**: メインのラッパースクリプト（全体の処理を制御）
2. **bigwig_peakprep_bamconvert.py**: BAMファイルからbigWigファイルへの変換を担当
3. **bigwig_peakprep_summary.py**: bigWigファイルからマルチサンプルカウントテーブルを生成

この3段階の処理フローにより、BAMファイルから正規化されたbigWigファイルを生成し、指定されたゲノム領域のシグナル値を集計して定量化します。各ステップを個別に実行することも、ラッパースクリプトを通して一連の処理として実行することも可能です。

デフォルトでは、bamCoverageは**RPGC (Reads Per Genomic Content)** 正規化を使用しますが、`--normalize_using`オプションで他の正規化方法（CPM, BPM, RPKM, RPGC, None）に変更可能です。

### 基本的な使用方法

```bash
python bigwig_peakprep.py \
    --bam_dir="/path/to/bam_files" \
    --annotation_file="/path/to/peaks.bed" \
    --output_dir="/path/to/output" \
    --tools_dir="/path/to/deeptools/bin" \
    --sample_delimiter="-H3K27ac.sorted.bam" \
    --threads=10
```

### 主要な引数

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

全てのオプションを確認するには、以下のコマンドを実行してください：
```bash
python bigwig_peakprep.py --help
```

より詳細な設定オプションや出力の説明については、[BigWig方式の詳細ドキュメント](./bigwig_README_ja.md)を参照してください。

## 出力フォーマット

> **重要な注意:** 入力アノテーションファイルとしてSAF形式（1-based座標）を使用した場合でも、**全ての出力テーブルの座標（Start列）は0-based（BED形式基準）に統一されます。**

### CPM方式の出力

1. **logCPMテーブル** (デフォルト: `logCPM.tsv`):
   - 各ピーク/領域におけるログ変換済みCPM値を含むテーブル
   - 列: PeakID、Chr、Start（0-based、BED形式基準）、End、各サンプルのlogCPM値

2. **中間ファイル**:
   - flagstat出力ファイル（総マップ済みリード数）
   - featureCounts出力ファイル（生のリードカウント）

### BigWig方式の出力

1. **bigWigファイル** (ステップ1の出力):
   - 各BAMファイルから生成された正規化されたbigWigファイル
   - ファイル名形式: `{サンプル名}.{正規化方法}.bigwig`

2. **カウントテーブル** (ステップ2の出力):
   - `multibigwig_counts.tsv`: 生のカウントテーブル
   - `multibigwig_counts_log2.tsv`: log2変換済みカウントテーブル

出力テーブルの例:
```
PeakID  Chr     Start   End     Sample1  Sample2  Sample3
peak1   chr1    999     2000    2.45     3.12     1.87
peak2   chr1    2999    4000    4.21     3.98     4.56
...
```


## トラブルシューティング

よくある問題と解決策：

1. **ツールが見つからないエラー**:
   - samtoolsやfeatureCounts、deeptoolsがインストールされているか確認
   - 完全なパスを対応するオプションで指定

2. **BAMファイルが見つからない**:
   - `--filename_pattern`引数を確認

3. **サンプル名抽出の問題**:
   - `--sample_delimiter`の値を適切に設定し、イコール記号（`=`）を使用して指定
   - ハイフンを含む値の場合は、必ずイコール記号を使用する: `--sample_delimiter="-H3K27ac.sorted.bam"`

4. **メモリ不足エラー**:
   - より多くのメモリを持つマシンで実行
   - 処理するファイル数を減らす、またはバッチに分けて処理

5. **処理が遅い場合**:
   - `--threads`オプションでスレッド数を増やす
   - インデックス付きのBAMファイルを使用する

## 引用

このツールを研究に使用する場合は、以下を引用してください：
**(論文は現在準備中です)**

## ライセンス

このプログラムはMITライセンスで公開されています。詳細については、[LICENSE](./LICENSE)をご覧ください。