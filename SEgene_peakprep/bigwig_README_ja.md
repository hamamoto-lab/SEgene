# SEgene_peakprep - BigWig方式の詳細ガイド

このドキュメントは、SEgene_peakprepの**BigWig方式**に特化した詳細な使用方法と実行例を提供します。基本的な概要や共通要件については、[メインREADME](./README_ja.md)を参照してください。

## BigWig方式の概要

BigWig方式は、`deeptools`の`bamCoverage`を使用してBAMファイルを正規化されたbigWigファイルに変換し、続いて`multiBigwigSummary`を使用して指定領域のカウント値を抽出します。この方法は次のような特長があります：

- ゲノム全体にわたるカバレッジデータの可視化が可能
- さまざまな正規化方法（RPGC, CPM, BPM, RPKM）に対応
- マルチサンプルデータの効率的な処理と比較

## パイプラインの構造

BigWig方式は以下の3つのスクリプトから構成されています：

1. **bigwig_peakprep.py**: メインのラッパースクリプト（全体の処理を制御）
2. **bigwig_peakprep_bamconvert.py**: BAMファイルからbigWigファイルへの変換を担当
3. **bigwig_peakprep_summary.py**: bigWigファイルからマルチサンプルカウントテーブルを生成

この3段階の処理フローにより、BAMファイルから正規化されたbigWigファイルを生成し、指定されたゲノム領域のシグナル値を集計して定量化します。各ステップを個別に実行することも、ラッパースクリプトを通して一連の処理として実行することも可能です。

## 全引数リスト

コマンドライン実行時に使用できる全ての引数とその詳細です：

| 引数名 | 説明 | デフォルト値 |
|--------|------|-------------|
| **必須引数** |||
| `--bam_dir` | BAMファイルが格納されているディレクトリのパス（BAM→bigWig変換時に必須） | - |
| `--annotation_file` | アノテーションファイル（BED, SAF, または merge_SE形式）のパス | - |
| `--output_dir` | 結果を保存するディレクトリ | - |
| `--bigwig_details_file` | サマリー処理のみを実行する場合（`--run_summary_only`）に必須のbam_to_bigwig_details.tsvファイルのパス | None |
| **ツールパス関連** |||
| `--tools_dir` | deeptoolsの実行ファイルを含むディレクトリ | None |
| `--bamcoverage_path` | bamCoverageの実行ファイルへのパス（`--tools_dir`指定時は無視） | "bamCoverage" |
| `--multibigwigsummary_path` | multiBigwigSummaryの実行ファイルへのパス（`--tools_dir`指定時は無視） | "multiBigwigSummary" |
| **ファイル/サンプル選択** |||
| `--filename_pattern` | BAMファイル名のワイルドカードパターン | "*.bam" |
| `--sample_delimiter` | BAMファイル名からサンプル名を抽出するための区切り文字列 | None |
| **bamCoverage関連** |||
| `--bigwig_dir` | bigWigファイルを保存するディレクトリ | output_dir/bigwig |
| `--bin_size` | bamCoverageのbinサイズ（ベース単位） | 50 |
| `--normalize_using` | bamCoverageの正規化メソッド（RPGC, CPMなど） | "RPGC" |
| `--effective_genome_size` | RPGC正規化用の有効ゲノムサイズ | 2913022398 |
| **multiBigwigSummary関連** |||
| `--summary_dir` | multiBigwigSummaryファイルを保存するディレクトリ | output_dir/summary |
| `--summary_basename` | multiBigwigSummary出力ファイルのベース名 | "bigwig_summary" |
| **データテーブル関連** |||
| `--force_chr_start_end_ids` | BEDファイルにname列があっても、chr_start_end形式をPeakIDとして使用 | False |
| `--pseudocount` | log2変換前に追加する擬似カウント | 1.0 |
| `--raw_counts_name` | 生カウントテーブルのファイル名 | "multibigwig_counts.tsv" |
| `--log2_counts_name` | log2変換カウントテーブルのファイル名 | "multibigwig_counts_log2.tsv" |
| **実行制御** |||
| `--threads` / `-T` | deeptoolsコマンドで使用するスレッド数 | 4 |
| `--force_overwrite` | 既存のbigWigファイルを上書き（指定しない場合はスキップ） | False |
| `--run_bamconvert_only` | BAM→bigWig変換のみを実行 | False |
| `--run_summary_only` | bigWigサマリー処理のみを実行 | False |
| **ログ設定** |||
| `--log_level` | ロギングレベル（DEBUG/INFO/WARNING/ERROR/CRITICAL） | "INFO" |
| `--script_log_file` | スクリプト実行ログのファイル名（空文字で無効化） | "pipeline_run.log" |

## 詳細な使用例

### ラッパースクリプトでの一括処理

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

### 2段階処理：BAM変換のみ実行

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

### 2段階処理：サマリー処理のみ実行

```bash
python bigwig_peakprep.py \
    --annotation_file="/path/to/peaks.bed" \
    --output_dir="/path/to/output" \
    --tools_dir="/path/to/deeptools/bin" \
    --bigwig_details_file="/path/to/output/bam_to_bigwig_details.tsv" \
    --threads=10 \
    --run_summary_only
```

### 個別ステップスクリプトを使った直接実行

個別スクリプトは、より細かな制御が必要な場合や、特定の処理ステップのみを実行したい場合に便利です。各スクリプトには独自のオプションがあります。詳細なオプションは各スクリプトに`--help`フラグを付けて実行することで確認できます。

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

## 出力ファイル詳細

### ステップ1（BAM to bigWig変換）の出力

1. **bigWigファイル**（`bigwig/`ディレクトリ内）:
   - BAMファイルごとに生成される正規化されたbigWigファイル
   - ファイル名形式: `{サンプル名}.{正規化方法}.bigwig`
   - 例: `sample1.RPGC.bigwig`

2. **bam_to_bigwig_details.tsv**:
   - ステップ1の詳細出力ファイル（ステップ2への入力）
   - コメント行にはパラメータ情報と実行時メタデータが含まれる
   - 列: `Sample_name`, `BAM_filename`, `BigWig_filename`, `BAM_fullpath`, `BigWig_fullpath`

### ステップ2（bigWigサマリーとカウントテーブル生成）の出力

1. **multiBigwigSummary出力**（`summary/`ディレクトリ内）:
   - `{summary_basename}.npz`: NumPy形式のバイナリデータ
   - `{summary_basename}.tab`: タブ区切りのテキストデータ

2. **カウントテーブル**:
   - `multibigwig_counts.tsv`: 生のカウントテーブル（PeakID情報付き）
   - `multibigwig_counts_log2.tsv`: log2変換済みカウントテーブル（PeakID情報付き）

カウントテーブルの形式例:
```
PeakID  Chr     Start   End     Sample1  Sample2  Sample3
peak1   chr1    999     2000    2.45     3.12     1.87
peak2   chr1    2999    4000    4.21     3.98     4.56
...
```

> **注意**: 出力ファイルの`Start`列は0-based（BED形式基準）座標になっています。

## 正規化方法の選択

`--normalize_using`オプションで以下の正規化方法を選択できます：

- **RPGC** (デフォルト): Reads Per Genomic Content。指定されたゲノムサイズを基準に正規化
- **CPM**: Counts Per Million mapped reads
- **BPM**: Bins Per Million mapped reads
- **RPKM**: Reads Per Kilobase per Million mapped reads
- **None**: 正規化なし

例えば、CPM正規化を使用する場合:
```bash
python bigwig_peakprep.py \
    --bam_dir="/path/to/bam_files" \
    --annotation_file="/path/to/peaks.bed" \
    --output_dir="/path/to/output" \
    --normalize_using="CPM" \
    --threads=10
```

## トラブルシューティング

よくある問題と解決策：

1. **ハイフンを含む引数値で発生するエラー**:
   - 引数のイコール渡し形式を使用してください: `--sample_delimiter="-H3K27ac.sorted.bam"`
   - スペース区切り形式は使わないでください: ~~`--sample_delimiter "-H3K27ac.sorted.bam"`~~

2. **bamCoverage/multiBigwigSummaryが見つからないエラー**:
   - `--tools_dir`オプションで正確なdeeptoolsのbinディレクトリを指定してください
   - 完全なパスを`--bamcoverage_path`と`--multibigwigsummary_path`で指定することもできます

3. **サンプル名抽出の問題**:
   - `--sample_delimiter`の値を適切に設定し、イコール記号（`=`）を使用して指定してください
   - BAMファイル名とサンプル名のマッピングを確認するには、ログレベルを`DEBUG`に設定してください

4. **メモリ不足エラー**:
   - 大量のBAMファイルや大きなゲノム領域を処理する場合は、より多くのメモリを持つマシンで実行してください
   - 処理するファイル数を減らすか、バッチに分けて処理することを検討してください

5. **処理が遅い場合**:
   - `--threads`オプションでスレッド数を増やしてください
   - インデックス付きのBAMファイルを使用してください
   - `--bin_size`値を大きくすると処理が速くなりますが、解像度は低下します

## 引用

このツールを研究に使用する場合は、以下を引用してください：
**(論文は現在準備中です)**

## ライセンス

このプログラムはMITライセンスで公開されています。詳細については、[LICENSE](./LICENSE)をご覧ください。