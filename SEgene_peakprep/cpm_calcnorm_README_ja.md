# cpm_peakprep - edgeR正規化CPMガイド

このドキュメントでは、cpm_peakprepの**edgeR正規化CPM**（calcnorm CPM）機能の概要と基本的な使用方法について説明します。基本的な概要やCPMメソッドの一般情報については、[CPM方式詳細ガイド](./cpm_README_ja.md)を参照してください。

## 目次

1. [概要](#概要)
2. [要件と依存関係](#要件と依存関係)
3. [処理フロー](#処理フロー)
4. [コマンドライン引数](#コマンドライン引数)
5. [出力ファイルの形式](#出力ファイルの形式)
6. [使用例](#使用例)
7. [関連ドキュメント](#関連ドキュメント)

## 概要

edgeR 正規化 CPM（calcnorm CPM）機能は、R 言語の **edgeR** パッケージを用いて ChIP‑seq ピークのカウントデータを正規化します。`--calcnorm` オプションを指定するとこの機能が有効になります。  
> **動作条件:** R ≥ 4.2、edgeR ≥ 3.40、rpy2 ≥ 3.5 がインストールされている必要があります。

### Standard log2-CPMとcalcnorm CPMの違い

1. **Standard log2-CPM（従来の方法）**:
   - samtoolsのflagstatで取得した総マップリード数に基づいて計算
   - 計算式: CPM = (count × 10⁶) / total_mapped_reads
   - その後、log₂(CPM + pseudocount)で対数変換

2. **edgeR正規化CPM（calcnorm CPM）**:
   - featureCountsで計数された該当領域のリードに基づく計算
   - edgeRの`calcNormFactors()`関数で算出した正規化係数を使用
   - 各サンプルのカウント分布特性に基づくスケーリングファクターを適用

## 要件と依存関係

この機能を使用するには、以下のソフトウェアが必要です：

1. **R言語**: バージョン4.2.0以上
2. **rpy2**: バージョン3.5.0以上（Python-R接続ライブラリ）
3. **edgeR**: バージョン3.40.0以上（Bioconductorパッケージ）

R/edgeRのインストールには複数の方法があります。例えば、conda-forgeを使用する場合：

```bash
# conda-forgeからのインストール例
conda install -c conda-forge r-base=4.2
conda install -c conda-forge r-essentials
conda install -c bioconda bioconductor-edger=3.40
conda install -c conda-forge rpy2

# 代替方法：pipでのインストール
# pip install rpy2>=3.5.0
```

各環境に合わせたインストール方法については、R/Bioconductorの公式ドキュメントを参照してください。

## 処理フロー

edgeR正規化CPMの処理は、以下のステップで行われます：

1. **データ前処理**:
   - featureCountsの出力ファイルからカウントデータを読み込み
   - サンプル名のクリーニング（BAMファイル名から意味のあるサンプル名を抽出）

2. **DGEList作成**:
   - edgeRの`DGEList`オブジェクトにカウントデータを変換

3. **フィルタリング**（オプション）:
   - CPM値に基づいて特定の基準を満たすピークを選択可能
   - パラメータ: `min_cpm`（最小CPM値）と`min_samples`（最小サンプル数）

4. **正規化係数の計算**:
   - 選択した方法で正規化係数を計算
   - ユーザーが選択可能なメソッド: upperquartile（デフォルト）、TMM、RLE、none

5. **CPM計算**（log変換しない）:
   - 正規化されたライブラリサイズに基づいてCPM値を算出

6. **結果出力**:
   - 標準フォーマット（PeakID, Chr, Start, End, 各サンプルの値）で保存
   - オプションで全メタデータ列を含む完全版も保存可能

### 正規化手法について

edgeRパッケージは、以下の主要な正規化手法を提供しています：

- **Upper Quartile**（デフォルト）: 各サンプルのカウント分布の75パーセンタイルに基づく正規化
- **TMM**: Trimmed Mean of M-values。サンプル間のログ比のトリム平均に基づく正規化
- **RLE**: Relative Log Expression。ピークごとの幾何平均を基準とした正規化
- **none**: 正規化係数を適用せず

各正規化手法の詳細な理論的背景については、edgeRのドキュメントを参照してください。

### フィルタリングについて

`min_cpm`と`min_samples`パラメータを使用すると、特定の基準を満たすピークのみを選択できます：

- 各ピークについて「CPM > min_cpm」となるサンプル数が「min_samples」以上の場合にそのピークが保持されます
- デフォルト設定は「min_cpm=1.0, min_samples=0」で、実質的にフィルタリングは行われません

## コマンドライン引数

edgeR正規化CPM機能で使用可能な主なコマンドライン引数は以下の通りです：

### 基本的な引数

| 引数名 | 説明 | デフォルト値 |
|--------|------|-------------|
| `--calcnorm` | edgeR正規化CPMを有効化するフラグ | `False` |
| `--calcnorm_method` | 正規化手法：upperquartile、TMM、RLE、none | `upperquartile` |
| `--calcnorm_output_name` | calcnorm CPM出力ファイル名。空文字列で保存しない | `calcnorm.tsv` |

### フィルタリング関連の引数

| 引数名 | 説明 | デフォルト値 |
|--------|------|-------------|
| `--min_cpm` | フィルタリングのための最小CPM閾値 | `1.0` |
| `--min_samples` | CPM閾値を超えるべき最小サンプル数。0=フィルタリングなし | `0` |

### サンプル名クリーニング関連の引数

| 引数名 | 説明 | デフォルト値 |
|--------|------|-------------|
| `--remove_extensions` | 一般的なBAMファイル拡張子（`.bam`, `.sorted.bam`など）を自動的に除去 | `False` |
| `--pattern_rules` | サンプル名クリーニングのための追加パターンルールを含むJSONファイル | `None` |

### 出力形式関連の引数

| 引数名 | 説明 | デフォルト値 |
|--------|------|-------------|
| `--full_metadata_output` | 標準出力に加えて全featureCountsメタデータ列を含む完全な出力を保存 | `None` |

## 出力ファイルの形式

### 標準出力ファイル（calcnorm.tsv）

標準的な出力形式は以下の構造を持ちます：

```
PeakID  Chr     Start   End     Sample1  Sample2  Sample3
peak1   chr1    999     2000    2.45     3.12     1.87
peak2   chr1    2999    4000    4.21     3.98     4.56
...
```

- **注意**: `Start`列は0-based（BED形式基準）で出力されます

### 完全メタデータ出力ファイル

`--full_metadata_output`を指定すると、featureCountsの全メタデータカラムを含む完全な出力が生成されます：

```
Geneid  Chr     Start   End     Strand  Length  Sample1  Sample2  Sample3
peak1   chr1    999     2000    +       1001    2.45     3.12     1.87
peak2   chr1    2999    4000    +       1001    4.21     3.98     4.56
...
```

この出力形式は、Strand情報やピーク長など、追加のメタデータが必要な場合に利用できます。

### サンプル名クリーニングのためのJSONファイル

複雑なサンプル名パターンがある場合、`--pattern_rules`オプションでJSONファイルを指定できます：

```json
[
  {"pattern": "\\.mLb\\.clN\\.sorted\\.bam$", "replace": ""},
  {"pattern": "^Sample_", "replace": ""},
  {"pattern": "-rep([0-9])", "replace": "_R$1"}
]
```

各ルールは、`pattern`（正規表現パターン）と`replace`（置換文字列）のペアで構成されます。

## 使用例

### 両方のCPM値を保存する基本例

```bash
python cpm_peakprep.py \
    --bam_dir /path/to/bam_files \
    --annotation_file peaks.bed \
    --output_dir results \
    --calcnorm \
    --threads 10
```

この例では、Standard log2-CPM（`logCPM.tsv`）とcalcnorm CPM（`calcnorm.tsv`）の両方が保存されます。

### TMM正規化を使用する例

```bash
python cpm_peakprep.py \
    --bam_dir /path/to/bam_files \
    --annotation_file peaks.bed \
    --output_dir results \
    --calcnorm \
    --calcnorm_method TMM \
    --threads 10
```

### フィルタリングを適用する例

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

### 完全メタデータ出力を保存する例

```bash
python cpm_peakprep.py \
    --bam_dir /path/to/bam_files \
    --annotation_file peaks.bed \
    --output_dir results \
    --calcnorm \
    --full_metadata_output \
    --threads 10
```

## 関連ドキュメント

- [メインREADME](./README_ja.md): ツール全体の概要
- [CPM方式詳細ガイド](./cpm_README_ja.md): CPM方式の基本的な使用方法と詳細情報
- [BigWig方式の詳細ガイド](./bigwig_README_ja.md): BigWig方式の詳細情報

## ライセンス

このプログラムはMITライセンスで公開されています。詳細については、[LICENSE](https://github.com/hamamoto-lab/SEgene/blob/main/LICENSE)をご覧ください。