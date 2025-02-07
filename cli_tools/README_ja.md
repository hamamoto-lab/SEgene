# CLI Tools for SE-Gene Correlation Analysis
# SE-Gene相関解析用CLIツール群

*(For the English version of this README, please see [README.md](https://github.com/hamamoto-lab/SEgene/blob/main/cli_tools/README.md).)*

このディレクトリには、ChIP-Seqデータやデータを用いて以下のような一連の解析を行うためのCLIスクリプトが含まれています。

1. featureCount用GTFファイルの生成
2. BAMファイルリストの生成
3. featureCountsジョブの実行
4. BAMファイル全体のリードカウント取得
5. featureCounts出力のマージとCPM計算
6. CPMデータの0-indexed変換
7. SE–遺伝子の相関解析

以下の手順は、あらかじめ[SEgene_package](https://github.com/hamamoto-lab/SEgene/blob/main/SE_to_gene_links/README_ja.md)（[English](https://github.com/hamamoto-lab/SEgene/tree/main/SE_to_gene_links)）で作成された`temp_full_df_info.pkl`や、RNAのTPMデータ、ChIP-SeqのBAMファイル群を前提としています。

`temp_full_df_info.pkl`は、[SEgene_package](https://github.com/hamamoto-lab/SEgene/blob/main/SE_to_gene_links/README_ja.md)のSEgeneJupyterインスタンスで以下のメソッドを実行することにより取得できます。

```python
SEgeneJupyterインスタンス.analyze_merge_SE_info(save_info_pkl="temp_full_df_info.pkl")
```

[SEgene_package](https://github.com/hamamoto-lab/SEgene/blob/main/SE_to_gene_links/README_ja.md)のより詳細な使用方法については、[チュートリアルノートブック](https://github.com/hamamoto-lab/SEgene/tree/main/SE_to_gene_links/notebooks)を参照してください。

## ファイル構成

```
cli_tools/
├── bam_count.sh                                    # BAMファイルのリードカウント取得
├── convert_merge_featurecount_CPM_zero_indexed.py  # CPMデータの0-indexed変換
├── correlation_analysis.py                         # SE-遺伝子間の相関解析
├── generate_file_list.sh                          # BAMファイルリスト生成
├── generate_gtf.py                                # GTFファイル生成
├── gtf_processing.py                              # GTF処理用ユーティリティ
├── merge_featurecount_CPM.py                      # featureCounts結果のマージとCPM計算
├── run_featurecounts_array.sh                     # featureCountsジョブ実行（アレイジョブ用）
└── run_featurecounts_local.sh                     # featureCountsジョブ実行（ローカル実行用）
```

## 動作環境

- **[featureCounts](http://subread.sourceforge.net/)** (v2.0.6)と **[samtools](http://www.htslib.org/)** (v1.15.1)がインストール済みで、パスが通っていること。
- **Python 3.10.8**および関連ライブラリ（`pandas`、`numpy`、`matplotlib`など）がインストールされていること。
- 各スクリプトは`python <script>.py --help`でヘルプを確認できます。

## 前提ファイル

- **`temp_full_df_info.pkl`**：[SEgene_package](https://github.com/hamamoto-lab/SEgene/blob/main/SE_to_gene_links/README_ja.md)のJupyter実行により生成されたpklファイル。
- **RNAのTPMファイル（CSV）**：相関解析用の遺伝子発現量。
- **ChIP-SeqのBAMファイル**：featureCounts用に必要。

BAMファイルやRNA CSVのフォーマットはあらかじめ適切に整形されている必要があります。RNA CSVは[P2GL](https://github.com/hamamoto-lab/SEgene/blob/main/peak_to_gene_links/README_ja.md)や[SEgene_package](https://github.com/hamamoto-lab/SEgene/blob/main/SE_to_gene_links/README_ja.md)で使用するものと同一です。

## 解析の流れ

以下、SEgeneリポジトリをクローンし、相関解析用の作業ディレクトリ（calculate_correlation）を作成して解析を行います。（なお、`/path/to/calculate_correlation`は任意の作業ディレクトリパスに置き換えてください）

```bash
# リポジトリのクローン
git clone https://github.com/hamamoto-lab/SEgene.git

# 相関解析用の作業ディレクトリを作成し、必要なスクリプトをコピー
mkdir -p /path/to/calculate_correlation
cp SEgene/cli_tools/* /path/to/calculate_correlation/
cd /path/to/calculate_correlation
```

以下の手順で解析を進めます。

### 1. featureCount用GTFファイルの生成

`temp_full_df_info.pkl`を使ってfeatureCounts用GTFを作成します。

```bash
python generate_gtf.py
```

出力例：
- `output/top_regions_df.pkl`：SE数（もしくはサンプル数等）がランキング上位のマージSE領域の情報を保存（デフォルトは上位20）。
- `output/unsorted_temp.gtf`：一時ファイル。
- `output/modified_for_featurecounts.gtf`：top_regions_df.pklに含まれるマージSE領域に対応するfeatureCounts用GTFファイル。

### 2. BAMファイルリストの生成

BAMファイルがあるディレクトリを指定し、ファイル一覧を作成します。

```bash
./generate_file_list.sh /path/to/bam_directory
```

- `output/bam_list.txt`にBAMファイルのフルパス一覧が出力されます。
- スクリプトに実行権限が必要になる場合があるため、`chmod +x generate_file_list.sh`を行ってください。

### 3. featureCountsジョブの実行

`output/bam_list.txt`を用いて配列ジョブ等で一括処理を行う一例を示します。

featureCountsの実行には以下の2つの方法があります。

#### 3.1 アレイジョブとして実行（計算機クラスタ向け）

```bash
qsub -t 1-$(wc -l < output/bam_list.txt) run_featurecounts_array.sh
```

提供されている`run_featurecounts_array.sh`は以下のようなOracle Grid Engine（OGE、旧Sun Grid Engine: SGE）向けのスクリプトです。

```bash
#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -jc your_job_class   # ジョブクラスの指定
#$ -o ./output/featurecounts_output/logs/
#$ -e ./output/featurecounts_output/logs/

# BAMファイルの取得（SGE_TASK_IDを使用）
bam_file=$(sed -n "${SGE_TASK_ID}p" output/bam_list.txt)
...
```

注意事項：
- この例はOracle Grid Engine向けのコマンドとスクリプトです。
- 使用するジョブスケジューラに応じて実行コマンドを変更してください。
  - SLURMの場合：`sbatch --array=1-$(wc -l < output/bam_list.txt) run_featurecounts_array.sh`
  - PBS/Torqueの場合：`qsub -J 1-$(wc -l < output/bam_list.txt) run_featurecounts_array.sh`
- スクリプトファイル自体も環境に応じて調整が必要です。
  - SLURMの場合：`#SBATCH`ディレクティブを使用し、`${SGE_TASK_ID}`を`${SLURM_ARRAY_TASK_ID}`に変更。
  - PBS/Torqueの場合：`#PBS`ディレクティブを使用し、`${SGE_TASK_ID}`を`${PBS_ARRAY_INDEX}`に変更。
- メモリ量、コア数、実行時間制限などのリソース指定は、使用する計算機環境に合わせて適切に調整してください。

#### 3.2 ローカルで逐次実行

単一マシンでの実行の場合は、`run_featurecounts_local.sh`を使用できます。

1. 実行権限を付与（必要な場合）
```bash
chmod +x run_featurecounts_local.sh
```

2. スクリプトを実行
```bash
./run_featurecounts_local.sh
```

- BAMファイルを1つずつ順番に処理します。
- `output/featurecounts_output/`配下に結果が出力されます。
- ログは`output/featurecounts_output/logs/`に保存されます。
- 必要に応じてスクリプト内のパラメータやリソース設定を編集してください。

### 4. BAMファイルのリードカウント取得

BAMファイルのリストを元に、各BAMファイルのリード数を取得し、CSVにまとめます。

```bash
./bam_count.sh output/bam_list.txt
```

実行結果：
- `output/bam_read_counts.csv`：各BAMのリード数。

### 5. featureCounts出力のマージとCPM計算

`merge_featurecount_CPM.py`を用いて、featureCountsの結果をマージし、リード数CSV（`bam_read_counts.csv`）を使ってCPMを計算します。

```bash
python merge_featurecount_CPM.py
```

出力例：
- `output/merge_featurecount_CPM.pkl`：マージ＆CPM計算後のpklファイル。

### 6. 0-indexedデータの作成

`convert_merge_featurecount_CPM_zero_indexed.py`を実行してデータを0-indexedに変換します。

```bash
python convert_merge_featurecount_CPM_zero_indexed.py --export_csv
```

- `--export_csv`を付けると同名ベースのCSVも出力されます。
- デフォルトで`output/merge_featurecount_CPM_zero_indexed.pkl`が生成されます。

### 7. 相関解析の実施

`correlation_analysis.py`を使い、SE（CPM）と遺伝子（TPM）のピアソン相関を計算してグラフを生成します。

```bash
python correlation_analysis.py \
  --top_regions_df ./output/top_regions_df.pkl \
  --se_count 10 \  # ランキング上位から解析するマージSE数を指定\
  --output_dir cli_cor_test \
  --merge_featurecount_pkl ./output/merge_featurecount_CPM_zero_indexed.pkl \
  --rna_csv path/to/gene_tpm.csv \
  --full_info_pkl temp_full_df_info.pkl \
  --use_all_samples
```

- `--se_count`：ランキング上位から解析するマージSE数を指定します。
- `--use_all_samples`を付けると、`SAMPLE_distinct`に依存せず共通サンプルを全て使用します。
- 結果は指定した`--output_dir`に保存されます。
- `--save_png`オプションを付ければPNG形式でも保存可能です。

## 補足

- 各スクリプトの引数や動作は`--help`オプションやソースコードのコメントを参照してください。
- BAMファイルやRNA CSVのフォーマットはあらかじめ適切に整形されている必要があります。RNA CSVは[P2GL](https://github.com/hamamoto-lab/SEgene/blob/main/peak_to_gene_links/README_ja.md)や[SEgene_package](https://github.com/hamamoto-lab/SEgene/blob/main/SE_to_gene_links/README_ja.md)で使用するものと同一です。
- ジョブ管理システム（SGE、SLURMなど）の設定に合わせて、スクリプト内のジョブ送信方法を調整してください。
- 各処理の途中で必要に応じて権限変更（`chmod +x`）やパス設定を行ってください。

以上の手順を踏むことで、ChIP-SeqおよびRNAデータを用いた一連の解析（GTF生成、featureCountsジョブ、CPM計算、相関解析）を実行することが可能です。