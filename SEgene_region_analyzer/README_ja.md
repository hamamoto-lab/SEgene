# SEgene_RegionAnalyzer

*(For the English version of this README, please see [README.md](https://github.com/hamamoto-lab/SEgene/blob/main/SEgene_region_analyzer/README.md).)*

**SEgene_RegionAnalyzer** は、**SEgene プロジェクト**のコンポーネントの一つであり、公共データベースの情報を利用して、ユーザーが関心を持つゲノム領域におけるスーパーエンハンサー（Super-Enhancer: SE）の活性や関連性を評価・分析するためのPythonツールパッケージです。

> **注意**: このツールは**現在開発中**であり、プログラムのファイル構成や使用方法はバージョンによって大きく変化する可能性があります。

SEgeneプロジェクトは主に以下のコンポーネントで構成されており、SEgene_RegionAnalyzerはその分析ツールチェーンを拡張するものです：
* [peak_to_gene_links](https://github.com/hamamoto-lab/SEgene/tree/main/peak_to_gene_links)：遺伝子発現とエンハンサーピーク間の相関情報を分析
* [SE_to_gene_links](https://github.com/hamamoto-lab/SEgene/tree/main/SE_to_gene_links)：相関情報を用いたスーパーエンハンサーの分析
* **SEgene_RegionAnalyzer**（本ツール）：SEgeneで同定されたSE領域の詳細な特性評価と公共データベースとの統合分析

本ツールは、SEgeneプロジェクトのワークフローを拡張し、以下のような解析を可能にします：

1. **入力データ**：
   * [SE_to_gene_links](https://github.com/hamamoto-lab/SEgene/tree/main/SE_to_gene_links) モジュールが生成したTSV形式のマージSE領域リスト
   * または一般的なBED形式で定義したゲノム領域

2. **解析プロセス**：
   * 指定された領域と公共データを統合
   * 組織タイプ別等のエンリッチメント解析を実行
   * 統計的有意性の評価

3. **出力結果**：
   * 詳細なTSVファイル形式の解析レポート
   * 視覚的に理解しやすいHTML形式のレポート
   * 領域ごとの統計情報と図表

現在は **ヒトゲノム hg38** および **SEdb 2.0** のデータに対応しています。今後、対応データベースや解析機能は拡張・アップデートされる予定です。

現在の対応データベース：
* [SEdb 2.0](http://www.licpathway.net/sedb/) - ヒトとマウスのスーパーエンハンサーの包括的データベース

## 目次
- [特徴](#特徴-features)
- [動作環境](#動作環境-requirements)
- [セットアップ](#セットアップ-setup)
- [使い方](#使い方-usage)
- [出力形式](#出力形式)
- [今後の開発予定](#今後の開発予定-future-development)
- [ライセンス](#ライセンス-license)
- [引用](#引用-citation)

## 特徴 (Features)

* **バッチ処理:** TSV または BED ファイルで定義された複数のゲノム領域を一括で解析
* **公共データベース連携:** スーパーエンハンサー（遺伝子発現制御に強力に影響を与える大規模なエンハンサー領域）データおよびサンプルメタデータを利用
* **組織エンリッチメント解析:** 指定領域に関連する SE が特定の組織タイプで濃縮されているかを統計的に評価
* **柔軟な入力:** SEgene 由来の TSV ファイル、または標準的な BED ファイルで解析対象領域を指定可能
* **豊富な出力:**
    * 各領域の解析結果（図、統計、関連データ）を個別のディレクトリに出力
    * 全領域の結果を集約した TSV ファイルを出力
    * 各領域の HTML レポートを生成（オプション）

## 動作環境 (Requirements)

* **Python:** 3.10 以降を推奨
* **主な Python ライブラリ:**
    * pandas: データ操作と解析
    * numpy: 数値計算 
    * scipy: 科学技術計算、統計解析
    * statsmodels: 統計モデル構築
    * pybedtools: BED ファイル操作 
    * matplotlib / seaborn: グラフ描画
    * Jinja2: HTML レポート生成用テンプレートエンジン
    * pyarrow : Parquet ファイルの読み書きに必要

    **注**: 現在は `requirements.txt` や `environment.yml` が提供されていません。将来的に提供される可能性はありますが、現時点では上記のライブラリを手動でインストールする必要があります。

* **外部ツール :**
    * [Bedtools](https://bedtools.readthedocs.io/): pybedtools が内部で使用します。

* **データ:**
    * [SEdb 2.0](http://www.licpathway.net/sedb/) からダウンロードした以下のファイルが必要です:
        * `human_sample_information_sedb2.txt` (サンプルメタデータ)
        * `SE_package_hg38.bed` (SE 定義 BED ファイル)

## セットアップ (Setup)

 1. **リポジトリのクローン:**
    ```bash
    git clone https://github.com/hamamoto-lab/SEgene.git
    cd SEgene/SEgene_region_analyzer # クローンしたリポジトリ内のSEgene_region_analyzerディレクトリに移動
    ```

2. **Python 環境の準備:**
    
    以下はcondaまたはpipを使用した環境構築の例です。必要なライブラリを手動でインストールする必要があります。

    **conda を使用する場合:**
    ```bash
    # 環境の作成
    conda create -n sedb_analyzer python=3.10
    # 環境のアクティベート
    conda activate sedb_analyzer
    # 主要ライブラリのインストール
    conda install pandas numpy scipy statsmodels matplotlib seaborn jinja2
    conda install -c bioconda pybedtools
    conda install -c conda-forge pyarrow fastparquet
    ```

    **pip を使用する場合:**
    ```bash
    # 仮想環境の作成
    python -m venv venv
    # 環境のアクティベート
    source venv/bin/activate  # Linux/macOS
    # または
    .\venv\Scripts\activate  # Windows
    # 主要ライブラリのインストール
    pip install pandas numpy scipy statsmodels matplotlib seaborn jinja2 pybedtools pyarrow fastparquet
    ```

    **注**: 将来的には `environment.yml` や `requirements.txt` が提供される可能性があります。

3. **SEdb データのダウンロードと準備:**
* [SEdb 2.0 Download ページ](http://www.licpathway.net/sedb/download.php)にアクセスします。
* "Sample information" セクションから `human_sample_information_sedb2.txt` をダウンロードします。
* "Package of SE, SE element and TE" セクションから Human (hg38) の "Super-enhancer information" (`SE_package_hg38.bed`) をダウンロードします。
* (サーバーが混雑している場合は時間をおいて試してください)
* **重要:** このツールは、`SEgene_region_analyzer` ディレクトリの **親ディレクトリ** に `data/SEdb` という構造でデータファイルが配置されていることを想定しています。ダウンロードしたファイルを以下のようなディレクトリ構造になるように配置してください (`your_project_root` は `SEgene_region_analyzer` が存在する場所です)。

```
your_project_root/
├── SEgene_region_analyzer/  <-- このリポジトリのコードがある場所
│   ├── run_batch_processor.py
│   └── ... (他のスクリプト)
└── data/                 <-- SEgene_region_analyzer と同じ階層に作成
    └── SEdb/
        ├── human_sample_information_sedb2.txt
        └── SE_package_hg38.bed 
```

`SEgene_region_analyzer` ディレクトリと同じ階層に `data/SEdb` ディレクトリを作成し、ダウンロードしたファイルを配置するコマンド例:
```bash
# SEgene_region_analyzer があるディレクトリに移動していると仮定
cd .. # 親ディレクトリへ移動
mkdir -p data/SEdb
mv /path/to/downloaded/human_sample_information_sedb2.txt ./data/SEdb/
mv /path/to/downloaded/SE_package_hg38.bed ./data/SEdb/
cd SEgene_region_analyzer # 元のディレクトリに戻る (任意)
```

## 使い方 (Usage)

解析は主にバッチスクリプト `run_batch_processor.py` を使用して行います。**このスクリプトは `SEgene_region_analyzer` ディレクトリ内から実行する必要があります。**

まず、`SEgene_region_analyzer` ディレクトリに移動します。
```bash
cd path/to/SEgene_region_analyzer
```

以下に基本的な実行例と主なオプションの説明を示します。

### 基本的な実行例 (TSV 入力)

*(データファイル (`--bed_file`, `--metadata_file`) は `SEgene_region_analyzer` の親ディレクトリにある `data/SEdb` を参照する例です)*

```bash
python run_batch_processor.py \
    --input_file <入力TSVファイルのパス> \
    --input_format tsv \
    --bed_file ../data/SEdb/SE_package_hg38_cleaned.bed \
    --metadata_file ../data/SEdb/human_sample_information_sedb2.txt \
    --output_dir <メイン出力ディレクトリパス> \
    # --max_rows <整数>         # (オプション)
```

### 基本的な実行例 (BED 入力)

*(データファイル (`--bed_file`, `--metadata_file`) は `SEgene_region_analyzer` の親ディレクトリにある `data/SEdb` を参照する例です)*

```bash
python run_batch_processor.py \
    --input_file <入力BEDファイルのパス> \
    --input_format bed \
    --bed_file ../data/SEdb/SE_package_hg38_cleaned.bed \
    --metadata_file ../data/SEdb/human_sample_information_sedb2.txt \
    --output_dir <メイン出力ディレクトリパス> \
    # --max_rows <整数>         # (オプション)
```
*(注: `<入力TSVファイルのパス>`, `<入力BEDファイルのパス>`, `<メイン出力ディレクトリパス>` は実際のパスに置き換えてください。データファイルのパス (`--bed_file`, `--metadata_file`) は、セットアップで配置した場所に合わせて調整してください。)*


### 主なオプション

* **`--input_file <パス>`**: (必須) 解析対象の領域リストを含むファイルのパス。
  * 領域の定義をTSVまたはBED形式で指定します。
  * TSV形式の場合は、SEgeneの出力形式を想定しています。

* **`--input_format <tsv|bed>`**: (オプション, デフォルト: tsv) `--input_file` の形式を指定。
  * `tsv`: タブ区切りファイル形式で、通常SEgeneの出力が対象
  * `bed`: 標準的なBED形式で、独自に定義したゲノム領域の指定が可能

* **`--bed_file <パス>`**: (必須) SEdb からダウンロードした SE BED ファイルのパス。
  * このファイルには、すべてのスーパーエンハンサー領域の定義が含まれています。
  * 例: `../data/SEdb/SE_package_hg38.bed` （`SEgene_region_analyzer` の親ディレクトリの `data/SEdb` にある場合）

* **`--metadata_file <パス>`**: (必須) SEdb からダウンロードしたサンプルメタデータファイルのパス。
  * このファイルには、各サンプルの詳細情報（組織タイプなど）が含まれています。
  * 例: `../data/SEdb/human_sample_information_sedb2.txt` （`SEgene_region_analyzer` の親ディレクトリの `data/SEdb` にある場合）

* **`--output_dir <パス>`**: (必須) すべての出力（中間結果、最終結果）を保存するメインディレクトリ。
  * このディレクトリ内に階層的に解析結果が保存されます。

* **`--max_rows <整数>`**: (オプション) テスト用に処理する入力ファイルの最大行数を指定。
  * 大規模なファイルを使用する前に、小さなサブセットでテストするのに便利です。

* **`--debug`**: (オプション) デバッグモードを有効化。
  * 通常の使用では不要ですが、問題がある場合は詳細なログ出力に役立ちます。

* **`--log_file <パス>`**: (オプション) このバッチスクリプト自体のログを保存するファイルパスを指定。
  * 指定しない場合はコンソールにのみログが出力されます。

### 入力ファイル形式

#### TSV 形式 (`--input_format tsv`)
* タブ区切りファイル。
* 最初の列に領域 ID が含まれている必要があります。現在のスクリプトは `chrX_START_END` 形式（例: `chr1_1109435_1174178`）を想定してパースします。
* 他の列は保持され、最終的な出力 TSV にも含まれます。

#### BED 形式 (`--input_format bed`)
* 標準的な BED 形式（タブ区切り）。
* 最初の3列（染色体、開始位置、終了位置）が領域情報として使用されます。ヘッダー行はないものとして処理されます。
* 現在の実装では BED ファイルの他の列（名前など）は最終結果には引き継がれません。

## 出力形式

`--output_dir` で指定したディレクトリ内に以下の構造で出力が生成されます。

```
<メイン出力ディレクトリ>/
├── 00_preprocessed_data/       # Preprocessor が生成した背景データパッケージ
│   ├── data/
│   │   ├── bed_data.parquet
│   │   ├── sample_info.parquet
│   │   └── cell_id_counts.parquet
│   ├── statistics/
│   │   ├── statistics.json
│   │   └── ... (他の統計ファイル)
│   └── metadata.json
├── 01_batch_region_results/    # 各領域の解析結果
│   ├── region_1_chr1_.../      # リージョン 1 の結果
│   │   ├── extractor_output/   # Extractor の出力
│   │   │   ├── bed_data.parquet
│   │   │   ├── sample_info.parquet
│   │   │   ├── distributions/
│   │   │   │   ├── tissue_counts.parquet
│   │   │   │   └── ...
│   │   │   ├── summary.json
│   │   │   └── metadata.json
│   │   └── analyzer_output/    # Analyzer の出力
│   │       ├── machine_readable/
│   │       │   ├── data/
│   │       │   │   ├── enrichment_results.parquet
│   │       │   │   └── ...
│   │       │   └── metadata.json
│   │       ├── human_readable/
│   │       │   ├── figures/
│   │       │   │   ├── tissue_enrichment.png/svg
│   │       │   │   └── ...
│   │       │   ├── tables/
│   │       │   │   ├── enrichment_results.csv/txt
│   │       │   │   └── ...
│   │       │   └── summary.txt
│   │       ├── summary.json
│   │       └── index.html (HTMLレポート)
│   └── region_2_chr7_.../      # リージョン 2 の結果
│       └── ...
└── batch_analysis_results_YYYYMMDD_HHMMSS.tsv  # ★最終的な集約結果 TSV
```

### 最終TSVファイルの主な追加カラム

`batch_analysis_results_YYYYMMDD_HHMMSS.tsv` には、以下のような情報が追加されます：

* **`se_unique_sample_count`**: Extractor が領域内で見つけた SE に関連するユニークなサンプル ID の数。
* **`linked_sample_count`**: 上記 ID のうち、メタデータと正常に紐付けられたサンプルの数。
* **`significant_tissues`**: Analyzer によって有意にエンリッチしていると判定された組織タイプのリスト（カンマ区切り）。
* **`significant_tissues_fc`**: 対応する Fold Change のリスト（カンマ区切り）。
* **`significant_tissues_pvalue`**: 対応する P 値のリスト（カンマ区切り）。
* **`significant_tissues_fdr`**: 対応する補正済み P 値 (FDR) のリスト（カンマ区切り）。
* **`processing_error`**: 処理中にエラーが発生した場合のエラーメッセージ（一部）。

## 今後の開発予定 (Future Development)

SEgene_RegionAnalyzerは現在も開発が進められており、以下のような拡張が予定されています：

* **対応データベースの拡張**: 現在はSEdb 2.0のみ対応していますが、他のスーパーエンハンサーデータベースにも順次対応予定です。
* **解析機能の強化**: より多様な統計解析手法の実装や視覚化機能の強化を計画しています。
* **他のゲノムアセンブリへの対応**: mm10(マウス)への対応も検討中です。

## ライセンス (License)

このプログラムはMITライセンスで公開されています。詳細については、[LICENSE](https://github.com/hamamoto-lab/SEgene/blob/main/LICENSE)をご覧ください。

## 引用 / Citation

このツールを研究に使用する場合は、以下を引用してください：

Shinkai, N., Asada, K., Machino, H., Takasawa, K., Takahashi, S., Kouno, N., Komatsu, M., Hamamoto, R., & Kaneko, S. (2025). SEgene identifies links between super enhancers and gene expression across cell types. *npj Systems Biology and Applications*, 11(1), 49. https://doi.org/10.1038/s41540-025-00533-x