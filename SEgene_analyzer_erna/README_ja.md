# SEgene_analyzer_erna (開発版)

*(For the English version of this README, please see [README.md](README.md).)*

**SEgene_analyzer_erna** は、eRNAbase データを使用したenhancer RNA領域解析のための **eRNAbase専用** 開発版コンポーネントです。強化された機能と現代的な Python パッケージングを備えています。

> **⚠️ 開発状況**: これは **eRNAbase解析専用の開発版** です。SEdbベースの解析については [SEgene_analyzer](https://github.com/hamamoto-lab/SEgene/tree/main/SEgene_analyzer) の使用をご検討ください。

[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## プロジェクトコンテキスト

SEgene_analyzer_ernaは、eRNAbaseデータ解析に特化して設計された**SEgeneプロジェクト**エコシステムの一部です。

### 関連コンポーネント

- **親プロジェクト**: [SEgene](https://github.com/hamamoto-lab/SEgene) - 完全なスーパーエンハンサー解析プラットフォーム
- **SEdb版**: [SEgene_analyzer](https://github.com/hamamoto-lab/SEgene/tree/main/SEgene_analyzer) - SEdbベースの領域アナライザー
- **データソース**: [eRNAbase](https://bio.liclab.net/eRNAbase/index.php) - エンハンサーRNAデータベース

## 概要

SEgene_analyzer_ernaは、ゲノム領域と**eRNAbaseデータベース**の重複解析を行うバイオインフォマティクスツールです。指定したゲノム領域に対して、eRNAbaseに登録されているenhancer RNA（eRNA）の重複を検出し、組織特異性や細胞タイプ別の統計解析を実行します。

## eRNAbaseについて

本解析ツールは **eRNAbase** (https://bio.liclab.net/eRNAbase/index.php) より取得したデータを用いて実行します。eRNAbaseは、複数の種や組織タイプにわたる実験的に検証されたenhancer RNA（eRNA）のデータセットを提供する包括的なデータベースです。

### 必要なデータ

**必須のeRNA BEDファイル:**
- [eRNAbaseダウンロードページ](https://bio.liclab.net/eRNAbase/download.php)から `download_peak.zip` をダウンロード
- BEDファイルをディレクトリに展開（例: `data/peak/`）

**必須のメタデータファイル:**
- `metadata/` 内のメタデータファイルを使用:
  - `eRNAbase_data.csv` - 人間が読みやすい形式
  - `eRNAbase_data.parquet` - 圧縮形式（推奨）

**データバージョン:** このツールは **2025年7月22日時点** のeRNAbaseデータに対して開発・テストされています。

## 主な機能

- **eRNAbase重複解析**: 指定ゲノム領域とeRNAbaseデータの重複検出
- **組織特異性解析**: 組織・細胞タイプ別のエンリッチメント解析
- **統計的検定**: Fisher's exact testによる有意性検定
- **可視化機能**: 
  - リード分布プロット
  - エンリッチメント火山プロット
  - 組織/細胞タイプ分布チャート
- **HTMLレポート生成**: 包括的な解析結果レポート
- **バッチ処理**: 複数領域の一括解析

## ディレクトリ構造

```
SEgene_analyzer_erna/
├── src/                    # ソースコード
│   ├── __init__.py
│   ├── erna_analyzer.py    # メインのeRNA解析クラス
│   └── sedb_analyzer.py    # 基底クラス
├── docs/                   # ドキュメント
├── examples/               # 使用例
│   ├── basic_usage.py     # 基本的な使用例
│   ├── batch_analysis.py  # バッチ処理例
│   └── README.md          # 使用例の説明
├── tests/                  # テストファイル
├── requirements.txt        # Python依存関係
├── setup.py               # インストール設定
└── README.md              # このファイル
```

## インストール

### 必要条件

- Python 3.11以上
- bedtools（pybedtoolsが使用）

### 推奨: Conda/Miniforge を使用

バイオインフォマティクスワークフローには、conda/miniforge の使用を強く推奨します：

```bash
# conda環境の作成
conda create -n segene-analyzer-erna python=3.11
conda activate segene-analyzer-erna

# bedtoolsのインストール
conda install -c bioconda bedtools pybedtools

# リポジトリのクローンとインストール
git clone https://github.com/hamamoto-lab/SEgene_analyzer_erna.git
cd SEgene_analyzer_erna
pip install -e .
```

### 1. リポジトリのクローン

```bash
git clone [repository_url]
cd SEgene_analyzer_erna
```

### 2. 依存パッケージのインストール

```bash
# 開発モードでインストール（推奨）
pip install -e .

# または requirements.txt を使用
pip install -r requirements.txt
```

### 3. bedtoolsのインストール

Ubuntu/Debian:
```bash
sudo apt-get install bedtools
```

macOS:
```bash
brew install bedtools
```

## SEgene_analyzerとのデータ要件比較

| 項目 | SEgene_analyzer_erna | SEgene_analyzer |
|--------|----------------------|---------------------|
| **データソース** | eRNAbase | SEdb |
| **入力形式** | 複数BEDファイル | 単一BEDファイル |
| **メタデータ** | Parquet/CSV | TSV |
| **サンプル規模** | 数百サンプル | 数千サンプル |
| **用途** | eRNA特化解析 | 汎用SE解析 |

## 使用方法

### コマンドラインインターフェース（CLI）

インストール後、`erna-analyzer` コマンドが利用可能になります：

```bash
# ヘルプの表示
erna-analyzer --help

# 単一領域の解析
erna-analyzer single -b data/peak -m metadata/eRNAbase_data.parquet \
  --chr chr7 --start 1000000 --end 2000000 \
  --species human -o results/

# バッチ解析
erna-analyzer batch -b data/peak -m metadata/eRNAbase_data.parquet \
  -r regions.tsv -o results/ --max-rows 10

# 統計情報の事前計算（高速化のため）
erna-analyzer prepare-stats -b data/peak -m metadata/eRNAbase_data.parquet \
  -o cache/stats.pkl

# メタデータレポートの生成
erna-analyzer report -m metadata/eRNAbase_data.parquet -o report/
```

### Python API

```python
from src.erna_analyzer import ERNARegionAnalyzer
# またはインストール後: from erna_analyzer import ERNARegionAnalyzer

# アナライザーの初期化
analyzer = ERNARegionAnalyzer(results_dir='output_directory')

# eRNAbaseメタデータの読み込み
analyzer.load_erna_metadata('path/to/metadata.csv')  # CSVまたはParquet形式

# BEDファイルの読み込み
analyzer.get_bed_files('path/to/bed_files_directory/')

# 領域解析の実行
analyzer.analyze_region(
    chrom='chr1',
    start=1000000,
    end=2000000,
    region_name='my_region'
)

# HTMLレポートの生成
analyzer.generate_region_analysis_report()
```

### バッチ処理

複数の領域を一括で解析する場合：

```python
# TSVファイルから複数領域を解析
# TSVファイル形式: region_name\tchrom\tstart\tend
analyzer.batch_analyze_from_tsv('regions.tsv')
```

### 出力ファイル

解析結果は指定したディレクトリに以下の形式で保存されます：

```
output_directory/
├── region_name/
│   ├── overlapping_ernas.csv
│   ├── tissue_enrichment.csv
│   ├── plots/
│   │   ├── read_distribution.png
│   │   ├── enrichment_volcano.png
│   │   └── tissue_distribution.png
│   └── report.html
└── batch_summary.csv  # バッチ処理の場合
```

## 解析内容

### 1. eRNA重複検出
- 指定領域と重複するeRNAを特定
- 重複の長さと割合を計算

### 2. 組織特異性解析
- 各組織/細胞タイプでのeRNA出現頻度を集計
- Fisher's exact testによるエンリッチメント解析
- FDR補正による多重検定の調整

### 3. 可視化
- リード分布のスタックプロット
- エンリッチメントスコアの火山プロット
- 組織別eRNA数の棒グラフ

### 使用例

詳細な使用例は `examples/` ディレクトリを参照してください：

- `examples/basic_usage.py` - 基本的な単一領域解析
- `examples/batch_analysis.py` - 複数領域のバッチ処理
- `examples/README.md` - 詳細な説明

## 開発

### 開発環境のセットアップ

```bash
# 開発モードでインストール
pip install -e .

# 開発に必要な追加パッケージ
pip install pytest black flake8
```

### コード品質

```bash
# コードフォーマット
black src/

# リンター
flake8 src/
```

### テスト

```bash
# 基本的なテストの実行
python tests/test_erna_analyzer.py

# サンプルデータでのテスト
python examples/create_sample_data.py --num-samples 10
python examples/production_workflow.py --bed-dir sample_data/peak --metadata sample_data/eRNAbase_data.parquet
```

## 関連プロジェクト

- [SEgene](https://github.com/hamamoto-lab/SEgene) - 親プロジェクト
- [SEgene_analyzer](https://github.com/hamamoto-lab/SEgene/tree/main/SEgene_analyzer) - SEdb解析ツール（本ツールの基盤）

## ライセンス

このプログラムは MIT ライセンスで公開されています。詳細については、LICENSE ファイルを参照してください。

## 引用

このツールを研究に使用する場合は、以下を引用してください：

Shinkai, N., Asada, K., Machino, H., Takasawa, K., Takahashi, S., Kouno, N., Komatsu, M., Hamamoto, R., & Kaneko, S. (2025). SEgene identifies links between super enhancers and gene expression across cell types. *npj Systems Biology and Applications*, 11(1), 49. https://doi.org/10.1038/s41540-025-00533-x

詳細な引用情報と追加の参考文献については、[CITATION](https://github.com/hamamoto-lab/SEgene/blob/main/CITATION) ファイルを参照してください。

