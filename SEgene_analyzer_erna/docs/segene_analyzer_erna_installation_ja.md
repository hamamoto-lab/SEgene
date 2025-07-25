# インストールガイド

## システム要件

### オペレーティングシステム
- Linux (Ubuntu 18.04+ 推奨)
- macOS (10.14+)
- Windows (WSL2 推奨)

### ソフトウェア要件
- Python 3.11+
- pip 20.0+
- Git
- bedtools

## インストール手順

### 1. bedtools のインストール

#### Ubuntu/Debian
```bash
sudo apt-get update
sudo apt-get install bedtools
```

#### macOS (Homebrew)
```bash
brew install bedtools
```

#### CentOS/RHEL
```bash
sudo yum install bedtools
```

#### Conda 環境
```bash
conda install -c bioconda bedtools
```

### 2. SEgene_analyzer_erna のインストール

#### 方法1: Conda/Miniforge を使用（推奨）

```bash
# miniforge のインストール（未インストールの場合）
# 参照: https://github.com/conda-forge/miniforge

# bedtools を含む conda 環境の作成
conda create -n segene-analyzer-erna python=3.11
conda activate segene-analyzer-erna

# bedtools のインストール
conda install -c bioconda bedtools pybedtools

# リポジトリのクローン
git clone https://github.com/your-org/SEgene_analyzer_erna.git
cd SEgene_analyzer_erna

# Python 依存関係のインストール
pip install pandas matplotlib seaborn scipy statsmodels networkx japanize-matplotlib natsort Jinja2

# 開発モードでインストール
pip install -e .
```

#### 方法2: pip を使用（代替案）

```bash
# リポジトリのクローン
git clone https://github.com/your-org/SEgene_analyzer_erna.git
cd SEgene_analyzer_erna

# 依存関係のインストール
pip install -r requirements.txt

# 開発モードでインストール
pip install -e .
```

### 3. eRNAbase データの準備

SEgene_analyzer_erna には eRNAbase データベースファイルが必要です：

1. **eRNAbase BED ファイル**: 複数の BED ファイルを含むディレクトリ（サンプルごとに1つ）
2. **メタデータファイル**: サンプルメタデータを含む Parquet または CSV ファイル
3. **領域リストファイル**: 解析対象のゲノム領域を含む TSV ファイル

#### データ構造の例
```
data/
├── peak/                           # BED ファイルディレクトリ
│   ├── Sample-01-0001.bed
│   ├── Sample-01-0002.bed
│   └── ...
├── eRNAbase_data.parquet          # メタデータファイル
└── regions.tsv                    # 解析対象領域
```

### 4. インストールの確認

```bash
# インポートテスト（開発モード）
python -c "from erna_analyzer import ERNARegionAnalyzer; print('✅ インストール成功')"

# またはソースインポートテスト
python -c "from src.erna_analyzer import ERNARegionAnalyzer; print('✅ インストール成功')"
```

## 使用例

### 基本的な使用方法
```python
from src.erna_analyzer import ERNARegionAnalyzer
# または pip install 後: from erna_analyzer import ERNARegionAnalyzer

# アナライザーの初期化
analyzer = ERNARegionAnalyzer(results_dir='results/erna_analysis')

# eRNAbase メタデータの読み込み
analyzer.load_erna_metadata('data/eRNAbase_data.parquet', species='human')

# TSV ファイルから領域を解析
result_file = analyzer.batch_analyze_regions_from_tsv(
    bed_directory='data/peak',
    metadata_file='data/eRNAbase_data.parquet',
    regions_tsv_file='data/regions.tsv',
    data_output_dir='results/data',
    report_output_dir='results/reports',
    max_rows=10,
    species='human'
)
```

## トラブルシューティング

### よくある問題

1. **インポートエラー: bedtools が見つからない**
   - 解決策: システムパッケージマネージャーまたは conda を使用して bedtools をインストール

2. **大規模データセットでのメモリ問題**
   - 解決策: システムメモリを増やすか、`max_rows` パラメータを使用してバッチサイズを削減

3. **権限エラー**
   - 解決策: 出力ディレクトリに対する書き込み権限があることを確認

