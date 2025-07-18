# インストールガイド

## システム要件

### OS
- Linux (Ubuntu 18.04以降推奨)
- macOS (10.14以降)
- Windows (WSL2推奨)

### ソフトウェア要件
- Python 3.11以上
- pip 20.0以上
- Git
- bedtools（pybedtoolsの動作に必要）

### ハードウェア要件
- メモリ: 8GB以上（大規模解析では16GB以上推奨）
- ストレージ: 10GB以上の空き容量

## インストール手順

### 1. bedtoolsのインストール

pybedtoolsを使用するため、システムにbedtoolsがインストールされている必要があります。

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

#### Conda環境
```bash
conda install -c bioconda bedtools
```

### 2. SEgene_analyzerのインストール

#### 方法1: Conda/Miniforgeを使用（推奨）

バイオインフォマティクスワークフローではconda/miniforgeの使用を強く推奨します：

```bash
# miniforgeのインストール（未インストールの場合）
# 参照: https://github.com/conda-forge/miniforge

# conda環境の作成
conda create -n segene-analyzer python=3.11
conda activate segene-analyzer

# bedtoolsをconda経由でインストール（推奨）
conda install -c bioconda bedtools pybedtools

# リポジトリをクローン
git clone https://github.com/hamamoto-lab/SEgene.git
cd SEgene/SEgene_analyzer

# Python依存関係をインストール
pip install pandas matplotlib seaborn scipy statsmodels networkx japanize-matplotlib

# 開発モードでインストール
pip install -e .
```

#### 方法2: pipを使用（代替手段）

```bash
# リポジトリをクローン
git clone https://github.com/hamamoto-lab/SEgene.git
cd SEgene/SEgene_analyzer

# 依存関係をインストール
pip install -r requirements.txt

# 開発モードでインストール
pip install -e .
```

### 3. SEdbデータのダウンロード

SEgene_analyzerにはSEdb 2.0データベースファイルが必要です：

1. **アクセス**: [SEdb 2.0 Download Page](http://www.licpathway.net/sedb/download.php)
2. **必要なファイルをダウンロード**:
   - `human_sample_information_sedb2.txt`（サンプル情報）
   - `SE_package_hg38.bed`（ヒトhg38用スーパーエンハンサーパッケージ）

3. **データディレクトリを作成**:
   ```bash
   mkdir -p data/SEdb
   mv human_sample_information_sedb2.txt data/SEdb/
   mv SE_package_hg38.bed data/SEdb/
   ```

### 4. インストールの確認

```bash
# sedb-analyzerが使用可能か確認
sedb-analyzer --help

# バージョン情報をテスト
sedb-analyzer --version

# 依存関係を確認
python -c "import pandas, numpy, scipy, pybedtools; print('All dependencies installed successfully')"
```

## インストール確認

### 基本機能テスト

簡単なテストでインストールを確認します:

```bash
# テストディレクトリを作成
mkdir -p test_output cache

# 統計情報の準備をテスト（SEdbデータが必要）
sedb-analyzer prepare-stats \
  -b data/SEdb/SE_package_hg38.bed \
  -s data/SEdb/human_sample_information_sedb2.txt \
  -o cache/test_stats.pkl

# 成功すればエラーなしで出力が表示されます
```

### サンプルデータでのテスト

```bash
# サンプル領域ファイルを作成
cat > test_regions.tsv << EOF
test_region1	chr1	1000000	2000000
test_region2	chr2	3000000	4000000
EOF

# 小規模なバッチ解析を実行
sedb-analyzer batch \
  -b data/SEdb/SE_package_hg38.bed \
  -s data/SEdb/human_sample_information_sedb2.txt \
  -r test_regions.tsv \
  -o test_output/ \
  --max-rows 2
```

## トラブルシューティング

### 1. pybedtoolsのインストールエラー

```
Error: Unable to find bedtools executable
```

**解決方法**: bedtoolsが正しくインストールされ、PATHに含まれているか確認

```bash
# bedtoolsの確認
which bedtools
bedtools --version

# PATHに追加（必要な場合）
export PATH=$PATH:/path/to/bedtools/bin
```

### 2. Pythonバージョンの問題
**エラー**: `Python 3.11+ required`

**解決方法**:
```bash
# Pythonバージョンを確認
python --version

# 古いバージョンを使用している場合、Python 3.11以上をインストール
# condaを使用:
conda install python=3.11
```

### 3. 権限エラー
**エラー**: `Permission denied when installing`

**解決方法**:
```bash
# conda環境を使用
conda create -n segene-analyzer python=3.11
conda activate segene-analyzer
pip install -r requirements.txt
pip install -e .
```

### 4. メモリ不足エラー
**エラー**: `Out of memory during analysis`

**解決方法**:
```bash
# 小さなバッチサイズを使用
sedb-analyzer batch --max-rows 50 ...

# 統計情報を事前計算
sedb-analyzer prepare-stats ...
sedb-analyzer batch --use-cached-stats ...
```

### 5. SEdbデータダウンロードの問題
**エラー**: `Cannot download SEdb data`

**解決方法**:
- インターネット接続を確認
- 異なる時間帯にダウンロードを試行（サーバーが混雑している場合がある）
- 利用可能な場合は代替ダウンロード方法を使用

### パフォーマンス最適化

#### 1. メモリ使用量
```bash
# メモリ使用量を監視
htop
# または
top

# メモリ使用量を削減
sedb-analyzer batch --max-rows 100 --no-peak-analysis ...
```

#### 2. 処理速度
```bash
# 統計キャッシュを使用
sedb-analyzer prepare-stats ...
sedb-analyzer batch --use-cached-stats ...

# 出力を最小化
sedb-analyzer batch --figure-formats png --no-tables ...
```

## 高度なインストール

### 開発者向けインストール

開発者がコードを修正したい場合：

```bash
# 開発用依存関係を含めてクローン
git clone https://github.com/hamamoto-lab/SEgene.git
cd SEgene/SEgene_analyzer

# 開発モードでインストール
pip install -e .

# テスト用依存関係をインストール
pip install pytest pytest-cov

# テストを実行
python tests/run_tests.py
```


## SEgene_analyzerの更新

### Gitからの更新

```bash
cd SEgene/SEgene_analyzer
git pull origin main
pip install -r requirements.txt
pip install -e .
```

### 依存関係の更新

```bash
# 全依存関係を更新
pip install --upgrade -r requirements.txt

# 特定のパッケージを更新
pip install --upgrade pandas numpy scipy
```

## アンインストール

### パッケージの削除

```bash
# SEgene_analyzerをアンインストール
pip uninstall segene-analyzer

# conda環境を削除（使用していた場合）
conda env remove -n segene-analyzer

# データファイルを削除（オプション）
rm -rf data/SEdb
```

### クリーンインストール

```bash
# キャッシュと一時ファイルを削除
rm -rf cache/ results/ output/

# __pycache__ディレクトリを削除
find . -name "__pycache__" -type d -exec rm -rf {} +
```


## 関連ドキュメント

- [使用方法ガイド](usage_ja.md) - 完全な使用説明書
- [使用方法ガイド](usage.md)（英語）- 包括的なコマンドリファレンス
- [使用例](../examples/README_ja.md) - 実世界のワークフロー
- [メインドキュメント](../README_ja.md) - プロジェクト概要