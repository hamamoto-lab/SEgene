# 使用例

このディレクトリには SEgene_analyzer_erna の使用方法を示すサンプルスクリプトが含まれています。

## 利用可能な例

### 1. コマンドラインインターフェース (`cli_usage.sh`)

`erna-analyzer` CLI ツールの包括的な使用例：
- 単一領域解析
- キャッシュあり/なしのバッチ処理
- パフォーマンス向上のための統計事前計算
- メタデータレポート生成
- 様々な出力形式と閾値

```bash
# 例を確認
cat cli_usage.sh

# 特定の例を実行
erna-analyzer single -b data/peak -m metadata/eRNAbase_data.parquet --chr chr7 --start 1000000 --end 2000000
```

### 2. 基本的な使用方法 (`basic_usage.py`)

単一のゲノム領域を解析するための基本的なワークフローを示します：
- ERNARegionAnalyzer の初期化
- eRNAbase メタデータの読み込み
- BED ファイルの読み込み
- 特定領域の解析
- レポート生成

```bash
python basic_usage.py
```

### 3. バッチ解析 (`batch_analysis.py`)

TSV 入力ファイルを使用して複数のゲノム領域を解析する方法を示します：
- 領域 TSV ファイルの作成または使用
- 複数領域のバッチ解析実行
- すべての領域に対する包括的なレポート生成

```bash
python batch_analysis.py
```

### 4. 本番用ワークフロー (`production_workflow.py`)

エラーハンドリングとログ出力を備えた本番対応スクリプト：
- コマンドライン引数の解析
- 包括的なエラーハンドリング
- 進捗追跡
- 詳細なログ出力

```bash
python production_workflow.py --bed-dir data/peak --metadata metadata/eRNAbase_data.parquet --regions regions.tsv
```

### 5. サンプルデータ生成 (`create_sample_data.py`)

開発とテスト用のテストデータを生成します：
- 合成 BED ファイルの作成
- CSV と Parquet 両形式でのメタデータ生成
- テスト領域リストの作成

```bash
python create_sample_data.py --num-samples 10 --num-regions 5
```

## 入力データ要件

### eRNAbase メタデータ
- 形式: CSV または Parquet
- 必須カラム: サンプル情報、組織タイプ、細胞タイプなど

### BED ファイル
- 標準的な BED 形式（ゲノム座標）
- 異なるサンプル用の複数 BED ファイルを含むディレクトリ

### 領域 TSV（バッチ解析用）
- タブ区切りファイル（カラム: `region_name`, `chrom`, `start`, `end`）
- 例:
  ```
  region_name	chrom	start	end
  region_1	chr1	1000000	2000000
  region_2	chr2	5000000	6000000
  ```

## 例の実行方法

1. **データの準備**: eRNAbase メタデータと BED ファイルがあることを確認
2. **パスの更新**: 例中のファイルパスを実際のデータを指すように修正
3. **依存関係のインストール**: 必要なパッケージがすべてインストールされているか確認
4. **スクリプトの実行**: このディレクトリから例を実行

## 出力

例を実行すると、以下を含む出力ディレクトリが作成されます：
- 重複 eRNA データ（CSV 形式）
- 統計解析結果
- 可視化プロット（PNG 形式）
- 包括的な解析を含む HTML レポート