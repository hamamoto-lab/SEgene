# Changelog (日本語)

本プロジェクトの変更点は、このファイルに記載します。

形式は [Keep a Changelog](https://keepachangelog.com/ja/1.1.0/) をベースとし、
[Semantic Versioning](https://semver.org/spec/v2.0.0.html) に準拠しています。

## [1.2.0] - 2025-03-11

### Added
- **SEランク分布解析のための新たなコア機能**:
  - 複数サンプル間で遺伝子に連結されたスーパーエンハンサーを検索する機能を追加
  - データセット内でのスーパーエンハンサーのランキング解析機能を実装
  - パーセンタイル計算および分布解析のための統計ツールを開発
  - サンプル間のSEランク分布を可視化するメソッドを作成
- **新機能を実演するチュートリアルノートブック**:
  - 英語版 (`tutorial_book_SE_rank_disribution.ipynb`)
  - 日本語版 (`tutorial_book_SE_rank_disribution_ja.ipynb`)
  - これらのノートブックは、特定の遺伝子に対応するスーパーエンハンサーのランキング状況（順位とパーセンタイル）を解析する方法をガイド

### Notes
- READMEファイルを更新し、新しい解析機能とノートブックに関する情報を追加

## [1.1.0] - 2025-02-07

### Added
- **新たな CLI ツール (`cli_tools/`)**: featureCounts 処理・CPM 計算・SE–遺伝子相関解析を効率化
  - `generate_gtf.py`, `generate_file_list.sh`, `run_featurecounts_array.sh` などのスクリプトを追加
  - `correlation_analysis.py` により、SE (CPM) と遺伝子 (TPM) のピアソン相関を計算
- **ドキュメントの更新**:
  - `/cli_tools/` ディレクトリに英語版 (`README.md`) と日本語版 (`README_ja.md`) のドキュメントを新設
  - リポジトリ全体の `README.md` および `/SE_to_gene_links/` のドキュメントを更新し、新しい CLI ツールとの連携方法を追加

### Changed
- `SEgene_package` から独立した CLI スクリプトを `cli_tools/` ディレクトリに整理し、リポジトリ構造を調整

### Notes
- 新規ドキュメントには、英語・日本語の両方でワークフローや使用方法を記載  
- 既存の `SE_to_gene_links` 機能と連携して解析可能

## [1.0.0] - 2025-01-15

### Added
- 初回リリース: Super-Enhancer と遺伝子のリンクを統計的手法で解析するプラットフォーム
  - `peak_to_gene_links`:
    - 遺伝子発現量とエンハンサーピークの相関情報を取得
    - データ解析と可視化支援
  - `SE_to_gene_links`:
    - Super-Enhancer と遺伝子のリンクを相関データで評価・解析
    - グラフ理論を用いた可視化
    - Jupyter Notebook での対話的解析
- GitHub リポジトリ内の詳細なインストール手順とドキュメント

[1.2.0]: https://github.com/hamamoto-lab/SEgene/releases/tag/v1.2.0
[1.1.0]: https://github.com/hamamoto-lab/SEgene/releases/tag/v1.1.0
[1.0.0]: https://github.com/hamamoto-lab/SEgene/releases/tag/v1.0.0
