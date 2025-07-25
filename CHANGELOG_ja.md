# Changelog (日本語)

本プロジェクトの変更点は、このファイルに記載します。

形式は [Keep a Changelog](https://keepachangelog.com/ja/1.1.0/) をベースとし、
[Semantic Versioning](https://semver.org/spec/v2.0.0.html) に準拠しています。

## [1.13.1] - 2025-07-25

### Added
- **SEgene_analyzer_erna メタデータドキュメント**:
  - `SEgene_analyzer_erna/metadata/README.md` を追加 - eRNAbaseメタデータファイルのドキュメント
  - `SEgene_analyzer_erna/metadata/README_ja.md` を追加 - eRNAbaseメタデータファイルの日本語ドキュメント
  - 両ファイルにeRNAbaseへのデータソースリンクとバージョン情報（2025年7月22日時点）を含む
  - 冗長な技術的説明を削除し、必要な情報に集約してドキュメントを簡素化

## [1.13.0] - 2025-07-23

### Added
- **新しい開発版: SEgene_analyzer_erna**:
  - SEgene_region_packageから派生したeRNAbase専用解析ツール
  - 4つの専用コマンドを持つコマンドラインインターフェース: `report`, `single`, `batch`, `prepare-stats`
  - `erna-analyzer`システム全体コマンドによるpipインストール対応パッケージ
  - eRNAbase特有の複数BEDファイルとParquetメタデータ形式のサポート
  - eRNAbaseデータセット用の生物種フィルタリング機能（human/mouse）
  - 包括的な二言語対応ドキュメント（英語・日本語）
  - ゲノム領域解析のための強化された染色体フィルタリングと検証機能
  - ヘッドレス環境対応のWSL互換matplotlibバックエンド設定
  - 実eRNAbaseデータによる完全なテストスイート（858サンプル検証）
  - Fisher's exact testとFDR補正による統計的エンリッチメント解析
  - 複数出力形式: PNG, SVG, PDF, CSV, HTMLレポート
  - SEgene_analyzerのeRNAbase特化版カウンターパートとして位置づけ

### Changed
- **ドキュメントの更新**:
  - メインREADMEファイルの「開発版プログラム」セクションにSEgene_analyzer_ernaを追加
  - 新しい開発版ドキュメントへのリンクを含むよう使用方法リンクを更新

## [1.12.0] - 2025-07-18

### Added
- **新しい開発版: SEgene_analyzer**:
  - 強化されたCLIインターフェースと現代的なPythonパッケージングを持つSEgene_region_analyzerの開発版
  - 4つの専用コマンドを持つコマンドラインインターフェース: `prepare-stats`, `batch`, `single`, `report`
  - `sedb-analyzer`システム全体コマンドによるpipインストール対応パッケージ
  - バッチ処理性能向上のための高度なキャッシュシステム
  - 包括的な二言語対応ドキュメント（英語・日本語）
  - 中核機能をカバーする29テストケースのテストスイート
  - WSL環境用の強化されたmatplotlibバックエンド設定
  - 複数出力形式: PNG, SVG, PDF, EPS
  - 安定版SEgene_region_analyzerと並行する実験的開発版として位置づけ

### Changed
- **ドキュメント構造の更新**:
  - メインREADMEファイルに「開発版プログラム」セクションを追加
  - プログラム構造セクションに開発版プログラムの利用可能性を記載
  - 開発版プログラム情報の太字強調によるREADMEファイルの強化
  - 開発版ドキュメントへの直接リンクによるナビゲーション改善

## [1.11.0] - 2025-05-27

### Published
- **🎉 SEgene論文正式出版**:
  - Shinkai, N., Asada, K., Machino, H., Takasawa, K., Takahashi, S., Kouno, N., Komatsu, M., Hamamoto, R., & Kaneko, S. (2025). SEgene identifies links between super enhancers and gene expression across cell types. *npj Systems Biology and Applications*, 11(1), 49. https://doi.org/10.1038/s41540-025-00533-x
  - テストデータと補足資料をFigshareで公開: https://doi.org/10.6084/m9.figshare.28171127

### Changed
- **論文出版に伴うドキュメント更新**:
  - 全READMEファイルに正式な論文引用情報を追加
  - 全ての「論文準備中」記載を出版済み論文情報に置き換え
  - CITATIONファイルを完全な出版詳細とBibTeX形式で更新
  - 日本語READMEの引用セクションタイトルを「引用 / Citation」に統一
  - ドキュメント全体でFigshareリンクをプレースホルダーから公式DOIに更新

### Added
- **引用セクション**: 全コンポーネントに適切な引用情報が含まれるよう、不足していた引用セクションを追加

## [1.10.0] - 2025-05-02

### Added
- **SEgene_peakprepにedgeR正規化CPM機能を追加**:
  - Standard log2-CPMとBigWig方式に加え、第3の定量化方法としてedgeR正規化CPM（calcnorm CPM）を追加
  - rpy2を介してRのedgeRパッケージと連携し、ChIP-seqカウントデータの高度な正規化を実装
  - edgeRの複数の正規化手法をサポート：upperquartile（デフォルト）、TMM、RLE
  - サンプル間のCPM値に基づくオプショナルなフィルタリング機能を実装
  - 英語と日本語の包括的なドキュメントを追加
  - edgeR正規化CPM専用ガイド（cpm_calcnorm_README.md、cpm_calcnorm_README_ja.md）を作成

### Changed
- **SEgene_peakprepのコードとドキュメントを更新**:
  - CPM方式を再編成し、標準log2-CPMとedgeR正規化CPMの両方の計算をサポート
  - 正規化制御のための新しい`--calcnorm`パラメータと関連オプションを持つコマンドラインインターフェースを更新
  - 複雑なファイル命名規則との互換性を向上させるサンプル名処理を改善
  - 新機能を説明するために関連するすべてのドキュメント（README.md、README_ja.md、cpm_README.md、cpm_README_ja.md）を更新
- **依存関係**:
  - edgeR正規化CPM機能のための新しい依存関係としてR (≥4.2.0)、edgeR (≥3.40.0)、rpy2 (≥3.5.0) を追加
  - RとBioconductor環境のセットアップを説明するインストール手順を更新

## [1.9.0] - 2025-04-24

### Added
- **新コンポーネント: SEgene_geneprep**:
  - RNA-seqデータ（TPM-TSVファイル）をP2GL input用に領域情報を追加する前処理ツール
  - GTFファイルから遺伝子位置情報を抽出し、TPMファイルの発現量データと結合
  - nf-core/rnaseqパイプラインのsalmon出力を、SEgene解析パイプラインのRNAinput用CSVファイルに変換
  - 標準GTFファイルとsalmon TPMファイルを含む複数の入力形式をサポート
  - Illumina iGenomesのゲノムアノテーション（例：GRCh38）に対応
  - 英語と日本語の包括的なドキュメントを提供
- **ワークフロー可視化**:
  - SEgeneの完全なワークフローを視覚化するmermaid図を追加
  - P2GLデータ準備、P2GL相関解析、スーパーエンハンサー解析、領域評価分析の4ステッププロセスを明確に図示

### Changed
- **プログラム構造の再編成**:
  - より理解しやすいように、ワークフローを4つの主要コンポーネントに再編成
  - SEgene_peakprepとSEgene_geneprepを初期のP2GLデータ準備ステップとして位置づけ
  - 新しいワークフロー体系を反映するようドキュメント構造を更新
- **ドキュメントの更新**:
  - メインREADMEファイル（英語・日本語）の包括的な改訂
  - プログラム構造セクションをより明確に再構築
  - コンポーネント間の関係性とデータフローの説明を強化

## [1.8.0] - 2025-04-17

### Added
- **SEgene_peakprepの機能拡張**: 
  - bigWig形式のワークフローにおいて、ChIP-seqデータとInputコントロールデータの差分解析機能を追加
  - `deeptools`の`bamCompare`を利用し、ChIPシグナルとInputシグナルを比較することで濃縮領域を同定

### Changed
- **SEgene_peakprepのコードとドキュメント修正**: 
  - 新機能である差分解析の導入に伴い、`SEgene_peakprep`ディレクトリ内の関連スクリプト（`bigwig_peakprep_bamdiff.py`, `bigwig_peakprep_bamdiff_utils.py`, `bigwig_peakprep.py`等）を更新
  - あわせてドキュメント（`README_ja.md`, `README.md`, `bigwig_README_ja.md`, `bigwig_README.md`）を更新し、現在の機能と使用法を反映

## [1.7.0] - 2025-04-14

### Added
- **SEgene_peakprepにBigWig方式を追加**:
  - BigWigファイルを使用したピーク定量化の新たな代替アプローチを実装
  - BAMからBigWigへの変換とシグナル定量化のための`deeptools`連携機能
  - 複数の正規化オプション（RPGC、CPM、BPM、RPKM）を持つ`bamCoverage`機能を追加
  - 指定領域からシグナル値を抽出するための`multiBigwigSummary`処理を追加
  - シグナル値の自動log2変換機能を実装
  - BigWigファイルとゲノム座標の処理を行う包括的なユーティリティを作成
- **処理ワークフローの強化**:
  - 異なる処理ステップのための3段階ワークフローと個別スクリプトを追加
  - パイプライン全体をシームレスに実行するためのラッパースクリプトを追加
  - 個別の処理ステップを別々に実行するオプションを追加

### Changed
- **SEgene_peakprepの機能拡張**:
  - CPM方式とBigWig方式の両方を反映するようドキュメントを更新
  - 両方の処理アプローチをサポートするようコマンドラインインターフェースを拡張
  - ダウンストリーム互換性のために両方の方式間で出力形式を標準化

### Fixed
- **SEgene_peakprepの用語を修正**:
  - v1.6.1での"merge_SV.tsv"形式への言及を"merge_SE.tsv"（スーパーエンハンサー）形式に修正
  - 一貫した用語を使用するよう関連するすべてのドキュメントとコードコメントを更新
  - merge_SE.tsv形式アノテーションを指定するためのフラグ名を`--is_mergesv_format`から`--is_mergese_format`に修正

## [1.6.1] - 2025-04-11

### Added
- **SEgene_peakprepの機能強化**:
  - アノテーションファイルとしてBEDおよびmerge_SE.tsv形式を直接サポート
  - BEDおよびmerge_SE形式からSAF形式への内部自動変換機能
  - merge_SE.tsv形式のアノテーションを指定するための`--is_mergese_format`フラグを追加
  
  [注: このバージョンでは一部で誤って"merge_SV.tsv"と記載されていましたが、正しくは"merge_SE.tsv"です。これはバージョン1.7.0で修正されました。]

### Changed
- **SEgene_peakprepワークフローの最適化**:
  - より直感的なパラメータ名を持つコマンドラインインターフェースの整理
  - 実行中の可読性向上のためのログレベル出力の調整
  - 一時ファイル処理の改善
- **ドキュメントの更新**:
  - READMEファイル（英語・日本語両方）の包括的な改訂
  - 異なるアノテーションファイル形式のより明確な例を追加
  - 新しいパラメータ構造を反映するようコマンド例を更新
  - トラブルシューティングセクションの強化

## [1.6.0] - 2025-04-06

### Added
- **新コンポーネント: SEgene_peakprep**:
  - ChIP-seqデータの正規化のための新しいデータ前処理パイプライン
  - 複数サンプルのBAMファイルからのピーク領域情報の自動抽出と正規化テーブルの構築
  - `samtools flagstat`を使用した総マップ済みリード数の計算
  - `featureCounts`を使用したゲノム領域内のリードカウント
  - CPM値の計算とログ変換機能
  - 既存のワークフローコンポーネントとの統合
- **ドキュメント**:
  - SEgene_peakprepのドキュメント
  - ワークフロー説明の修正（SEgene_peakprepを初期データ準備ステップとして位置づけ）

### Changed
- プロジェクト構造の再編成（SEgene_peakprepをワークフローの最初のステップとして含める）
- SEgeneワークフローの説明の更新:
  - SEgene_peakprep → peak_to_gene_links → SE_to_gene_links → SEgene_RegionAnalyzer

## [1.5.0] - 2025-04-05

### Added
- **新コンポーネント: SEgene_RegionAnalyzer**:
  - 関心あるゲノム領域におけるスーパーエンハンサー活性を評価する新たな分析ツールを追加
  - 公共データベース（現在はSEdb 2.0）との連携機能
  - 組織特異的なスーパーエンハンサー関連性のエンリッチメント解析を実行
  - SEgene出力のTSVファイルと標準BEDフォーマット入力の両方をサポート
  - 包括的なTSVレポートと視覚的なHTML出力を生成
- **ドキュメント**:
  - SEgene_RegionAnalyzerの詳細なドキュメントを追加
  - プロジェクト構成に新コンポーネントを含めるようメインREADMEを更新

### Changed
- SEgene_RegionAnalyzerをオプショナルコンポーネントとして含むようプロジェクト構造を再編成
- 拡張された分析機能を反映するようワークフロー説明を更新

## [1.4.0] - 2025-03-19

### Added
- **ネットワーク可視化機能の強化**:
  - 詳細なグラフ描画メソッド（`draw_network_detailed`、`draw_subnetwork_detailed`、`draw_two_layer_subnetwork_detailed`）を追加
  - すべての可視化機能で複数の出力形式（SVG、PNG、EPS、PDF）をサポート
  - ノードスタイリング、エッジスタイリング、サブグラフ境界のカスタマイズオプションを追加
- **データエクスポート機能の改善**:
  - 可視化と共にグラフデータを保存する機能を追加
  - 再現性のためのログ出力機能を追加
  - 様々な形式（CSV、TSV、JSON）でのデータエクスポートをサポート
- **視覚化パラメータの拡張**:
  - ラスター形式出力の解像度（DPI）制御を追加
  - 図のサイズとスタイリングのカスタマイズオプションを追加
  - すべてのネットワーク可視化メソッド（`draw_network`、`draw_subnetwork`、`draw_two_layer_subnetwork`）にカスタムタイトルのサポートを追加

### Changed
- 既存の可視化メソッドを拡張して複数の形式オプションをサポート
- グラフレンダリングインターフェースに追加のカスタマイズパラメータを追加
- 出力ファイル命名規則等を改善

### Notes
- これらの視覚化機能の強化は、既存のネットワーク解析機能との互換性を維持するよう設計

## [1.3.0] - 2025-03-15

### Added
- **ROSE概要プロットの視覚化オプション強化**:
  - DPI、図のサイズ、出力形式をカスタマイズするパラメータを追加
  - 形式選択による図の保存機能を改善

### Fixed
- `p2gl_path_to_filter_df`関数の閾値比較を修正し、フィルタリング動作を正常化

### Notes
- ドキュメントにSEgene（「エスイージーン」と発音、英語では "S-E-gene"）の発音ガイドを追加

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
  - `cli_tools/` ディレクトリに英語版 (`README.md`) と日本語版 (`README_ja.md`) のドキュメントを新設
  - リポジトリ全体の `README.md` および `SE_to_gene_links/` のドキュメントを更新し、新しい CLI ツールとの連携方法を追加

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

[1.12.0]: https://github.com/hamamoto-lab/SEgene/compare/v1.11.0...v1.12.0
[1.11.0]: https://github.com/hamamoto-lab/SEgene/compare/v1.10.0...v1.11.0
[1.10.0]: https://github.com/hamamoto-lab/SEgene/releases/tag/v1.10.0
[1.9.0]: https://github.com/hamamoto-lab/SEgene/releases/tag/v1.9.0
[1.8.0]: https://github.com/hamamoto-lab/SEgene/releases/tag/v1.8.0
[1.7.0]: https://github.com/hamamoto-lab/SEgene/releases/tag/v1.7.0
[1.6.1]: https://github.com/hamamoto-lab/SEgene/releases/tag/v1.6.1
[1.6.0]: https://github.com/hamamoto-lab/SEgene/releases/tag/v1.6.0
[1.5.0]: https://github.com/hamamoto-lab/SEgene/releases/tag/v1.5.0
[1.4.0]: https://github.com/hamamoto-lab/SEgene/releases/tag/v1.4.0
[1.3.0]: https://github.com/hamamoto-lab/SEgene/releases/tag/v1.3.0
[1.2.0]: https://github.com/hamamoto-lab/SEgene/releases/tag/v1.2.0
[1.1.0]: https://github.com/hamamoto-lab/SEgene/releases/tag/v1.1.0
[1.0.0]: https://github.com/hamamoto-lab/SEgene/releases/tag/v1.0.0