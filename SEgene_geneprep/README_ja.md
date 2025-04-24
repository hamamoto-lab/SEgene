# SEgene_geneprep

*(For the English version of this README, please see [README.md](./README.md).)*

**SEgene_geneprep** は、GTFファイルから遺伝子の位置情報を抽出し、TPMファイルの発現量データと結合するためのPythonツールです。特にnf-core/rnaseqパイプラインのsalmon出力を、SEgene解析パイプラインのRNAinput用CSVファイルに変換することを主な目的としています。

## 概要

SEgene_geneprepは、[SEgeneプロジェクト](https://github.com/hamamoto-lab/SEgene)のユーティリティツールとして、以下の役割を果たします：

1. **入力**: 
   - GTFファイル（遺伝子アノテーション）：
     - nf-core/rnaseqの`--save_reference`オプションを使用して取得できます
     - 外部から取得する場合は、[Illumina iGenomes](https://support.illumina.com/sequencing/sequencing_software/igenome.html)のゲノムアノテーション（例：GRCh38）を使用することも可能です
     - いずれの場合も、標準的なGTF形式に準拠している必要があります
   - TPMファイル（遺伝子発現量データ）：`salmon.merged.gene_tpm.tsv`などのnf-core/rnaseq出力ファイル

2. **出力**: 遺伝子ID、ゲノム位置情報（染色体、開始位置、終了位置、ストランド）、および発現量を含む統合テーブル

3. **フロー**: この出力は、SEgeneの`peak_to_gene_links`モジュールへのRNAseq入力として利用され、後続の`SE_to_gene_links`解析に繋がります


nf-core/rnaseqの詳細については、公式リポジトリ（https://github.com/nf-core/rnaseq）を参照してください。


## 要件と依存関係

### Python 環境
- Python 3.6 以上

### Python ライブラリ
- pandas（バージョン 1.5 以上を推奨）

## 使用方法

### 基本的な使用方法

```bash
python geneprep.py --gtf path/to/genes.gtf --tpm path/to/salmon.merged.gene_tpm.tsv
```

### すべてのオプション

```
usage: geneprep.py [-h] --gtf GTF --tpm TPM [--output OUTPUT] [--output-dir OUTPUT_DIR] 
                   [--unmatched UNMATCHED] [--mismatch MISMATCH] 
                   [--undefined-strand UNDEFINED_STRAND] [--join {inner,left,outer}] 
                   [--log LOG] [--no-create-dirs] [--verbose]

GTFファイルから遺伝子位置情報を抽出し、TPMデータとマージする

optional arguments:
  -h, --help            ヘルプメッセージを表示して終了
  --gtf GTF             GTFファイルのパス
  --tpm TPM             TPMファイルのパス
  --output OUTPUT       出力ファイル名（デフォルト: gene_tpm_with_position.csv）
  --output-dir OUTPUT_DIR
                        全ファイルの出力ディレクトリ（個別のファイルパスを上書き）
  --unmatched UNMATCHED マッチしなかった行の出力ファイル（デフォルト: unmatched_gtf_lines.tsv）
  --mismatch MISMATCH   gene_idとgene_nameが不一致の行の出力ファイル（デフォルト: id_name_mismatches.tsv）
  --undefined-strand UNDEFINED_STRAND
                        未定義のストランド('.')を持つ行の出力ファイル（デフォルト: undefined_strand_lines.tsv）
  --join {inner,left,outer}
                        結合方法:
                        inner（両方のファイルに存在する遺伝子のみ）,
                        left（GTFファイルに存在するすべての遺伝子）,
                        outer（どちらかのファイルに存在するすべての遺伝子）
                        （デフォルト: inner）
  --log LOG             ログファイルのパス（指定した場合、ログがこのファイルに保存されます）
  --no-create-dirs      出力ディレクトリが存在しない場合でも自動作成しない
                        （デフォルトでは自動作成されます）
  --verbose, -v         詳細出力モードを有効にする
```

### 結合方法（--join）の説明

- **inner**（内部結合）: デフォルト設定。両方のファイル（GTFとTPM）に存在する遺伝子のみを結果に含めます。
- **left**（左結合）: GTFファイルに存在するすべての遺伝子を結果に含めます。TPMファイルに対応するデータがない場合、発現量の値は欠損値（NaN）として出力されます。
- **outer**（外部結合）: どちらかのファイルに存在するすべての遺伝子を結果に含めます。片方のファイルにしか存在しない遺伝子の場合、もう片方のデータは欠損値として扱われます。

## 入力ファイル形式

### GTFファイル

nf-core/rnaseqパイプラインで使用される標準的なGTF（Gene Transfer Format）ファイル形式に対応しています。以下は例です：

```
chr1    BestRefSeq      exon    11874   12227   .       +       .       gene_id "DDX11L1"; gene_name "DDX11L1"; transcript_id "rna0"; tss_id "TSS31672";
chr1    BestRefSeq      exon    12613   12721   .       +       .       gene_id "DDX11L1"; gene_name "DDX11L1"; transcript_id "rna0"; tss_id "TSS31672";
chr1    BestRefSeq      exon    13221   14409   .       +       .       gene_id "DDX11L1"; gene_name "DDX11L1"; transcript_id "rna0"; tss_id "TSS31672";
```

各エクソンから遺伝子単位で位置情報が統合されます。

**GTFファイルの取得方法**:
- **nf-core/rnaseq経由**: パイプライン実行時に`--save_reference`オプションを使用すると、使用されたレファレンスGTFファイルが保存されます
- **Illumina iGenomes経由**: 標準的なゲノムアノテーション（例：GRCh38、GRCm39など）は、[Illumina iGenomes](https://support.illumina.com/sequencing/sequencing_software/igenome.html)から取得可能です
- **注意**: どの方法で取得した場合でも、**exonフィーチャーが含まれている**GTFファイルであることを確認してください

**重要な仕様**:
- 処理対象: **exonフィーチャーのみ**（CDS、UTR等は除外）
- ストランド情報（7列目）が'.'（未定義）の行はスキップされます。これらの行は別途 `undefined_strand_lines.tsv` ファイルに記録されます
- gene_idとgene_nameの属性が必須となります

### TPMファイル

nf-core/rnaseqパイプラインの出力である、タブ区切りのTPM（Transcripts Per Million）データファイル。基本的に`salmon.merged.gene_tpm.tsv`を対象としていますが、同様の形式であれば他のsalmon出力ファイルも変換可能です。以下は例です：

```
gene_id gene_name       SRX2370497
A1BG    A1BG    0.067802
A1BG-AS1        A1BG-AS1        1.100025
A1CF    A1CF    0
```

`gene_id`カラムは必須で、これをキーにGTFファイルとの結合が行われます。

## 出力ファイル

### 主出力ファイル（デフォルト: gene_tpm_with_position.csv）

```
symbol,chrom,strand,start,end,サンプル1,サンプル2,...
A1BG,chr19,+,58345178,58353492,0.067802,0.123456,...
A1BG-AS1,chr19,-,58353066,58358438,1.100025,2.345678,...
```

カラム説明:
- **symbol**: 遺伝子ID（GTFファイルのgene_id）
- **chrom**: 染色体名
- **strand**: 転写方向（+, -）※ストランドが未定義（'.'）の行は出力に含まれません
- **start**: 開始位置（0-based）
- **end**: 終了位置
- **サンプル1, サンプル2, ...**: TPMファイルから取得した各サンプルの発現量

### その他の出力ファイル

- **unmatched_gtf_lines.tsv**: GTFファイルでパースに失敗した行
- **id_name_mismatches.tsv**: gene_idとgene_nameが一致しない遺伝子のリスト
- **undefined_strand_lines.tsv**: ストランドが未定義（'.'）の行のリスト

## 実行例

### nf-core/rnaseq出力を使用した基本的な実行

```bash
python geneprep.py --gtf /path/to/nfcore/references/genes.gtf --tpm /path/to/nfcore/results/salmon.merged.gene_tpm.tsv
```

### 出力ディレクトリの一括指定

```bash
# 全ての出力ファイルを同じディレクトリに保存
python geneprep.py --gtf genes.gtf --tpm salmon.merged.gene_tpm.tsv --output-dir /path/to/output
```

### 詳細出力モードとログファイル出力

```bash
python geneprep.py --gtf genes.gtf --tpm salmon.merged.gene_tpm.tsv --verbose --log processing.log
```

### 外部結合を使用した実行（GTFにはあるがTPMにない遺伝子も含める）

```bash
python geneprep.py --gtf genes.gtf --tpm salmon.merged.gene_tpm.tsv --join outer --output all_genes.csv
```

### ストランド未定義の行の出力先を変更

```bash
python geneprep.py --gtf genes.gtf --tpm salmon.merged.gene_tpm.tsv --undefined-strand undefined_strands.txt
```

## 詳細統計の表示

`--verbose`オプションを使用すると、以下のような詳細な統計情報が表示されます：

```
統計情報:
GTFの遺伝子数: 25631
TPMファイルの遺伝子数: 23486
マッチした遺伝子数: 22978
GTFにあってTPMにない遺伝子数: 2653
TPMにあってGTFにない遺伝子数: 508

詳細統計:
染色体分布（GTF内、上位10）:
  chr1: 2103 遺伝子
  chr2: 1283 遺伝子
  chr3: 1078 遺伝子
  ...

遺伝子長の統計:
  平均: 32456.1 bp
  中央値: 16782.0 bp
  最小: 102.0 bp
  最大: 2304562.0 bp
```

## トラブルシューティング

よくある問題と解決策：

1. **ファイルが見つからないエラー**:
   - GTFやTPMファイルのパスが正しいか確認してください
   - ファイル権限が適切に設定されているか確認してください

2. **ストランド未定義の行が多い**:
   - `undefined_strand_lines.tsv`ファイルを確認して、どの遺伝子がスキップされたかを確認してください
   - 必要に応じて、別のGTFファイルを使用することを検討してください

3. **メモリ不足エラー**:
   - 大きなGTFファイルを処理する場合は、より多くのメモリを持つマシンで実行してください

## 注意事項

- 大きなGTFファイルを処理する場合は、十分なメモリを確保してください。
- TPMファイルには必ず`gene_id`カラムが含まれている必要があります。
- **GTFファイルのexonフィーチャーのみが処理対象となります**。CDS、UTR等のフィーチャーはスキップされます。
- **ストランドが'.'（未定義）の行は処理から除外され**、別ファイルに記録されます。
- 遺伝子の位置情報はエクソンの情報から統合され、最小の開始位置と最大の終了位置が使用されます。
- **出力ディレクトリが存在しない場合、デフォルトで自動的に作成されます**。この動作は`--no-create-dirs`オプションで無効化できます。

## 引用

このツールを研究に使用する場合は、以下を引用してください：
**(論文は現在準備中です)**

## ライセンス

このプログラムはMITライセンスで公開されています。詳細については、[LICENSE](https://github.com/hamamoto-lab/SEgene/blob/main/LICENSE)をご覧ください。