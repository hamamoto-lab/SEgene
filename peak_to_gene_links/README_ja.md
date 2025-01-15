
# P2GL: Peak-to-Gene Links ツール

このプログラムは、ChIP-seqデータのピークとRNA-seqデータの発現量を使用して、ピークから遺伝子へのリンクを構築するためのツールです。**R 4.2.2**および**Julia 1.8.3**の環境下で動作します。また、セットアップを簡略化するためのDockerコンテナも用意されています。

---

## 必要な環境

### R ライブラリ
以下のRパッケージが必要です：
- CRANパッケージ：
  - `BiocManager`
  - `data.table`
  - `openxlsx`
  - `optparse`
  - `pbmcapply`
  - `stringr`
- Bioconductorパッケージ：
  - `GenomicRanges`
  - `rhdf5`

インストールコマンド：
```r
# CRANパッケージ
install.packages(c("BiocManager", "data.table", "openxlsx", "optparse", "pbmcapply", "stringr"))

# Bioconductorパッケージ
BiocManager::install(c("GenomicRanges", "rhdf5"))
```

### Julia ライブラリ
以下のJuliaパッケージが必要です：
- `ArgParse`
- `HDF5`
- `RData`
- `StatsBase`

インストールコマンド：
```julia
using Pkg
Pkg.add(["ArgParse", "HDF5", "RData", "StatsBase"])
```

---

## Dockerを使用したクイックスタート

このプログラムはDockerコンテナ内で動作します。以下の手順に従ってテストを実行してください。

### 1. テストデータのダウンロード
以下のテストファイルを[Figshare](url)からダウンロードしてください：
- `GSE156614_rna_tumor.csv`
- `GSE156614_ChIP_tumor.tsv`

次に、以下のローカルディレクトリを作成します：
- **データディレクトリ**: `/path/to/data`
- **出力ディレクトリ**: `/path/to/output`

ダウンロードしたファイルを`/path/to/data`に配置してください。

### 2. リポジトリのクローン
リポジトリをクローンし、適切なディレクトリに移動します：
```bash
git clone https://github.com/your-repo-name/P2GL.git
cd P2GL/SEgene_test/peak_to_gene_links
```

### 3. Dockerイメージのビルド
提供されているDockerfileを使用してDockerイメージをビルドします：
```bash
docker build -t p2gl:latest -f docker/Dockerfile .
```

### 4. Docker内でプログラムを実行
Dockerコンテナを起動し、プログラムを実行します：

1. 以下のコマンドでDockerコンテナに入ります：
```bash
docker run -it --rm \
  -v /path/to/data:/data \
  -v /path/to/output:/opt/P2GL/output \
  p2gl:latest
```

2. コンテナ内で以下のコマンドを実行してP2GLを実行します：
```bash
Rscript peakToGeneLinks.R \
  -E /data/GSE156614_rna_tumor.csv \
  -P /data/GSE156614_ChIP_tumor.tsv \
  -O sample_2000000 \
  --cores 1 \
  --txWidth 2000000
```

### 5. 出力結果を確認
成功した場合、以下のファイルが`/path/to/output`に生成されます：
- `sample_2000000.tsv`
- `sample_2000000.RData`

### 6. （オプション）可視化PDFの生成
結果を元に可視化PDFを生成するには、以下を実行してください：
```bash
Rscript peakToGeneLinksPlot.R \
  -I /opt/P2GL/output/sample_2000000 \
  -O sample_2000000_plot.pdf \
  -x 1:1000
```

この作業により、ユーザー側ローカルの`/path/to/output`フォルダに`sample_2000000_plot.pdf`が生成されます。

### 7. Dockerの終了
作業が完了したら、Dockerコンテナを終了します：
```bash
exit
```

---

## 注意事項
- Dockerがシステムにインストールされ、正しく設定されていることを確認してください。
- `/path/to/data`および`/path/to/output`は、ローカル環境に合わせて変更してください。
- 実行に成功した場合、`sample_2000000.tsv`は[Figshare](url)と同様のファイルとなります。ただし、randomseedは環境によって影響を受ける為全く同一のファイルとはなりません。あなたの環境で--seedオプションにてseedを固定することで出力を固定できますが、このseedの値もあなたの環境における一時的なものであり、他の環境で同じseedを用いても結果は微妙に変動することにご留意下さい。

---

## ライセンス
このプロジェクトはMITライセンスの下で公開されています。詳細については、[`LICENSE`ファイル](https://github.com/hamamoto-lab/SEgene/blob/main/LICENSE)をご覧ください。

---

## 引用
このツールを研究に使用する場合は、以下のCITATIONファイルを参照してください：
[CITATION](https://github.com/hamamoto-lab/SEgene/blob/main/CITATION)

