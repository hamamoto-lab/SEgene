# SEgene_package

*(For the English version of this README, please see [README.md](https://github.com/hamamoto-lab/SEgene/blob/main/SE_to_gene_links/README.md).)*

SEgene_packageは、スーパーエンハンサーと遺伝子の関係を解析するためのPythonパッケージです。このリポジトリでは、Dockerユーザーと非Dockerユーザーの両方がパッケージを利用できるよう、必要な環境構築手順を提供しています。

## 目次

- [特徴](#特徴)
- [必要要件](#必要要件)
- [インストール手順](#インストール手順)
  - [Dockerユーザー向け](#dockerユーザー向け)
  - [非Dockerユーザー向け](#非dockerユーザー向け)
- [使い方](#使い方)
- [ライセンス](#ライセンス)
- [CITATION](#CITATION)

## 特徴

- スーパーエンハンサーと遺伝子の関連性を解析
- グラフ理論を用いたデータの可視化
- Jupyter Notebookでのインタラクティブな解析

## 必要要件

本システムはDockerもしくはユーザーローカルに[miniforge3](https://github.com/conda-forge/miniforge)をインストールした環境で実行できます。
Dockerで運用する場合は、Docker単体でも動作は可能ですが、動作の簡略化のためにDocker Compose V2の使用をお勧めします。
Compose V2のインストール確認は、以下のコマンドを使用してください。

```bash
docker compose version
```

バージョン2のメッセージが表示されれば、インストール済みです。

**注意:** コマンドは `docker compose` です。`docker-compose` ではありません。

インストールは、Ubuntuの場合、以下のコマンドで実行できます。

```bash
sudo apt-get update
sudo apt-get install docker-compose-plugin
```

詳細は公式ドキュメント（[https://docs.docker.com/compose/install/linux/](https://docs.docker.com/compose/install/linux/)）を参照してください。

ローカルインストールの場合は、ユーザーローカルにminiforge3をインストールしてください。

## インストール手順

**共通の準備事項**
あなたの環境に合わせて、

`/path/to/data` フォルダを解析用データの保存フォルダとしてご用意ください。

`/path/to/work` フォルダを、解析作業を行い、中間のデータおよび解析結果を出力するフォルダとしてご用意ください。

特にDockerで動作させる場合は、これらのフォルダをあらかじめDockerにバインドさせる必要があります。

### Dockerユーザー向け

**Docker Composeを使用した起動手順です**

1. **リポジトリのクローン**

    ```bash
    git clone https://github.com/hamamoto-lab/SEgene.git
    ```

2. **ディレクトリに移動**

    ```bash
    cd SEgene/SE_to_gene_links
    ```

3. **docker-compose.ymlの編集**

    エディタで `docker-compose.yml` を開き、以下の部分を編集してください。

    ```yaml
      # Replace "/path/to/data" with your "data" folder
      - /path/to/data:/home/jovyan/data
      # Replace "/path/to/work" with your "work" folder
      - /path/to/work:/home/jovyan/work
    ```

4. **Dockerイメージのビルド**

    ```bash
    docker-compose build
    ```

5. **コンテナの起動**

    ```bash
    docker-compose up -d
    ```

6. **Jupyter Notebookへのアクセス**

    ブラウザで以下のURLにアクセスします。

    ```
    http://localhost:8888
    ```

7. **ノートブックの使用**

    `notebooks/` ディレクトリ内のノートブックを開いて作業を開始してください。

8. **コンテナの停止**

    作業が終了したら、以下のコマンドでコンテナを停止します。

    ```bash
    docker-compose down
    ```

### 非Dockerユーザー向け

1. **リポジトリのクローン**

    ```bash
    git clone https://github.com/hamamoto-lab/SEgene.git
    ```

2. **ディレクトリに移動**

    ```bash
    cd SEgene/SE_to_gene_links
    ```

3. **Conda環境の作成**

    `environment.yml` を使用して環境を作成します。

    ```bash
    conda env create -f environment.yml
    ```

4. **環境のアクティベート**

    ```bash
    conda activate se_gene
    ```

5. **パッケージのインストール**

    パッケージをインストールします。

    ```bash
    pip install .
    ```

6. **IPythonスタートアップスクリプトの設定（任意）**

    ```bash
    mkdir -p ~/.ipython/profile_default/startup/
    cp docker/startup/00-imports.py ~/.ipython/profile_default/startup/
    ```

7. **データディレクトリの作成（必要に応じて）**

    ```bash
    mkdir -p ~/data
    ```

8. **Jupyter Notebookの起動**

    ```bash
    jupyter notebook
    ```

## 使い方

- **ノートブックの使用**: `notebooks/` ディレクトリ内に[チュートリアル用のjupyter notebook](https://github.com/hamamoto-lab/SEgene/blob/main/SE_to_gene_links/notebooks/tutorial_book_ja.ipynb)があります。これを参考にして、自身のデータ解析を始めてください。
- **パッケージの機能**: SEgene_packageには、スーパーエンハンサーと遺伝子の関連解析に必要なツールが含まれています。
- **追加の解析ツール**: SEで同定された領域と遺伝子発現の相関解析を行うための[コマンドラインツール群](https://github.com/hamamoto-lab/SEgene/tree/main/cli_tools/README_ja.md)も提供しています。詳細は各ツールのドキュメントを参照してください。

## ライセンス

このプロジェクトはMITライセンスの下で公開されています。詳細については、[`LICENSE`](https://github.com/hamamoto-lab/SEgene/blob/main/LICENSE)をご覧ください。

## CITATION

このツールを研究に使用する場合は、以下のCITATIONファイルを参照してください：
[CITATION](https://github.com/hamamoto-lab/SEgene/blob/main/CITATION)
**(論文は現在投稿準備中です。)**
