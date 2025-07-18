# sedb_tools/batch_analyzer.py

import os
import pandas as pd
import re
import logging
from datetime import datetime
import traceback

class SEdbBatchAnalyzer:
    """
    SEdb のバッチ解析機能を担当するクラス。
    複数のゲノム領域を一括処理し、結果を統合したレポートを生成します。
    """
    
    def __init__(self, analyzer):
        """
        初期化
        
        Parameters:
        -----------
        analyzer : SEdbRegionAnalyzer
            元となるアナライザークラスのインスタンス
        """
        self.analyzer = analyzer
        self.logger = analyzer.logger or self._setup_default_logger()
        self.logger.info("Initializing SEdbBatchAnalyzer")
    
    def _setup_default_logger(self):
        """デフォルトロガーのセットアップ"""
        logger = logging.getLogger(f"{__name__}.{id(self)}")
        if not logger.handlers:
            logger.setLevel(logging.INFO)
            handler = logging.StreamHandler()
            formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
            handler.setFormatter(formatter)
            logger.addHandler(handler)
            logger.propagate = False
        return logger

    def batch_analyze_regions_from_tsv(self, 
                                    bed_file: str,  # ディレクトリではなく単一ファイル
                                    metadata_file: str,
                                    data_output_dir: str,
                                    report_output_dir: str, 
                                    regions_tsv_file: str,
                                    filter_chromosomes: bool = True,  # 染色体フィルタリングオプション
                                    count_cell_id: bool = True,       # cell_id カウントオプション
                                    max_rows: int = None,
                                    **kwargs):
        """
        ユーザー補助メソッド：TSVファイルから複数のゲノム領域を一括解析し、結果を新たなTSVとして出力します。
        
        このメソッドの動作：
        1. generate_sedb_report を実行して全体のメタデータ評価を行います
        2. 指定されたTSVファイルから読み込んだ各ゲノム領域に対して analyze_and_report_region を実行します
        3. 各領域の解析結果を元のTSVに追加カラムとして付加した新たなTSVファイルを生成します
        
        Parameters:
        -----------
        bed_file : str
            BEDファイルのパス
        metadata_file : str
            SEdbメタデータファイルのパス
        data_output_dir : str
            解析データの出力先ディレクトリ
        report_output_dir : str
            レポートファイルの出力先ディレクトリ
        regions_tsv_file : str
            解析対象のゲノム領域が記載されたTSVファイルのパス
        filter_chromosomes : bool, optional
            人間の標準染色体でフィルタリングするかどうか（デフォルト: True）
        count_cell_id : bool, optional
            cell_idの出現回数をカウントするかどうか（デフォルト: True）
        max_rows : int, optional
            TSVファイルから読み込む最大行数。Noneの場合は全行を読み込みます
        **kwargs : dict
            その他のパラメータ（analyze_and_report_regionに渡されます）
        
        Returns:
        --------
        str
            生成された結果TSVファイルのパス
        """
        
        # 出力ディレクトリが存在しない場合は作成
        os.makedirs(data_output_dir, exist_ok=True)
        os.makedirs(report_output_dir, exist_ok=True)
        
        # タイムスタンプを作成（ファイル名に使用）
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # BEDファイルとメタデータを読み込み
        self.logger.info(f"Loading BED file: {bed_file}")
        self.analyzer.set_se_bed_file(bed_file)
        self.analyzer.set_sample_info_file(metadata_file)
        self.analyzer.load_databases()
        
        # 染色体フィルタリング（オプション）
        if filter_chromosomes:
            self.logger.info("Filtering for standard human chromosomes")
            self.analyzer.filter_human_chromosomes()
        
        # cell_id出現回数カウント（オプション）
        if count_cell_id:
            self.logger.info("Counting cell_id occurrences")
            self.analyzer.count_cell_id_occurrences()
        
        # SEdbデータベース全体の解析を実行
        self.logger.info("Generating SEdb metadata report")
        db_report_path = self.analyzer.generate_sedb_report(
            output_dir=os.path.join(report_output_dir, "sedb_report")
        )
        self.logger.info(f"SEdb metadata report generated: {db_report_path}")
        
        # TSVファイルを読み込み
        self.logger.info(f"Reading regions from TSV file: {regions_tsv_file}")
        regions_df = pd.read_csv(regions_tsv_file, sep='\t')
        
        # 行数を制限（指定がある場合）
        if max_rows is not None and max_rows > 0:
            regions_df = regions_df.head(max_rows)
            self.logger.info(f"Processing first {max_rows} rows from TSV file")
        
        # 結果を格納するための新しいカラムを追加
        regions_df['has_overlaps'] = False
        regions_df['overlap_sample_count'] = 0
        regions_df['overlap_read_count'] = 0
        regions_df['significant_tissues'] = ''
        regions_df['significant_tissues_fc'] = ''  # 追加: 折り畳み変化値
        regions_df['significant_tissues_pvalue'] = ''  # 追加: p値
        regions_df['significant_tissues_fdr'] = ''  # 追加: FDR値  
        

        # 各リージョンに対して解析を実行
        for idx, row in regions_df.iterrows():
            # 最初の列をリージョン識別子として取得
            region_id = str(row.iloc[0])
            
            # リージョン文字列をパース（例: "chr1_123456_234567"）
            match = re.match(r'(chr\w+)_(\d+)_(\d+)', region_id)
            if match:
                chrom, start, end = match.groups()
                start = int(start)
                end = int(end)
                
                # リージョンの整形された名前
                region_name = f"Region_{idx+1}_{region_id}"
                region_dir = os.path.join(report_output_dir, region_name)
                
                self.logger.info(f"Analyzing region {idx+1}/{len(regions_df)}: {region_id} ({chrom}:{start}-{end})")
                
                # try:
                #     # リージョンの解析を実行
                #     results = self.analyzer.analyze_and_report_region(
                #         chrom=chrom,
                #         start=start,
                #         end=end,
                #         region_name=region_name,
                #         output_dir=region_dir,
                #         **kwargs
                #     )
                #     print("analyze_and_report_region実行")
                #     print("result:")
                #     print(results)

                try:
                    # リージョンの解析を実行
                    # extendに書き換え
                    results = self.analyzer.generate_region_report_extended(
                        chrom=chrom,
                        start=start,
                        end=end,
                        region_name=region_name,
                        output_dir=region_dir,
                        **kwargs
                    )
                    print("generate_region_report_extended実行")
                    print("result:")
                    print(results)


                    # # 結果をDataFrameに追加
                    # regions_df.at[idx, 'has_overlaps'] = results['has_overlaps']
                    # regions_df.at[idx, 'overlap_sample_count'] = results['overlap_sample_count']
                    # regions_df.at[idx, 'overlap_read_count'] = results['overlap_read_count']
                    
                    # # 有意な組織のリストをカンマ区切りの文字列に変換
                    # if results['significant_tissues']:
                    #     tissues = [item['tissue'] for item in results['significant_tissues']]
                    #     regions_df.at[idx, 'significant_tissues'] = ','.join(tissues)
                    


                    # 改善版:
                    regions_df.at[idx, 'has_overlaps'] = results.get('has_overlaps', False)
                    regions_df.at[idx, 'overlap_sample_count'] = results.get('overlap_sample_count', 0) 
                    regions_df.at[idx, 'overlap_read_count'] = results.get('overlap_read_count', 0)

                    if results.get('significant_tissues'):
                        # 組織名
                        tissues = [item.get('tissue', '') for item in results['significant_tissues']]
                        regions_df.at[idx, 'significant_tissues'] = ','.join(tissues)
                        
                        # Fold change値
                        fcs = [f"{item.get('fold_change', 0):.2f}" for item in results['significant_tissues']]
                        regions_df.at[idx, 'significant_tissues_fc'] = ','.join(fcs)
                        
                        # p値
                        pvals = [f"{item.get('pvalue', 1):.2e}" for item in results['significant_tissues']]
                        regions_df.at[idx, 'significant_tissues_pvalue'] = ','.join(pvals)
                        
                        # FDR値
                        fdrs = [f"{item.get('fdr', 1):.2e}" for item in results['significant_tissues']]
                        regions_df.at[idx, 'significant_tissues_fdr'] = ','.join(fdrs)



                    self.logger.info(f"Analysis completed for region {region_id}")
                    self.logger.info(f"Report saved to {region_dir}")
                    
                except Exception as e:
                    self.logger.error(f"Error analyzing region {region_id}: {str(e)}")
                    self.logger.error(traceback.format_exc())
            else:
                self.logger.warning(f"Row {idx+1}: Could not parse region ID '{region_id}'. Expected format: chr1_123456_234567. Skipping.")
        
        # 結果TSVファイルを生成
        output_tsv_path = os.path.join(data_output_dir, f"regions_analysis_results_{timestamp}.tsv")
        regions_df.to_csv(output_tsv_path, sep='\t', index=False)
        
        self.logger.info(f"Analysis completed for all regions")
        self.logger.info(f"Results saved to: {output_tsv_path}")
        
        return output_tsv_path