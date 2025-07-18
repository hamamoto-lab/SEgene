# SEgene_analyzer/batch/batch_analyzer_extended.py
from .batch_analyzer import SEdbBatchAnalyzer
from .file_converters import FileConverter

class SEdbBatchAnalyzerExtended(SEdbBatchAnalyzer):
    """
    SEdbBatchAnalyzer の拡張版
    BEDファイル対応など追加機能を提供
    """
    
    def __init__(self, analyzer):
        super().__init__(analyzer)
        self.file_converter = FileConverter(logger=self.logger)
    
    # def batch_analyze_regions_from_file(self, 
    #                                   bed_file: str,
    #                                   metadata_file: str,
    #                                   data_output_dir: str,
    #                                   report_output_dir: str, 
    #                                   regions_file: str,  # TSV または BED
    #                                   **kwargs):
    #     """
    #     TSVまたはBEDファイルから複数のゲノム領域を一括解析
        
    #     Parameters:
    #     -----------
    #     regions_file : str
    #         解析対象のゲノム領域ファイル（TSVまたはBED形式）
        
    #     その他のパラメータは batch_analyze_regions_from_tsv と同じ
    #     """
    #     # ファイル形式を判別
    #     file_format = self.file_converter.detect_file_format(regions_file)
    #     self.logger.info(f"Detected file format: {file_format}")
        
    #     if file_format == 'bed':
    #         # BED → 専用TSV形式に変換
    #         self.logger.info("Converting BED file to TSV format...")
    #         regions_tsv_file = self.file_converter.convert_bed_to_tsv(
    #             regions_file, data_output_dir
    #         )
    #     else:
    #         # すでにTSV形式
    #         regions_tsv_file = regions_file
    #         self.logger.info("Using TSV file directly")
        
    #     # 既存のTSV処理パイプラインを使用
    #     return self.batch_analyze_regions_from_tsv(
    #         bed_file=bed_file,
    #         metadata_file=metadata_file,
    #         data_output_dir=data_output_dir,
    #         report_output_dir=report_output_dir,
    #         regions_tsv_file=regions_tsv_file,
    #         **kwargs
    #     )



    def batch_analyze_regions_from_file(self, 
                                    bed_file: str,
                                    metadata_file: str,
                                    data_output_dir: str,
                                    report_output_dir: str, 
                                    regions_file: str = None,      # 新しい引数名
                                    regions_tsv_file: str = None,  # 後方互換性のため
                                    **kwargs):
        """
        TSVまたはBEDファイルから複数のゲノム領域を一括解析
        
        Parameters:
        -----------
        regions_file : str, optional
            解析対象のゲノム領域ファイル（TSVまたはBED形式）
        regions_tsv_file : str, optional
            解析対象のゲノム領域ファイル（後方互換性のため）
        
        その他のパラメータは batch_analyze_regions_from_tsv と同じ
        """
        # 引数の優先順位: regions_file > regions_tsv_file
        target_file = regions_file if regions_file is not None else regions_tsv_file
        
        if target_file is None:
            raise ValueError("Either 'regions_file' or 'regions_tsv_file' must be specified")
        
        # ファイル形式を判別
        file_format = self.file_converter.detect_file_format(target_file)
        self.logger.info(f"Detected file format: {file_format}")
        
        if file_format == 'bed':
            # BED → 専用TSV形式に変換
            self.logger.info("Converting BED file to TSV format...")
            regions_tsv_file_converted = self.file_converter.convert_bed_to_tsv(
                target_file, data_output_dir
            )
        else:
            # すでにTSV形式
            regions_tsv_file_converted = target_file
            self.logger.info("Using TSV file directly")
        
        # 既存のTSV処理パイプラインを使用
        return self.batch_analyze_regions_from_tsv(
            bed_file=bed_file,
            metadata_file=metadata_file,
            data_output_dir=data_output_dir,
            report_output_dir=report_output_dir,
            regions_tsv_file=regions_tsv_file_converted,
            **kwargs
        )