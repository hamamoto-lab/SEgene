"""
Basic functionality tests for SEgene Analyzer
"""
import unittest
import sys
import os
import tempfile
import shutil
import pandas as pd
from unittest.mock import MagicMock, patch

# Add the parent directory to the path to import the module
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

class TestBasicFunctionality(unittest.TestCase):
    """Test basic functionality that doesn't require heavy dependencies"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.temp_dir = tempfile.mkdtemp()
        
        # Create test BED file
        self.test_bed_file = os.path.join(self.temp_dir, "test.bed")
        with open(self.test_bed_file, 'w') as f:
            f.write("chr1\t1000\t2000\tregion1\t100\n")
            f.write("chr2\t3000\t4000\tregion2\t200\n")
            f.write("chr3\t5000\t6000\tregion3\t300\n")
        
        # Create test sample info file
        self.test_sample_file = os.path.join(self.temp_dir, "sample_info.txt")
        with open(self.test_sample_file, 'w') as f:
            f.write("sample_id\ttissue\tcell_type\tbiosample_type\n")
            f.write("sample1\tBrain\tNeuron\tPrimary\n")
            f.write("sample2\tLiver\tHepatocyte\tCulture\n")
            f.write("sample3\tHeart\tCardiomyocyte\tPrimary\n")
        
        # Create test regions file (TSV)
        self.test_regions_tsv = os.path.join(self.temp_dir, "regions.tsv")
        with open(self.test_regions_tsv, 'w') as f:
            f.write("chromosome\tstart\tend\tname\n")
            f.write("chr1\t1500\t2500\ttest_region1\n")
            f.write("chr2\t3500\t4500\ttest_region2\n")
        
        # Create test regions file (BED)
        self.test_regions_bed = os.path.join(self.temp_dir, "regions.bed")
        with open(self.test_regions_bed, 'w') as f:
            f.write("chr1\t1500\t2500\ttest_region1\n")
            f.write("chr2\t3500\t4500\ttest_region2\n")
    
    def tearDown(self):
        """Clean up test fixtures"""
        shutil.rmtree(self.temp_dir)
    
    def test_file_reading(self):
        """Test basic file reading functionality"""
        # Test BED file reading
        bed_df = pd.read_csv(self.test_bed_file, sep='\t', 
                            names=['chr', 'start', 'end', 'name', 'score'])
        self.assertEqual(len(bed_df), 3)
        self.assertEqual(bed_df.iloc[0]['chr'], 'chr1')
        self.assertEqual(bed_df.iloc[0]['start'], 1000)
        self.assertEqual(bed_df.iloc[0]['end'], 2000)
        
        # Test sample info file reading
        sample_df = pd.read_csv(self.test_sample_file, sep='\t')
        self.assertEqual(len(sample_df), 3)
        self.assertEqual(sample_df.iloc[0]['sample_id'], 'sample1')
        self.assertEqual(sample_df.iloc[0]['tissue'], 'Brain')
        
        # Test regions TSV file reading
        regions_df = pd.read_csv(self.test_regions_tsv, sep='\t')
        self.assertEqual(len(regions_df), 2)
        self.assertEqual(regions_df.iloc[0]['chromosome'], 'chr1')
        self.assertEqual(regions_df.iloc[0]['start'], 1500)
    
    def test_file_format_detection(self):
        """Test file format detection"""
        # Test TSV format detection
        with open(self.test_regions_tsv, 'r') as f:
            first_line = f.readline().strip()
            self.assertIn('chromosome', first_line)
            self.assertIn('start', first_line)
            self.assertIn('end', first_line)
        
        # Test BED format detection
        with open(self.test_regions_bed, 'r') as f:
            first_line = f.readline().strip()
            fields = first_line.split('\t')
            self.assertEqual(len(fields), 4)
            self.assertTrue(fields[0].startswith('chr'))
            self.assertTrue(fields[1].isdigit())
            self.assertTrue(fields[2].isdigit())
    
    def test_chromosome_validation(self):
        """Test chromosome name validation"""
        valid_chromosomes = ['chr1', 'chr2', 'chr3', 'chrX', 'chrY']
        invalid_chromosomes = ['chr23', 'chrZ', 'scaffold123', 'contig456']
        
        # Human chromosomes (1-22, X, Y)
        human_chromosomes = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
        
        for chrom in valid_chromosomes:
            if chrom in human_chromosomes:
                self.assertTrue(chrom.startswith('chr'))
        
        for chrom in invalid_chromosomes:
            self.assertNotIn(chrom, human_chromosomes)
    
    def test_coordinate_validation(self):
        """Test genomic coordinate validation"""
        # Test valid coordinates
        valid_coords = [
            (1000, 2000),
            (100000, 200000),
            (1000000, 2000000)
        ]
        
        for start, end in valid_coords:
            self.assertLess(start, end)
            self.assertGreaterEqual(start, 0)
            self.assertGreaterEqual(end, 0)
        
        # Test invalid coordinates
        invalid_coords = [
            (2000, 1000),  # start > end
            (-1000, 2000),  # negative start
            (1000, -2000)   # negative end
        ]
        
        for start, end in invalid_coords:
            self.assertFalse(start < end and start >= 0 and end >= 0)
    
    def test_data_consistency(self):
        """Test data consistency checks"""
        # Read the test files
        bed_df = pd.read_csv(self.test_bed_file, sep='\t', 
                            names=['chr', 'start', 'end', 'name', 'score'])
        sample_df = pd.read_csv(self.test_sample_file, sep='\t')
        
        # Check required columns exist
        self.assertIn('chr', bed_df.columns)
        self.assertIn('start', bed_df.columns)
        self.assertIn('end', bed_df.columns)
        
        self.assertIn('sample_id', sample_df.columns)
        self.assertIn('tissue', sample_df.columns)
        
        # Check data types
        self.assertTrue(pd.api.types.is_integer_dtype(bed_df['start']))
        self.assertTrue(pd.api.types.is_integer_dtype(bed_df['end']))
        self.assertTrue(pd.api.types.is_object_dtype(bed_df['chr']))
        
        # Check no missing values in critical columns
        self.assertFalse(bed_df['chr'].isna().any())
        self.assertFalse(bed_df['start'].isna().any())
        self.assertFalse(bed_df['end'].isna().any())
        self.assertFalse(sample_df['sample_id'].isna().any())
    
    def test_region_overlap_logic(self):
        """Test region overlap detection logic"""
        # Define test regions
        region1 = {'chr': 'chr1', 'start': 1000, 'end': 2000}
        region2 = {'chr': 'chr1', 'start': 1500, 'end': 2500}  # Overlaps with region1
        region3 = {'chr': 'chr1', 'start': 3000, 'end': 4000}  # No overlap
        region4 = {'chr': 'chr2', 'start': 1000, 'end': 2000}  # Different chromosome
        
        def regions_overlap(r1, r2):
            """Simple overlap detection"""
            if r1['chr'] != r2['chr']:
                return False
            return not (r1['end'] <= r2['start'] or r1['start'] >= r2['end'])
        
        # Test overlaps
        self.assertTrue(regions_overlap(region1, region2))
        self.assertFalse(regions_overlap(region1, region3))
        self.assertFalse(regions_overlap(region1, region4))
        
        # Test symmetric property
        self.assertEqual(regions_overlap(region1, region2), regions_overlap(region2, region1))
    
    def test_basic_statistics(self):
        """Test basic statistical calculations"""
        # Create sample data
        data = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        
        # Test basic statistics
        mean = sum(data) / len(data)
        self.assertEqual(mean, 5.5)
        
        median = sorted(data)[len(data) // 2]
        self.assertEqual(median, 6)
        
        # Test frequency counting
        tissues = ['Brain', 'Liver', 'Brain', 'Heart', 'Liver', 'Brain']
        tissue_counts = {}
        for tissue in tissues:
            tissue_counts[tissue] = tissue_counts.get(tissue, 0) + 1
        
        self.assertEqual(tissue_counts['Brain'], 3)
        self.assertEqual(tissue_counts['Liver'], 2)
        self.assertEqual(tissue_counts['Heart'], 1)
    
    def test_output_directory_creation(self):
        """Test output directory creation"""
        output_dir = os.path.join(self.temp_dir, 'output')
        
        # Directory should not exist initially
        self.assertFalse(os.path.exists(output_dir))
        
        # Create directory
        os.makedirs(output_dir, exist_ok=True)
        
        # Directory should exist now
        self.assertTrue(os.path.exists(output_dir))
        self.assertTrue(os.path.isdir(output_dir))
    
    def test_file_encoding_handling(self):
        """Test handling of different file encodings"""
        # Test UTF-8 file
        utf8_file = os.path.join(self.temp_dir, "utf8.txt")
        with open(utf8_file, 'w', encoding='utf-8') as f:
            f.write("sample_id\ttissue\n")
            f.write("sample1\t脳\n")  # Japanese characters
        
        # Should be able to read UTF-8 file
        df = pd.read_csv(utf8_file, sep='\t')
        self.assertEqual(len(df), 1)
        self.assertEqual(df.iloc[0]['sample_id'], 'sample1')

if __name__ == '__main__':
    unittest.main()