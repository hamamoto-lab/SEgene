"""
Test suite for CLI module
"""
import unittest
import sys
import os
from unittest.mock import patch, MagicMock
import tempfile
import shutil

# Add the parent directory to the path to import the module
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from segene_analyzer.cli import create_parser, validate_file_exists, setup_logging

class TestCLI(unittest.TestCase):
    """Test cases for CLI functionality"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.temp_dir = tempfile.mkdtemp()
        self.test_bed_file = os.path.join(self.temp_dir, "test.bed")
        self.test_sample_file = os.path.join(self.temp_dir, "test_sample.txt")
        
        # Create test files
        with open(self.test_bed_file, 'w') as f:
            f.write("chr1\t1000\t2000\tregion1\n")
        with open(self.test_sample_file, 'w') as f:
            f.write("sample_id\ttissue\tcell_type\n")
    
    def tearDown(self):
        """Clean up test fixtures"""
        shutil.rmtree(self.temp_dir)
    
    def test_create_parser(self):
        """Test argument parser creation"""
        parser = create_parser()
        self.assertIsNotNone(parser)
        self.assertEqual(parser.prog, 'sedb-analyzer')
        
        # Test help output contains expected commands
        help_text = parser.format_help()
        self.assertIn('prepare-stats', help_text)
        self.assertIn('batch', help_text)
        self.assertIn('single', help_text)
        self.assertIn('report', help_text)
    
    def test_validate_file_exists(self):
        """Test file validation function"""
        # Test with existing file
        result = validate_file_exists(self.test_bed_file)
        self.assertEqual(result, self.test_bed_file)
        
        # Test with non-existing file
        with self.assertRaises(FileNotFoundError):
            validate_file_exists("/non/existent/file.txt")
    
    def test_setup_logging(self):
        """Test logging setup"""
        # Test default logging
        logger = setup_logging(verbose=False)
        self.assertEqual(logger.level, 20)  # INFO level
        
        # Test verbose logging
        logger = setup_logging(verbose=True)
        self.assertEqual(logger.level, 10)  # DEBUG level
    
    def test_prepare_stats_parser(self):
        """Test prepare-stats subcommand parser"""
        parser = create_parser()
        
        # Test valid arguments
        args = parser.parse_args([
            'prepare-stats',
            '-b', self.test_bed_file,
            '-s', self.test_sample_file,
            '-o', 'output.pkl'
        ])
        
        self.assertEqual(args.command, 'prepare-stats')
        self.assertEqual(args.bed_file, self.test_bed_file)
        self.assertEqual(args.sample_info, self.test_sample_file)
        self.assertEqual(args.output, 'output.pkl')
    
    def test_batch_parser(self):
        """Test batch subcommand parser"""
        parser = create_parser()
        regions_file = os.path.join(self.temp_dir, "regions.tsv")
        with open(regions_file, 'w') as f:
            f.write("chr1\t1000\t2000\tregion1\n")
        
        # Test valid arguments
        args = parser.parse_args([
            'batch',
            '-b', self.test_bed_file,
            '-s', self.test_sample_file,
            '-r', regions_file,
            '-o', 'output_dir',
            '--max-rows', '100'
        ])
        
        self.assertEqual(args.command, 'batch')
        self.assertEqual(args.bed_file, self.test_bed_file)
        self.assertEqual(args.sample_info, self.test_sample_file)
        self.assertEqual(args.regions_file, regions_file)
        self.assertEqual(args.output_dir, 'output_dir')
        self.assertEqual(args.max_rows, 100)
    
    def test_single_parser(self):
        """Test single subcommand parser"""
        parser = create_parser()
        
        # Test valid arguments
        args = parser.parse_args([
            'single',
            '-b', self.test_bed_file,
            '-s', self.test_sample_file,
            '--chr', 'chr1',
            '--start', '1000000',
            '--end', '2000000',
            '--region-name', 'test_region'
        ])
        
        self.assertEqual(args.command, 'single')
        self.assertEqual(args.bed_file, self.test_bed_file)
        self.assertEqual(args.sample_info, self.test_sample_file)
        self.assertEqual(args.chr, 'chr1')
        self.assertEqual(args.start, 1000000)
        self.assertEqual(args.end, 2000000)
        self.assertEqual(args.region_name, 'test_region')
    
    def test_report_parser(self):
        """Test report subcommand parser"""
        parser = create_parser()
        
        # Test valid arguments
        args = parser.parse_args([
            'report',
            '-b', self.test_bed_file,
            '-s', self.test_sample_file,
            '-o', 'report_dir',
            '--top-n', '25',
            '--dpi', '300'
        ])
        
        self.assertEqual(args.command, 'report')
        self.assertEqual(args.bed_file, self.test_bed_file)
        self.assertEqual(args.sample_info, self.test_sample_file)
        self.assertEqual(args.output_dir, 'report_dir')
        self.assertEqual(args.top_n, 25)
        self.assertEqual(args.dpi, 300)
    
    def test_common_arguments(self):
        """Test common arguments across all subcommands"""
        parser = create_parser()
        
        # Test verbose flag
        args = parser.parse_args([
            'prepare-stats',
            '-b', self.test_bed_file,
            '-s', self.test_sample_file,
            '--verbose'
        ])
        self.assertTrue(args.verbose)
        
        # Test no-peak-analysis flag
        args = parser.parse_args([
            'batch',
            '-b', self.test_bed_file,
            '-s', self.test_sample_file,
            '-r', self.test_bed_file,
            '--no-peak-analysis'
        ])
        self.assertTrue(args.no_peak_analysis)

if __name__ == '__main__':
    unittest.main()