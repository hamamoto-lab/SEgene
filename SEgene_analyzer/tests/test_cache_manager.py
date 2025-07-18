"""
Test suite for cache manager module
"""
import unittest
import sys
import os
import tempfile
import shutil
import pickle
import json
from unittest.mock import MagicMock, patch

# Add the parent directory to the path to import the module
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from segene_analyzer.cache_manager import SEdbCacheManager

class TestSEdbCacheManager(unittest.TestCase):
    """Test cases for SEdbCacheManager"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.temp_dir = tempfile.mkdtemp()
        self.cache_manager = SEdbCacheManager()
        
        # Create test files
        self.test_bed_file = os.path.join(self.temp_dir, "test.bed")
        self.test_sample_file = os.path.join(self.temp_dir, "test_sample.txt")
        self.cache_file = os.path.join(self.temp_dir, "test_cache.pkl")
        
        with open(self.test_bed_file, 'w') as f:
            f.write("chr1\t1000\t2000\tregion1\n")
        with open(self.test_sample_file, 'w') as f:
            f.write("sample_id\ttissue\tcell_type\n")
    
    def tearDown(self):
        """Clean up test fixtures"""
        shutil.rmtree(self.temp_dir)
    
    def test_calculate_file_hash(self):
        """Test file hash calculation"""
        hash1 = self.cache_manager.calculate_file_hash(self.test_bed_file)
        hash2 = self.cache_manager.calculate_file_hash(self.test_bed_file)
        
        # Same file should produce same hash
        self.assertEqual(hash1, hash2)
        self.assertIsInstance(hash1, str)
        self.assertEqual(len(hash1), 32)  # MD5 hash length
    
    def test_save_statistics(self):
        """Test statistics saving"""
        # Create mock analyzer
        mock_analyzer = MagicMock()
        mock_analyzer.sample_info = {'sample1': 'data1', 'sample2': 'data2'}
        mock_analyzer.bed_df = [1, 2, 3, 4, 5]  # Mock data with length
        mock_analyzer._tissue_distribution = {'tissue1': 10, 'tissue2': 20}
        mock_analyzer._cell_id_counts = {'cell1': 5, 'cell2': 15}
        
        # Save statistics
        metadata = self.cache_manager.save_statistics(
            mock_analyzer, 
            self.cache_file, 
            self.test_bed_file, 
            self.test_sample_file
        )
        
        # Check that cache file was created
        self.assertTrue(os.path.exists(self.cache_file))
        
        # Check metadata
        self.assertIsInstance(metadata, dict)
        self.assertEqual(metadata['total_samples'], 2)
        self.assertEqual(metadata['total_se_records'], 5)
        self.assertIn('created_at', metadata)
        self.assertIn('bed_file_hash', metadata)
        self.assertIn('sample_info_hash', metadata)
        
        # Check that metadata JSON file was created
        metadata_json_file = self.cache_file.replace('.pkl', '_metadata.json')
        self.assertTrue(os.path.exists(metadata_json_file))
        
        # Verify the cached data can be loaded
        with open(self.cache_file, 'rb') as f:
            cached_data = pickle.load(f)
        
        self.assertIn('metadata', cached_data)
        self.assertIn('statistics', cached_data)
        self.assertEqual(cached_data['statistics']['tissue_distribution'], {'tissue1': 10, 'tissue2': 20})
        self.assertEqual(cached_data['statistics']['cell_id_counts'], {'cell1': 5, 'cell2': 15})
    
    def test_load_statistics(self):
        """Test statistics loading"""
        # First save some statistics
        mock_analyzer = MagicMock()
        mock_analyzer.sample_info = {'sample1': 'data1'}
        mock_analyzer.bed_df = [1, 2, 3]
        mock_analyzer._tissue_distribution = {'tissue1': 10}
        mock_analyzer._cell_id_counts = {'cell1': 5}
        
        self.cache_manager.save_statistics(
            mock_analyzer, 
            self.cache_file, 
            self.test_bed_file, 
            self.test_sample_file
        )
        
        # Create new analyzer for loading
        new_analyzer = MagicMock()
        
        # Load statistics
        success = self.cache_manager.load_statistics(
            new_analyzer,
            self.cache_file,
            verify_files=True,
            bed_file=self.test_bed_file,
            sample_info_file=self.test_sample_file
        )
        
        self.assertTrue(success)
        self.assertEqual(new_analyzer._tissue_distribution, {'tissue1': 10})
        self.assertEqual(new_analyzer._cell_id_counts, {'cell1': 5})
    
    def test_load_statistics_file_verification_failure(self):
        """Test statistics loading with file verification failure"""
        # Create cache with different file
        mock_analyzer = MagicMock()
        mock_analyzer.sample_info = {'sample1': 'data1'}
        mock_analyzer.bed_df = [1, 2, 3]
        mock_analyzer._tissue_distribution = {'tissue1': 10}
        mock_analyzer._cell_id_counts = {'cell1': 5}
        
        self.cache_manager.save_statistics(
            mock_analyzer, 
            self.cache_file, 
            self.test_bed_file, 
            self.test_sample_file
        )
        
        # Try to load with different file (should fail verification)
        different_file = os.path.join(self.temp_dir, "different.bed")
        with open(different_file, 'w') as f:
            f.write("chr2\t3000\t4000\tregion2\n")
        
        new_analyzer = MagicMock()
        success = self.cache_manager.load_statistics(
            new_analyzer,
            self.cache_file,
            verify_files=True,
            bed_file=different_file,
            sample_info_file=self.test_sample_file
        )
        
        self.assertFalse(success)
    
    def test_load_statistics_without_verification(self):
        """Test statistics loading without file verification"""
        # Create cache
        mock_analyzer = MagicMock()
        mock_analyzer.sample_info = {'sample1': 'data1'}
        mock_analyzer.bed_df = [1, 2, 3]
        mock_analyzer._tissue_distribution = {'tissue1': 10}
        mock_analyzer._cell_id_counts = {'cell1': 5}
        
        self.cache_manager.save_statistics(
            mock_analyzer, 
            self.cache_file, 
            self.test_bed_file, 
            self.test_sample_file
        )
        
        # Load without verification
        new_analyzer = MagicMock()
        success = self.cache_manager.load_statistics(
            new_analyzer,
            self.cache_file,
            verify_files=False
        )
        
        self.assertTrue(success)
        self.assertEqual(new_analyzer._tissue_distribution, {'tissue1': 10})
        self.assertEqual(new_analyzer._cell_id_counts, {'cell1': 5})
    
    def test_load_statistics_nonexistent_file(self):
        """Test statistics loading with non-existent cache file"""
        new_analyzer = MagicMock()
        success = self.cache_manager.load_statistics(
            new_analyzer,
            "/non/existent/file.pkl",
            verify_files=False
        )
        
        self.assertFalse(success)
    
    def test_default_logger(self):
        """Test default logger setup"""
        cache_manager = SEdbCacheManager()
        self.assertIsNotNone(cache_manager.logger)
        self.assertEqual(cache_manager.logger.level, 20)  # INFO level

if __name__ == '__main__':
    unittest.main()