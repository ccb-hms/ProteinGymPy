#!/usr/bin/env python3
"""
Tests for ProteinGym data pipeline functions.

Tests the Python equivalents of the R scripts from ProteinGymR.
"""

import unittest
import tempfile
import os
import shutil
import pandas as pd
from unittest.mock import patch, MagicMock, mock_open
import sys

# Add src to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from pyproteingym import (
    get_dms_substitution_data,
    get_dms_metadata,
    get_alphamissense_proteingym_data,
    get_alphamissense_summary_stats,
    get_supervised_scores_data,
    get_supervised_model_list,
    get_zero_shot_scores_data,
    get_zero_shot_model_list,
    get_zero_shot_benchmark_data,
    get_benchmark_summary_stats,
    get_top_models_by_metric,
    create_complete_metadata_table,
    save_metadata_csv
)


class TestDataPipelines(unittest.TestCase):
    """Test suite for ProteinGym data pipeline functions."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        self.cache_dir = os.path.join(self.temp_dir, "cache")
        os.makedirs(self.cache_dir, exist_ok=True)
    
    def tearDown(self):
        """Clean up test fixtures."""
        shutil.rmtree(self.temp_dir)
    
    def test_get_supervised_model_list(self):
        """Test that supervised model list is returned correctly."""
        models = get_supervised_model_list()
        
        self.assertIsInstance(models, list)
        self.assertGreater(len(models), 0)
        self.assertIn("ProteinNPT", models)
        self.assertIn("MSA_Transformer", models)
        
        # Check no duplicates
        self.assertEqual(len(models), len(set(models)))
    
    def test_get_zero_shot_model_list(self):
        """Test that zero-shot model list is returned correctly."""
        models = get_zero_shot_model_list()
        
        self.assertIsInstance(models, list)
        self.assertGreater(len(models), 0)
        self.assertIn("ESM_2_t48_15B_UR50D", models)
        self.assertIn("AlphaMissense", models)
        
        # Check no duplicates
        self.assertEqual(len(models), len(set(models)))
    
    def test_create_complete_metadata_table(self):
        """Test metadata table creation."""
        metadata_df = create_complete_metadata_table()
        
        self.assertIsInstance(metadata_df, pd.DataFrame)
        self.assertGreater(len(metadata_df), 0)
        
        # Check required columns
        required_cols = ['Title', 'Description', 'SourceType', 'SourceUrl', 'DataProvider']
        for col in required_cols:
            self.assertIn(col, metadata_df.columns)
        
        # Check for expected datasets
        titles = metadata_df['Title'].tolist()
        self.assertTrue(any('AlphaMissense' in title for title in titles))
        self.assertTrue(any('DMS' in title for title in titles))
    
    def test_save_metadata_csv(self):
        """Test saving metadata to CSV."""
        output_path = os.path.join(self.temp_dir, "test_metadata.csv")
        
        save_metadata_csv(output_path)
        
        self.assertTrue(os.path.exists(output_path))
        
        # Read back and verify
        df = pd.read_csv(output_path)
        self.assertGreater(len(df), 0)
        self.assertIn('Title', df.columns)
    
    def test_get_alphamissense_summary_stats(self):
        """Test AlphaMissense summary statistics calculation."""
        # Create mock DataFrame
        test_data = pd.DataFrame({
            'DMS_id': ['TEST1_assay', 'TEST2_assay', 'TEST1_assay'],
            'Uniprot_ID': ['P12345', 'P67890', 'P12345'],
            'variant_id': ['A1V', 'G2D', 'L3F'],
            'AlphaMissense': [0.1, 0.8, 0.3]
        })
        
        stats = get_alphamissense_summary_stats(test_data)
        
        self.assertIsInstance(stats, dict)
        self.assertEqual(stats['total_variants'], 3)
        self.assertEqual(stats['unique_dms_assays'], 2)
        self.assertEqual(stats['unique_proteins'], 2)
        self.assertAlmostEqual(stats['alphamissense_mean'], 0.4, places=2)
    
    @patch('pyproteingym.make_dms_substitutions.requests.get')
    @patch('pyproteingym.make_dms_substitutions.zipfile.ZipFile')
    def test_get_dms_substitution_data_cached(self, mock_zipfile, mock_requests):
        """Test DMS data loading with cached file."""
        # Create a fake zip file
        zip_path = os.path.join(self.cache_dir, "DMS_ProteinGym_substitutions.zip")
        with open(zip_path, 'w') as f:
            f.write("fake zip")
        
        # Mock zipfile behavior
        mock_zip_instance = MagicMock()
        mock_zipfile.return_value.__enter__.return_value = mock_zip_instance
        mock_zip_instance.namelist.return_value = ['test_assay.csv']
        
        # Mock CSV data
        csv_data = "mutant,mutated_sequence,DMS_score,DMS_score_bin\nA1V,MVKL...,1.2,1\n"
        mock_zip_instance.open.return_value.__enter__.return_value = mock_open(read_data=csv_data).return_value
        
        # Test function
        result = get_dms_substitution_data(cache_dir=self.cache_dir)
        
        # Verify no download occurred (cached)
        mock_requests.assert_not_called()
        
        # Verify structure
        self.assertIsInstance(result, dict)
    
    def test_get_benchmark_summary_stats_empty(self):
        """Test benchmark summary stats with empty data."""
        empty_benchmarks = {}
        stats = get_benchmark_summary_stats(empty_benchmarks)
        
        self.assertIsInstance(stats, dict)
        self.assertEqual(len(stats), 0)
    
    def test_get_benchmark_summary_stats_with_data(self):
        """Test benchmark summary stats with sample data."""
        # Create mock benchmark data
        sample_data = pd.DataFrame({
            'DMS_id': ['ASSAY1', 'ASSAY2'],
            'Model1': [0.5, 0.7],
            'Model2': [0.6, 0.8]
        })
        
        benchmarks = {
            'Spearman': sample_data,
            'AUC': sample_data.copy()
        }
        
        stats = get_benchmark_summary_stats(benchmarks)
        
        self.assertIsInstance(stats, dict)
        self.assertEqual(stats['num_metrics'], 2)
        self.assertEqual(stats['num_assays'], 2)
        self.assertEqual(stats['num_models'], 2)
        self.assertIn('Spearman', stats['metrics_available'])
    
    def test_get_top_models_by_metric(self):
        """Test getting top models by performance metric."""
        # Create mock benchmark data
        sample_data = pd.DataFrame({
            'DMS_id': ['ASSAY1', 'ASSAY2'],
            'ModelA': [0.5, 0.7],
            'ModelB': [0.9, 0.8],
            'ModelC': [0.3, 0.4]
        })
        
        benchmarks = {'Spearman': sample_data}
        
        top_models = get_top_models_by_metric(benchmarks, 'Spearman', top_k=2)
        
        self.assertIsInstance(top_models, pd.Series)
        self.assertEqual(len(top_models), 2)
        
        # ModelB should be first (highest average: 0.85)
        self.assertEqual(top_models.index[0], 'ModelB')
    
    def test_get_top_models_by_metric_missing(self):
        """Test getting top models for missing metric."""
        benchmarks = {'AUC': pd.DataFrame()}
        
        top_models = get_top_models_by_metric(benchmarks, 'Spearman', top_k=5)
        
        self.assertIsInstance(top_models, pd.Series)
        self.assertEqual(len(top_models), 0)
    
    @patch('pyproteingym.make_supervised_scores.zipfile.ZipFile')
    def test_get_supervised_scores_data_structure(self, mock_zipfile):
        """Test supervised scores data structure without actual download."""
        # Mock zipfile to return empty results for now
        mock_zip_instance = MagicMock()
        mock_zipfile.return_value.__enter__.return_value = mock_zip_instance
        mock_zip_instance.namelist.return_value = []
        
        supervised_data, summary_df = get_supervised_scores_data("random_5", cache_dir=self.cache_dir)
        
        self.assertIsInstance(supervised_data, dict)
        self.assertIsInstance(summary_df, pd.DataFrame)
    
    def test_invalid_fold_type(self):
        """Test invalid fold type raises error."""
        with self.assertRaises(ValueError):
            get_supervised_scores_data("invalid_fold", cache_dir=self.cache_dir)
    
    @patch('pyproteingym.make_alphamissense_supplementary.requests.get')
    def test_get_alphamissense_proteingym_data_download_error(self, mock_requests):
        """Test AlphaMissense data download error handling."""
        # Mock failed request
        mock_requests.side_effect = Exception("Network error")
        
        with self.assertRaises(Exception):
            get_alphamissense_proteingym_data(cache_dir=self.cache_dir)
    
    def test_model_lists_consistency(self):
        """Test that model lists are consistent and reasonable."""
        supervised_models = get_supervised_model_list()
        zeroshot_models = get_zero_shot_model_list()
        
        # Should be different sets (supervised vs zero-shot)
        self.assertNotEqual(set(supervised_models), set(zeroshot_models))
        
        # Should have reasonable sizes
        self.assertGreaterEqual(len(supervised_models), 5)
        self.assertGreaterEqual(len(zeroshot_models), 10)
        
        # All should be strings
        self.assertTrue(all(isinstance(m, str) for m in supervised_models))
        self.assertTrue(all(isinstance(m, str) for m in zeroshot_models))


class TestDataPipelineIntegration(unittest.TestCase):
    """Integration tests for data pipeline functions."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
    
    def tearDown(self):
        """Clean up test fixtures."""
        shutil.rmtree(self.temp_dir)
    
    def test_metadata_generation_complete(self):
        """Test complete metadata generation workflow."""
        # Generate metadata
        metadata_df = create_complete_metadata_table()
        
        # Save to file
        output_path = os.path.join(self.temp_dir, "complete_metadata.csv")
        save_metadata_csv(output_path)
        
        # Verify file was created and has content
        self.assertTrue(os.path.exists(output_path))
        
        # Read back and verify structure
        df = pd.read_csv(output_path)
        self.assertEqual(len(df), len(metadata_df))
        self.assertEqual(list(df.columns), list(metadata_df.columns))
    
    def test_model_list_functions(self):
        """Test all model list functions work together."""
        supervised = get_supervised_model_list()
        zeroshot = get_zero_shot_model_list()
        
        # Should have different models
        overlap = set(supervised) & set(zeroshot)
        total_unique = len(set(supervised) | set(zeroshot))
        
        # Some overlap is expected (models that can do both)
        self.assertGreaterEqual(total_unique, max(len(supervised), len(zeroshot)))


if __name__ == '__main__':
    # Create test suite
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    # Add test classes
    suite.addTests(loader.loadTestsFromTestCase(TestDataPipelines))
    suite.addTests(loader.loadTestsFromTestCase(TestDataPipelineIntegration))
    
    # Run tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    # Exit with error code if tests failed
    exit(0 if result.wasSuccessful() else 1)