#!/usr/bin/env python3
"""
Tests for model_corr_plot module.
"""

import unittest
import sys
import os
import pandas as pd
import numpy as np

# Add src to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from proteingympy.model_corr_plot import (
    _filter_alphamissense_table,
    _filter_model_table,
    _merge_model_tables,
    _calculate_spearman_correlation,
    get_available_models
)


class TestModelCorrPlot(unittest.TestCase):
    """Test suite for model_corr_plot module."""

    def setUp(self):
        """Set up test fixtures."""
        # Create mock AlphaMissense data
        self.mock_am_data = pd.DataFrame({
            'Uniprot_ID': ['P12345', 'P12345', 'P12345'],
            'DMS_id': ['TEST_HUMAN_2023', 'TEST_HUMAN_2023', 'TEST_HUMAN_2023'],
            'variant_id': ['A1G', 'A2C', 'A3D'],
            'AlphaMissense': [0.1, 0.5, 0.9]
        })

        # Create mock model data
        self.mock_model_data = {
            'TEST_HUMAN_2023': pd.DataFrame({
                'UniProt_id': ['P12345', 'P12345', 'P12345'],
                'DMS_id': ['TEST_HUMAN_2023', 'TEST_HUMAN_2023', 'TEST_HUMAN_2023'],
                'mutant': ['A1G', 'A2C', 'A3D'],
                'GEMME': [0.2, 0.6, 0.8],
                'ESM_2_t48_15B_UR50D': [0.15, 0.55, 0.85]
            })
        }

    def test_filter_alphamissense_table(self):
        """Test filtering AlphaMissense table by UniProt ID."""
        # Test with provided table
        filtered = _filter_alphamissense_table(
            am_table=self.mock_am_data.rename(columns={
                'Uniprot_ID': 'UniProt_id',
                'variant_id': 'mutant'
            }),
            uniprot_id='P12345'
        )

        self.assertEqual(len(filtered), 3)
        self.assertTrue('UniProt_id' in filtered.columns)
        self.assertTrue('mutant' in filtered.columns)

    def test_filter_alphamissense_table_not_found(self):
        """Test filtering with non-existent UniProt ID."""
        renamed_data = self.mock_am_data.rename(columns={
            'Uniprot_ID': 'UniProt_id',
            'variant_id': 'mutant'
        })

        with self.assertRaises(ValueError) as context:
            _filter_alphamissense_table(
                am_table=renamed_data,
                uniprot_id='NOTFOUND'
            )

        self.assertIn('No AlphaMissense information found', str(context.exception))

    def test_filter_model_table(self):
        """Test filtering model table by UniProt ID."""
        filtered = _filter_model_table(
            model_table=self.mock_model_data,
            uniprot_id='P12345'
        )

        self.assertEqual(len(filtered), 3)
        self.assertTrue('UniProt_id' in filtered.columns)
        self.assertTrue('GEMME' in filtered.columns)

    def test_filter_model_table_not_found(self):
        """Test filtering model table with non-existent UniProt ID."""
        with self.assertRaises(ValueError) as context:
            _filter_model_table(
                model_table=self.mock_model_data,
                uniprot_id='NOTFOUND'
            )

        self.assertIn('No ProteinGym assay found', str(context.exception))

    def test_merge_model_tables(self):
        """Test merging two model tables."""
        model_df1 = pd.DataFrame({
            'UniProt_id': ['P12345', 'P12345', 'P12345'],
            'DMS_id': ['TEST', 'TEST', 'TEST'],
            'mutant': ['A1G', 'A2C', 'A3D'],
            'AlphaMissense': [0.1, 0.5, 0.9]
        })

        model_df2 = pd.DataFrame({
            'UniProt_id': ['P12345', 'P12345', 'P12345'],
            'DMS_id': ['TEST', 'TEST', 'TEST'],
            'mutant': ['A1G', 'A2C', 'A3D'],
            'GEMME': [0.2, 0.6, 0.8]
        })

        merged = _merge_model_tables(
            model_df1=model_df1,
            model_df2=model_df2,
            model1='AlphaMissense',
            model2='GEMME'
        )

        self.assertEqual(len(merged), 3)
        self.assertTrue('mean_model1' in merged.columns)
        self.assertTrue('mean_model2' in merged.columns)
        self.assertEqual(merged['mean_model1'].tolist(), [0.1, 0.5, 0.9])

    def test_merge_model_tables_different_uniprot(self):
        """Test merging tables with different UniProt IDs raises error."""
        model_df1 = pd.DataFrame({
            'UniProt_id': ['P12345'],
            'mutant': ['A1G'],
            'AlphaMissense': [0.1]
        })

        model_df2 = pd.DataFrame({
            'UniProt_id': ['P99999'],
            'mutant': ['A1G'],
            'GEMME': [0.2]
        })

        with self.assertRaises(ValueError) as context:
            _merge_model_tables(
                model_df1=model_df1,
                model_df2=model_df2,
                model1='AlphaMissense',
                model2='GEMME'
            )

        self.assertIn('UniProt IDs must be the same', str(context.exception))

    def test_calculate_spearman_correlation(self):
        """Test Spearman correlation calculation."""
        merged_table = pd.DataFrame({
            'mean_model1': [0.1, 0.2, 0.3, 0.4, 0.5],
            'mean_model2': [0.15, 0.25, 0.32, 0.45, 0.52]
        })

        correlation, pvalue = _calculate_spearman_correlation(merged_table)

        self.assertIsInstance(correlation, (float, np.floating))
        self.assertIsInstance(pvalue, (float, np.floating))
        self.assertGreater(correlation, 0.9)  # Should be highly correlated
        self.assertLess(pvalue, 0.1)

    def test_get_available_models(self):
        """Test getting available models."""
        models = get_available_models()

        self.assertIn('zero_shot', models)
        self.assertIn('supervised', models)
        self.assertIn('other', models)

        self.assertIsInstance(models['zero_shot'], list)
        self.assertIsInstance(models['supervised'], list)
        self.assertIsInstance(models['other'], list)

        self.assertGreater(len(models['zero_shot']), 0)
        self.assertGreater(len(models['supervised']), 0)
        self.assertIn('AlphaMissense', models['other'])

    def test_merge_with_averaging(self):
        """Test that merging averages scores across multiple studies."""
        model_df1 = pd.DataFrame({
            'UniProt_id': ['P12345', 'P12345', 'P12345', 'P12345'],
            'DMS_id': ['TEST1', 'TEST1', 'TEST2', 'TEST2'],
            'mutant': ['A1G', 'A2C', 'A1G', 'A2C'],
            'AlphaMissense': [0.1, 0.5, 0.2, 0.6]
        })

        model_df2 = pd.DataFrame({
            'UniProt_id': ['P12345', 'P12345', 'P12345', 'P12345'],
            'DMS_id': ['TEST1', 'TEST1', 'TEST2', 'TEST2'],
            'mutant': ['A1G', 'A2C', 'A1G', 'A2C'],
            'GEMME': [0.2, 0.6, 0.3, 0.7]
        })

        merged = _merge_model_tables(
            model_df1=model_df1,
            model_df2=model_df2,
            model1='AlphaMissense',
            model2='GEMME'
        )

        # Should average across studies
        self.assertEqual(len(merged), 2)  # Only 2 unique mutants

        # Check A1G average: (0.1 + 0.2) / 2 = 0.15
        a1g_row = merged[merged['mutant'] == 'A1G']
        self.assertAlmostEqual(a1g_row['mean_model1'].values[0], 0.15)
        self.assertAlmostEqual(a1g_row['mean_model2'].values[0], 0.25)


class TestModelCorrPlotIntegration(unittest.TestCase):
    """Integration tests for model_corr_plot (require actual data)."""

    @unittest.skip("Requires actual data download - run manually")
    def test_model_corr_plot_full(self):
        """Test full model_corr_plot workflow with real data."""
        from proteingympy.model_corr_plot import model_corr_plot
        import matplotlib.pyplot as plt

        # This would require actual data
        fig = model_corr_plot(
            uniprot_id="Q9NV35",
            model1="AlphaMissense",
            model2="GEMME"
        )

        self.assertIsNotNone(fig)
        plt.close(fig)


if __name__ == '__main__':
    unittest.main()
