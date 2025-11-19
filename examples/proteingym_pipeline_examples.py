#!/usr/bin/env python3
"""
ProteinGym Data Pipeline Examples

This script demonstrates how to use the Python equivalents of the R scripts 
from ProteinGymR to load and process ProteinGym datasets.

Based on the R scripts in rscripts/:
- make_DMS_substitutions.Rmd -> make_dms_substitutions.py
- make_AM_supplementary.Rmd -> make_alphamissense_supplementary.py  
- make_metadata.R -> make_metadata.py
- make_supervised_scores.Rmd -> make_supervised_scores.py
- make_zero_shot_substitutions.Rmd -> make_zero_shot_substitutions.py
- make_zeroshot_DMS_subs.Rmd -> make_zeroshot_dms_benchmarks.py
"""

import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from proteingympy import (
    get_dms_substitution_data,
    get_dms_metadata, 
    get_alphamissense_proteingym_data,
    get_alphamissense_summary_stats,
    get_supervised_substitution_data,
    available_supervised_models,
    get_zero_shot_substitution_data,
    available_zero_shot_models,
    get_zero_shot_benchmark_data,
    get_benchmark_summary_stats,
    get_top_models_by_metric,
    create_complete_metadata_table
)


def example_dms_substitutions():
    """
    Example 1: Load DMS substitution data (equivalent to make_DMS_substitutions.Rmd)
    """
    print("=" * 60)
    print("Example 1: DMS Substitution Data")
    print("=" * 60)
    
    # Load all 217 DMS substitution assays
    print("Loading DMS substitution data...")
    dms_data = get_dms_substitution_data(cache_dir=".cache")
    
    print(f"✓ Loaded {len(dms_data)} DMS assays")
    
    # Show sample from first assay
    if dms_data:
        first_assay = next(iter(dms_data))
        sample_df = dms_data[first_assay]
        print(f"\nSample from {first_assay}:")
        print(f"  Shape: {sample_df.shape}")
        print(f"  Columns: {list(sample_df.columns)}")
        print(f"  First 3 rows:")
        print(sample_df.head(3))
    
    # Load DMS metadata/reference file
    print(f"\nLoading DMS metadata...")
    metadata = get_dms_metadata(cache_dir=".cache")
    print(f"✓ Loaded metadata for {len(metadata)} DMS assays")
    print(f"  Metadata columns: {list(metadata.columns)}")
    

def example_alphamissense():
    """
    Example 2: Load AlphaMissense data (equivalent to make_AM_supplementary.Rmd)
    """
    print("\n" + "=" * 60)
    print("Example 2: AlphaMissense Supplementary Data")
    print("=" * 60)
    
    # Load AlphaMissense pathogenicity scores for ProteinGym variants
    print("Loading AlphaMissense ProteinGym data...")
    try:
        am_data = get_alphamissense_proteingym_data(cache_dir=".cache")
        
        # Get summary statistics
        stats = get_alphamissense_summary_stats(am_data)
        
        print("✓ AlphaMissense data loaded successfully")
        print(f"  Total variants: {stats['total_variants']:,}")
        print(f"  DMS assays: {stats['unique_dms_assays']}")
        print(f"  Proteins: {stats['unique_proteins']}")
        print(f"  AlphaMissense score range: {stats['alphamissense_min']:.3f} - {stats['alphamissense_max']:.3f}")
        
    except Exception as e:
        print(f"⚠ AlphaMissense data not available: {e}")


def example_supervised_scores():
    """
    Example 3: Load supervised model scores (equivalent to make_supervised_scores.Rmd)
    """
    print("\n" + "=" * 60)
    print("Example 3: Supervised Model Scores")
    print("=" * 60)
    
    # Load supervised scores for random_5 fold
    print("Loading supervised model scores (random_5 fold)...")
    try:
        supervised_data, summary_metrics = get_supervised_substitution_data("random_5", cache_dir=".cache")
        
        print(f"✓ Loaded supervised data for {len(supervised_data)} DMS assays")
        
        if supervised_data:
            first_assay = next(iter(supervised_data))
            sample_df = supervised_data[first_assay]
            print(f"  Sample assay: {first_assay}")
            print(f"  Shape: {sample_df.shape}")
            print(f"  Model columns: {len(sample_df.columns) - 6}")  # Subtract core columns
        
        if not summary_metrics.empty:
            print(f"  Summary metrics shape: {summary_metrics.shape}")
        
        # Show available models
        models = available_supervised_models()
        print(f"  Available models ({len(models)}): {', '.join(models[:5])}...")
        
    except Exception as e:
        print(f"⚠ Supervised data not available: {e}")


def example_zero_shot_scores():
    """
    Example 4: Load zero-shot model scores (equivalent to make_zero_shot_substitutions.Rmd)
    """
    print("\n" + "=" * 60)
    print("Example 4: Zero-Shot Model Scores")
    print("=" * 60)
    
    print("Loading zero-shot model scores...")
    try:
        zeroshot_data = get_zero_shot_substitution_data(cache_dir=".cache")
        
        print(f"✓ Loaded zero-shot data for {len(zeroshot_data)} DMS assays")
        
        if zeroshot_data:
            first_assay = next(iter(zeroshot_data))
            sample_df = zeroshot_data[first_assay]
            print(f"  Sample assay: {first_assay}")
            print(f"  Shape: {sample_df.shape}")
            print(f"  Model columns: {len(sample_df.columns) - 6}")  # Subtract core columns
        
        # Show available models
        models = available_zero_shot_models()
        print(f"  Available models ({len(models)}): {', '.join(models[:5])}...")
        
    except Exception as e:
        print(f"⚠ Zero-shot data not available: {e}")


def example_zero_shot_benchmarks():
    """
    Example 5: Load zero-shot benchmarks (equivalent to make_zeroshot_DMS_subs.Rmd)
    """
    print("\n" + "=" * 60)
    print("Example 5: Zero-Shot Benchmarking Metrics")
    print("=" * 60)
    
    print("Loading zero-shot benchmark data...")
    try:
        benchmark_data = get_zero_shot_benchmark_data(cache_dir=".cache")
        
        # Get summary stats
        stats = get_benchmark_summary_stats(benchmark_data)
        
        print(f"✓ Loaded benchmark data")
        print(f"  Metrics: {stats.get('num_metrics', 0)}")
        print(f"  DMS assays: {stats.get('num_assays', 0)}")
        print(f"  Models: {stats.get('num_models', 0)}")
        print(f"  Available metrics: {', '.join(stats.get('metrics_available', []))}")
        
        # Show top models for Spearman correlation
        if "Spearman" in benchmark_data and not benchmark_data["Spearman"].empty:
            top_models = get_top_models_by_metric(benchmark_data, "Spearman", 5)
            print(f"\n  Top 5 models by Spearman correlation:")
            for model, score in top_models.items():
                print(f"    {model}: {score:.4f}")
        
    except Exception as e:
        print(f"⚠ Benchmark data not available: {e}")


def example_metadata():
    """
    Example 6: Generate metadata (equivalent to make_metadata.R)
    """
    print("\n" + "=" * 60)
    print("Example 6: Dataset Metadata")
    print("=" * 60)
    
    print("Creating ProteinGym metadata table...")
    
    # Create comprehensive metadata table
    metadata_df = create_complete_metadata_table()
    
    print(f"✓ Created metadata for {len(metadata_df)} datasets")
    print(f"  Columns: {list(metadata_df.columns)}")
    
    print(f"\n  Sample metadata entries:")
    print(metadata_df[['Title', 'SourceType', 'DataProvider']].head())


def main():
    """Run all examples."""
    print("ProteinGym Data Pipeline Examples")
    print("Python equivalents of ProteinGymR scripts")
    print()
    
    # Note: Most examples will show structure but not actually download large datasets
    # since the actual URLs and file structures would need to be verified
    
    example_dms_substitutions()
    example_alphamissense() 
    example_supervised_scores()
    example_zero_shot_scores()
    example_zero_shot_benchmarks()
    example_metadata()
    
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    print()
    print("All functions are available through: from proteingympy import <function_name>")
    print()
    print("Note: To actually download data, ensure you have proper URLs and network access.")
    print("The current implementations include placeholder URLs that would need to be updated")
    print("with the actual ProteinGym Zenodo repository links.")


if __name__ == "__main__":
    main()