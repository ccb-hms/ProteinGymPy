"""
make_zeroshot_dms_benchmarks.py - Python equivalent of make_zeroshot_DMS_subs.Rmd

Downloads and processes ProteinGym zero-shot benchmarking metrics for DMS substitutions.
Handles 5 performance metrics: Spearman, AUC, MCC, NDCG, and Top_recall.
"""

import os
import pandas as pd
import requests
import tempfile
import zipfile
from typing import Dict, List, Optional, Any
import numpy as np


def get_zero_shot_metrics(cache_dir: str = ".cache") -> Dict[str, pd.DataFrame]:
    """
    Download and process ProteinGym zero-shot benchmarking metrics.
    
    This loads performance metrics for zero-shot models across 217 DMS assays.
    The benchmarking uses 5 metrics to evaluate model performance in predicting
    experimental DMS measurements without training on the specific assay labels.
    
    Metrics included:
    1. Spearman's rank correlation coefficient (primary metric)
    2. Area Under the ROC Curve (AUC) 
    3. Matthews Correlation Coefficient (MCC) for bimodal measurements
    4. Normalized Discounted Cumulative Gains (NDCG) for identifying top variants
    5. Top K Recall (top 10% of DMS values)
    
    Args:
        cache_dir: Directory to cache downloaded files
        
    Returns:
        Dictionary with 5 entries (one per metric), each containing a DataFrame with:
        - Rows: 217 DMS assays
        - Columns: Model performance scores (79 models in v1.2)
    """
    os.makedirs(cache_dir, exist_ok=True)
    
    # Option 1: Load from GitHub (older approach with 62 models)
    # benchmark_data = _load_from_github()
    
    # Option 2: Load from Zenodo v1.2 (79 models)
    benchmark_data = _load_from_zenodo_v12(cache_dir)
    
    return benchmark_data


def _load_from_github() -> Dict[str, pd.DataFrame]:
    """
    Load benchmark data directly from ProteinGym GitHub repository.
    This loads the original 62-model benchmark data.
    
    Returns:
        Dictionary with 5 DataFrames for each metric
    """
    base_url = "https://raw.githubusercontent.com/OATML-Markslab/ProteinGym/main/benchmarks/DMS_zero_shot/substitutions"
    
    metrics = ["AUC", "MCC", "NDCG", "Spearman", "Top_recall"]
    score_list = {}
    
    print("Loading zero-shot benchmarks from GitHub...")
    
    for metric in metrics:
        url = f"{base_url}/{metric}/DMS_substitutions_{metric}_DMS_level.csv"
        print(f"Loading {metric} scores...")
        
        try:
            df = pd.read_csv(url)
            # Clean column names
            df.columns = df.columns.str.replace('.', '_').str.replace('__', '_').str.rstrip('_')
            score_list[metric] = df
        except Exception as e:
            print(f"Error loading {metric}: {e}")
            score_list[metric] = pd.DataFrame()
    
    # Verify consistency across metrics
    if score_list:
        _verify_benchmark_consistency(score_list)
    
    return score_list


def _load_from_zenodo_v12(cache_dir: str) -> Dict[str, pd.DataFrame]:
    """
    Load benchmark data from Zenodo v1.2 repository (79 models).
    
    Args:
        cache_dir: Directory to cache downloaded files
        
    Returns:
        Dictionary with 5 DataFrames for each metric
    """
    zip_path = os.path.join(cache_dir, "DMS_benchmarks_performance.zip")
    
    if not os.path.exists(zip_path):
        # URL from ProteinGym Zenodo v1.2
        url = "https://zenodo.org/records/14997691/files/DMS_benchmark_performance.zip?download=1"
        print(f"Downloading benchmarks from {url}...")
        response = requests.get(url, stream=True)
        response.raise_for_status()
        with open(zip_path, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
        print("Download complete.")
    
    # Extract and load benchmark files
    score_list = {}
    benchmark_dir = "benchmarks/DMS_zero_shot/substitutions"
    
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        # Find DMS_level.csv files in each metric directory
        for metric in ["AUC", "MCC", "NDCG", "Spearman", "Top_recall"]:
            metric_file = f"{benchmark_dir}/{metric}/DMS_substitutions_{metric}_DMS_level.csv"
            
            #import pdb; pdb.set_trace()
            
            if metric_file in zip_ref.namelist():
                print(f"Loading {metric} scores...")
                with zip_ref.open(metric_file) as f:
                    df = pd.read_csv(f)
                
                # Clean column names
                df = _clean_benchmark_columns(df)
                score_list[metric] = df
            else:
                print(f"Warning: {metric_file} not found in zip")
                score_list[metric] = pd.DataFrame()
    
    # Verify consistency and update model names if needed
    if score_list:
        _verify_benchmark_consistency(score_list)
        score_list = _standardize_model_names(score_list)
    
    return score_list


def _clean_benchmark_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Clean column names in benchmark DataFrames.
    
    Args:
        df: Raw benchmark DataFrame
        
    Returns:
        DataFrame with cleaned column names
    """
    df_clean = df.copy()
    
    # Replace punctuation with underscores
    df_clean.columns = df_clean.columns.str.replace(r'[^\w]', '_', regex=True)
    
    # Remove spaces
    df_clean.columns = df_clean.columns.str.replace(' ', '_')
    
    # Replace multiple underscores with single underscore
    df_clean.columns = df_clean.columns.str.replace('__+', '_', regex=True)
    
    # Remove trailing underscores
    df_clean.columns = df_clean.columns.str.rstrip('_')
    
    return df_clean


def _verify_benchmark_consistency(score_list: Dict[str, pd.DataFrame]) -> bool:
    """
    Verify that all benchmark DataFrames have consistent structure.
    
    Args:
        score_list: Dictionary of benchmark DataFrames
        
    Returns:
        True if consistent, False otherwise
    """
    if not score_list:
        return False
    
    # Get reference DataFrame (first non-empty one)
    reference_df = None
    for metric, df in score_list.items():
        if not df.empty:
            reference_df = df
            break
    
    if reference_df is None:
        return False
    
    # Check consistency
    for metric, df in score_list.items():
        if df.empty:
            continue
            
        # Check same columns
        if not list(df.columns) == list(reference_df.columns):
            print(f"Warning: Column mismatch in {metric}")
            return False
        
        # Check same DMS IDs (first column)
        if not df.iloc[:, 0].equals(reference_df.iloc[:, 0]):
            print(f"Warning: DMS ID mismatch in {metric}")
            return False
    
    print("Benchmark data consistency verified.")
    return True


def _standardize_model_names(score_list: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
    """
    Standardize model names across benchmark DataFrames.
    
    Args:
        score_list: Dictionary of benchmark DataFrames
        
    Returns:
        Dictionary with standardized model names
    """
    # This would map to standard model names used elsewhere in the package
    # For now, we'll use the cleaned names as-is
    standardized_list = {}
    
    for metric, df in score_list.items():
        if df.empty:
            standardized_list[metric] = df
            continue
        
        df_std = df.copy()
        
        # Apply any specific model name standardizations here
        # For example, ensuring consistency with zero-shot model list
        
        standardized_list[metric] = df_std
    
    return standardized_list


def get_benchmark_summary_stats(score_list: Dict[str, pd.DataFrame]) -> Dict[str, Any]:
    """
    Generate summary statistics for benchmark data.
    
    Args:
        score_list: Dictionary of benchmark DataFrames
        
    Returns:
        Dictionary with summary statistics
    """
    if not score_list:
        return {}
    
    # Get basic info from first non-empty DataFrame
    sample_df = None
    for df in score_list.values():
        if not df.empty:
            sample_df = df
            break
    
    if sample_df is None:
        return {'error': 'No valid benchmark data found'}
    
    num_assays = len(sample_df)
    num_models = len(sample_df.columns) - 1  # Subtract DMS_id column
    
    stats = {
        'num_metrics': len([df for df in score_list.values() if not df.empty]),
        'num_assays': num_assays,
        'num_models': num_models,
        'metrics_available': list(score_list.keys()),
        'model_names': list(sample_df.columns[1:]) if num_models > 0 else []
    }
    
    # Add basic statistics for each metric
    for metric, df in score_list.items():
        if df.empty:
            continue
        
        # Get numeric columns (skip DMS_id column)
        numeric_cols = df.select_dtypes(include=[np.number]).columns
        if len(numeric_cols) > 0:
            stats[f'{metric}_mean'] = df[numeric_cols].mean().mean()
            stats[f'{metric}_std'] = df[numeric_cols].std().mean()
    
    return stats


def get_top_models_by_metric(
    score_list: Dict[str, pd.DataFrame], 
    metric: str = "Spearman", 
    top_k: int = 10
) -> pd.Series:
    """
    Get top-performing models for a specific metric.
    
    Args:
        score_list: Dictionary of benchmark DataFrames
        metric: Metric to use for ranking ("Spearman", "AUC", etc.)
        top_k: Number of top models to return
        
    Returns:
        Series with top models and their average performance
    """
    if metric not in score_list or score_list[metric].empty:
        return pd.Series()
    
    df = score_list[metric]
    
    # Get numeric columns (skip DMS_id column)
    numeric_cols = df.select_dtypes(include=[np.number]).columns
    
    if len(numeric_cols) == 0:
        return pd.Series()
    
    # Calculate mean performance across all assays
    model_means = df[numeric_cols].mean()
    
    # Return top k models
    return model_means.nlargest(top_k)


def filter_benchmarks_by_models(
    score_list: Dict[str, pd.DataFrame], 
    models: List[str]
) -> Dict[str, pd.DataFrame]:
    """
    Filter benchmark data to include only specified models.
    
    Args:
        score_list: Dictionary of benchmark DataFrames
        models: List of model names to include
        
    Returns:
        Filtered dictionary
    """
    filtered_list = {}
    
    for metric, df in score_list.items():
        if df.empty:
            filtered_list[metric] = df
            continue
        
        # Keep first column (DMS_id) plus specified models
        cols_to_keep = [df.columns[0]]  # First column
        model_cols = [col for col in df.columns[1:] if col in models]
        cols_to_keep.extend(model_cols)
        
        filtered_df = df[cols_to_keep].copy()
        filtered_list[metric] = filtered_df
    
    return filtered_list


if __name__ == "__main__":
    # Example usage
    print("Loading zero-shot benchmarking data...")
    
    # This would normally download and load the data
    print("Note: Actual download not implemented in this example")
    
    benchmark_data = {}  # Placeholder
    
    if benchmark_data:
        print(f"Loaded benchmark data for {len(benchmark_data)} metrics")
        
        # Show summary stats
        stats = get_benchmark_summary_stats(benchmark_data)
        print(f"\nBenchmark Summary:")
        for key, value in stats.items():
            if not key.endswith('_names'):
                print(f"{key}: {value}")
        
        # Show top models for Spearman correlation
        if "Spearman" in benchmark_data:
            top_models = get_top_models_by_metric(benchmark_data, "Spearman", 5)
            print(f"\nTop 5 models by Spearman correlation:")
            for model, score in top_models.items():
                print(f"  {model}: {score:.4f}")
    
    # Show what metrics are available
    metrics = ["AUC", "MCC", "NDCG", "Spearman", "Top_recall"]
    print(f"\nBenchmarking metrics available:")
    for metric in metrics:
        print(f"  - {metric}")
    
    print(f"\nMetrics explanation:")
    print(f"  - Spearman: Rank correlation (primary metric)")
    print(f"  - AUC: Area under ROC curve") 
    print(f"  - MCC: Matthews correlation coefficient")
    print(f"  - NDCG: Normalized discounted cumulative gains")
    print(f"  - Top_recall: Top 10% recall")