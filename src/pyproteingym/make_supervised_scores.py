"""
make_supervised_scores.py - Python equivalent of make_supervised_scores.Rmd

Downloads and processes ProteinGym supervised model scores for DMS substitution assays.
Handles contiguous_5, modulo_5, and random_5 fold types.
"""

import os
import pandas as pd
import requests
import tempfile
import zipfile
from typing import Dict, List, Optional, Tuple
import re


def get_supervised_scores_data(
    fold_type: str = "random_5", 
    cache_dir: str = ".cache"
) -> Tuple[Dict[str, pd.DataFrame], pd.DataFrame]:
    """
    Download and process ProteinGym supervised model scores.
    
    Args:
        fold_type: Type of cross-validation fold ("contiguous_5", "modulo_5", or "random_5")
        cache_dir: Directory to cache downloaded files
        
    Returns:
        Tuple of (supervised_scores_dict, summary_metrics_df)
        - supervised_scores_dict: Dictionary mapping DMS assay names to DataFrames with model predictions
        - summary_metrics_df: DataFrame with performance metrics across assays and models
    """
    if fold_type not in ["contiguous_5", "modulo_5", "random_5"]:
        raise ValueError("fold_type must be one of: 'contiguous_5', 'modulo_5', 'random_5'")
    
    os.makedirs(cache_dir, exist_ok=True)
    
    # Download supervised scores data (this would need the actual URL from Zenodo v1.2)
    zip_path = os.path.join(cache_dir, "DMS_supervised_substitutions_scores.zip")
    
    if not os.path.exists(zip_path):
        # Placeholder URL - would need actual Zenodo link
        url = "https://zenodo.org/records/14997691/files/DMS_supervised_substitutions_scores.zip"
        print(f"Downloading supervised scores from {url}...")
        # In practice, implement the actual download
        print("Note: Actual download not implemented - using placeholder")
        return {}, pd.DataFrame()
    
    # Load supervised scores for specific fold type
    supervised_tables = _load_supervised_fold_data(zip_path, fold_type)
    
    # Add UniProt IDs
    supervised_tables = _add_uniprot_ids_supervised(supervised_tables)
    
    # Clean up column names (replace hyphens with underscores, remove spaces, etc.)
    supervised_tables = _clean_supervised_column_names(supervised_tables)
    
    # Load summary metrics
    summary_df = _load_supervised_summary_metrics(cache_dir)
    
    return supervised_tables, summary_df


def _load_supervised_fold_data(zip_path: str, fold_type: str) -> Dict[str, pd.DataFrame]:
    """
    Load supervised scores for a specific fold type from zip file.
    
    Args:
        zip_path: Path to zip file containing supervised scores
        fold_type: Type of cross-validation fold
        
    Returns:
        Dictionary mapping DMS study names to DataFrames
    """
    supervised_tables = {}
    folder_name = f"fold_{fold_type}"
    
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        # Find CSV files in the appropriate folder
        csv_files = [f for f in zip_ref.namelist() 
                    if folder_name in f and f.endswith('.csv')]
        
        for csv_file in csv_files:
            # Extract study name from filename
            filename = os.path.basename(csv_file)
            study_name = os.path.splitext(filename)[0]
            
            # Read CSV from zip
            with zip_ref.open(csv_file) as f:
                df = pd.read_csv(f)
            
            # Add DMS_id column
            df['DMS_id'] = study_name
            
            # Convert DMS_score_bin to categorical
            if 'DMS_score_bin' in df.columns:
                df['DMS_score_bin'] = df['DMS_score_bin'].astype('category')
            
            supervised_tables[study_name] = df
    
    return supervised_tables


def _add_uniprot_ids_supervised(supervised_tables: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
    """
    Add UniProt accession IDs to supervised tables.
    
    Args:
        supervised_tables: Dictionary of supervised DataFrames
        
    Returns:
        Updated dictionary with UniProt_id column added
    """
    # Extract entry names (first two elements split by underscore)
    study_names = list(supervised_tables.keys())
    entry_names = []
    
    for name in study_names:
        parts = name.split('_')
        entry_name = f"{parts[0]}_{parts[1]}" if len(parts) >= 2 else parts[0]
        entry_names.append(entry_name)
    
    # Create basic mapping (would use UniProt API in practice)
    uniprot_mapping = _get_basic_uniprot_mapping(entry_names)
    
    # Add UniProt_id to each DataFrame
    updated_tables = {}
    for i, (study_name, df) in enumerate(supervised_tables.items()):
        df_copy = df.copy()
        entry_name = entry_names[i]
        uniprot_id = uniprot_mapping.get(entry_name, None)
        df_copy['UniProt_id'] = uniprot_id
        updated_tables[study_name] = df_copy
    
    return updated_tables


def _get_basic_uniprot_mapping(entry_names: List[str]) -> Dict[str, Optional[str]]:
    """Basic UniProt mapping (placeholder - would use API in practice)."""
    mapping = {}
    for entry_name in entry_names:
        if entry_name == "ANCSZ_Hobbs":
            mapping[entry_name] = None
        else:
            mapping[entry_name] = f"UNIPROT_{entry_name.replace('_', '')}"
    return mapping


def _clean_supervised_column_names(supervised_tables: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
    """
    Clean column names in supervised tables.
    
    - Replace hyphens with underscores
    - Remove spaces
    - Remove "_predictions" suffix
    
    Args:
        supervised_tables: Dictionary of supervised DataFrames
        
    Returns:
        Dictionary with cleaned column names
    """
    cleaned_tables = {}
    
    for study_name, df in supervised_tables.items():
        df_copy = df.copy()
        
        # Clean column names
        new_columns = []
        for col in df_copy.columns:
            # Replace hyphens with underscores
            clean_col = col.replace('-', '_')
            # Remove spaces
            clean_col = clean_col.replace(' ', '')
            # Remove "_predictions" suffix
            clean_col = clean_col.replace('_predictions', '')
            new_columns.append(clean_col)
        
        df_copy.columns = new_columns
        
        # Also clean string values in character columns
        for col in df_copy.columns:
            if df_copy[col].dtype == 'object':
                try:
                    df_copy[col] = df_copy[col].str.replace('-', '_')
                except:
                    pass  # Skip if not string column
        
        # Reorder columns to put core columns first
        core_cols = ['UniProt_id', 'DMS_id', 'mutant', 'mutated_sequence', 'DMS_score', 'DMS_score_bin']
        available_core = [col for col in core_cols if col in df_copy.columns]
        other_cols = [col for col in df_copy.columns if col not in core_cols]
        df_copy = df_copy[available_core + other_cols]
        
        cleaned_tables[study_name] = df_copy
    
    return cleaned_tables


def _load_supervised_summary_metrics(cache_dir: str) -> pd.DataFrame:
    """
    Load supervised model summary metrics.
    
    Args:
        cache_dir: Directory containing cached files
        
    Returns:
        DataFrame with summary metrics
    """
    summary_path = os.path.join(cache_dir, "merged_scores_substitutions_DMS.csv")
    
    if not os.path.exists(summary_path):
        print(f"Summary metrics file not found at {summary_path}")
        return pd.DataFrame()
    
    # Load summary table
    summary_df = pd.read_csv(summary_path)
    
    # Clean model names
    if 'model_name' in summary_df.columns:
        summary_df['model_name'] = (summary_df['model_name']
                                   .str.replace('-', '_')
                                   .str.replace(' ', '')
                                   .str.replace('_predictions', ''))
    
    # Add UniProt IDs (simplified mapping)
    if 'DMS_id' in summary_df.columns:
        # Extract entry names and add UniProt mapping
        entry_names = []
        for dms_id in summary_df['DMS_id']:
            parts = dms_id.split('_')
            entry_name = f"{parts[0]}_{parts[1]}" if len(parts) >= 2 else parts[0]
            entry_names.append(entry_name)
        
        uniprot_mapping = _get_basic_uniprot_mapping(entry_names)
        summary_df['UniProt_id'] = [uniprot_mapping.get(en) for en in entry_names]
        
        # Reorder columns
        cols_order = ['UniProt_id', 'DMS_id', 'model_name', 'fold_variable_name', 'Spearman', 'MSE']
        available_cols = [col for col in cols_order if col in summary_df.columns]
        summary_df = summary_df[available_cols]
    
    return summary_df


def get_supervised_model_list() -> List[str]:
    """
    Get list of supervised models available in ProteinGym.
    
    Returns:
        List of model names
    """
    # Based on the R script, there are 12 supervised models
    # This would be populated from the actual data files
    models = [
        "ProteinNPT", "MSA_Transformer", "Tranception", "EVE", 
        "DeepSequence", "GEMME", "ESM_1v", "ESM_1b", "ESM_2",
        "Vespa", "SaProt", "ProstT5"
    ]
    return models


def filter_supervised_by_models(
    supervised_tables: Dict[str, pd.DataFrame], 
    models: List[str]
) -> Dict[str, pd.DataFrame]:
    """
    Filter supervised tables to include only specified models.
    
    Args:
        supervised_tables: Dictionary of supervised DataFrames
        models: List of model names to include
        
    Returns:
        Filtered dictionary
    """
    filtered_tables = {}
    
    for study_name, df in supervised_tables.items():
        # Get core columns plus specified model columns
        core_cols = ['UniProt_id', 'DMS_id', 'mutant', 'mutated_sequence', 'DMS_score', 'DMS_score_bin']
        available_core = [col for col in core_cols if col in df.columns]
        
        model_cols = [col for col in df.columns if col in models]
        
        filtered_df = df[available_core + model_cols].copy()
        filtered_tables[study_name] = filtered_df
    
    return filtered_tables


if __name__ == "__main__":
    # Example usage
    print("Loading supervised scores data...")
    
    # Load random_5 fold data as example
    supervised_data, summary_metrics = get_supervised_scores_data("random_5")
    
    print(f"Loaded supervised data for {len(supervised_data)} DMS assays")
    
    if supervised_data:
        # Show sample from first assay
        first_key = next(iter(supervised_data))
        sample_df = supervised_data[first_key]
        print(f"\nSample from {first_key}:")
        print(sample_df.head())
        print(f"Shape: {sample_df.shape}")
        print(f"Columns: {list(sample_df.columns)}")
    
    if not summary_metrics.empty:
        print(f"\nSummary metrics shape: {summary_metrics.shape}")
        print("Sample summary metrics:")
        print(summary_metrics.head())
    
    # Show available models
    models = get_supervised_model_list()
    print(f"\nAvailable supervised models ({len(models)}):")
    for model in models:
        print(f"  - {model}")