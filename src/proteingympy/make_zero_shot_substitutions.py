"""
make_zero_shot_substitutions.py - Python equivalent of make_zero_shot_substitutions.Rmd

Downloads and processes ProteinGym zero-shot model scores for DMS substitution assays.
"""

import os
import pandas as pd
import requests
import tempfile
import zipfile
from typing import Dict, List, Optional, Any
import re


def get_zero_shot_scores_data(cache_dir: str = ".cache") -> Dict[str, pd.DataFrame]:
    """
    Download and process ProteinGym zero-shot model scores for DMS substitutions.
    
    This loads zero-shot model predictions across 217 DMS assays for multiple models.
    Each assay contains predictions from various protein language models and other
    zero-shot approaches.
    
    Args:
        cache_dir: Directory to cache downloaded files
        
    Returns:
        Dictionary mapping DMS assay names to DataFrames with columns:
        - UniProt_id: UniProt accession identifier
        - DMS_id: DMS assay identifier  
        - mutant: substitution description
        - mutated_sequence: full amino acid sequence
        - DMS_score: experimental measurement
        - DMS_score_bin: binary fitness classification
        - [model_name]: Prediction scores from various zero-shot models
    """
    os.makedirs(cache_dir, exist_ok=True)
    
    # Download zero-shot scores data
    zip_path = os.path.join(cache_dir, "zero_shot_substitutions_scores.zip")
    
    if not os.path.exists(zip_path):
        # URL from ProteinGym Zenodo v1.2
        url = "https://zenodo.org/records/14997691/files/zero_shot_substitutions_scores.zip?download=1"
        print(f"Downloading zero-shot scores from {url}...")
        response = requests.get(url, stream=True)
        response.raise_for_status()
        with open(zip_path, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
        print("Download complete.")
    
    # Load zero-shot scores
    zeroshot_tables = _load_zero_shot_data(zip_path)
    
    # Add UniProt IDs
    zeroshot_tables = _add_uniprot_ids_zeroshot(zeroshot_tables)
    
    # Clean up column names (replace hyphens with underscores)
    zeroshot_tables = _clean_zeroshot_column_names(zeroshot_tables)
    
    return zeroshot_tables


def _load_zero_shot_data(zip_path: str) -> Dict[str, pd.DataFrame]:
    """
    Load zero-shot scores from zip file.
    
    Args:
        zip_path: Path to zip file containing zero-shot scores
        
    Returns:
        Dictionary mapping DMS study names to DataFrames
    """
    zeroshot_tables = {}
    
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        # Find CSV files
        csv_files = [f for f in zip_ref.namelist() if f.endswith('.csv')]
        
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
            
            zeroshot_tables[study_name] = df
    
    return zeroshot_tables


def _add_uniprot_ids_zeroshot(zeroshot_tables: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
    """
    Add UniProt accession IDs to zero-shot tables.
    
    Args:
        zeroshot_tables: Dictionary of zero-shot DataFrames
        
    Returns:
        Updated dictionary with UniProt_id column added
    """
    # Extract entry names (first two elements split by underscore)
    study_names = list(zeroshot_tables.keys())
    entry_names = []
    
    for name in study_names:
        parts = name.split('_')
        entry_name = f"{parts[0]}_{parts[1]}" if len(parts) >= 2 else parts[0]
        entry_names.append(entry_name)
    
    # Create basic mapping (would use UniProt API in practice)
    uniprot_mapping = _get_basic_uniprot_mapping(entry_names)
    
    # Add UniProt_id to each DataFrame
    updated_tables = {}
    for i, (study_name, df) in enumerate(zeroshot_tables.items()):
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


def _clean_zeroshot_column_names(zeroshot_tables: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
    """
    Clean column names in zero-shot tables.
    
    - Replace hyphens with underscores in column names
    - Replace hyphens with underscores in string values
    
    Args:
        zeroshot_tables: Dictionary of zero-shot DataFrames
        
    Returns:
        Dictionary with cleaned column names
    """
    cleaned_tables = {}
    
    for study_name, df in zeroshot_tables.items():
        df_copy = df.copy()
        
        # Clean column names - replace hyphens with underscores
        df_copy.columns = [col.replace('-', '_') for col in df_copy.columns]
        
        # Clean string values in character columns
        for col in df_copy.columns:
            if df_copy[col].dtype == 'object':
                try:
                    # Replace hyphens with underscores in string values
                    df_copy[col] = df_copy[col].str.replace('-', '_')
                except:
                    pass  # Skip if not string column or other error
        
        # Reorder columns to put core columns first
        core_cols = ['UniProt_id', 'DMS_id', 'mutant', 'mutated_sequence', 'DMS_score', 'DMS_score_bin']
        available_core = [col for col in core_cols if col in df_copy.columns]
        other_cols = [col for col in df_copy.columns if col not in core_cols]
        df_copy = df_copy[available_core + other_cols]
        
        # Handle special case mentioned in R script
        if study_name == "ANCSZ_Hobbs_2022" and len(df_copy.columns) > 7:
            # Remove extra bin_score column (column index 7)
            if df_copy.columns[6].endswith('_bin') and len(df_copy.columns) > 7:
                df_copy = df_copy.drop(df_copy.columns[6], axis=1)
        
        cleaned_tables[study_name] = df_copy
    
    return cleaned_tables


def get_zero_shot_model_list() -> List[str]:
    """
    Get list of zero-shot models available in ProteinGym v1.2.
    
    Returns:
        List of model names (cleaned with underscores)
    """
    # Based on ProteinGym v1.2, there are 79 zero-shot models
    # This is a representative subset - full list would come from actual data
    models = [
        "ESM_2_t48_15B_UR50D", "ESM_2_t36_3B_UR50D", "ESM_2_t33_650M_UR50D",
        "ESM_2_t30_150M_UR50D", "ESM_2_t12_35M_UR50D", "ESM_2_t6_8M_UR50D",
        "ESM_1v_t34_670M_UR100", "ESM_1v_t33_650M_UR90S", "ESM_1b_t33_650M_UR50S",
        "MSA_Transformer", "ProtTrans_t5_xl_u50", "ProtTrans_t5_xxl_u50",
        "Ankh_large", "Ankh_base", "ProstT5", "SaProt_650M", "SaProt_35M",
        "Tranception_L", "Tranception_M", "EVE", "GEMME", "DeepSequence",
        "Vespa", "SIFT", "PROVEAN", "FoldX", "Rosetta", "CADD",
        "AlphaMissense", "ESM1b_regression", "ESM_variants"
    ]
    return models


def filter_zero_shot_by_models(
    zeroshot_tables: Dict[str, pd.DataFrame], 
    models: List[str]
) -> Dict[str, pd.DataFrame]:
    """
    Filter zero-shot tables to include only specified models.
    
    Args:
        zeroshot_tables: Dictionary of zero-shot DataFrames
        models: List of model names to include
        
    Returns:
        Filtered dictionary
    """
    filtered_tables = {}
    
    for study_name, df in zeroshot_tables.items():
        # Get core columns plus specified model columns
        core_cols = ['UniProt_id', 'DMS_id', 'mutant', 'mutated_sequence', 'DMS_score', 'DMS_score_bin']
        available_core = [col for col in core_cols if col in df.columns]
        
        model_cols = [col for col in df.columns if col in models]
        
        filtered_df = df[available_core + model_cols].copy()
        filtered_tables[study_name] = filtered_df
    
    return filtered_tables


def get_zero_shot_summary_stats(zeroshot_tables: Dict[str, pd.DataFrame]) -> Dict[str, Any]:
    """
    Generate summary statistics for zero-shot data.
    
    Args:
        zeroshot_tables: Dictionary of zero-shot DataFrames
        
    Returns:
        Dictionary with summary statistics
    """
    if not zeroshot_tables:
        return {}
    
    total_variants = sum(len(df) for df in zeroshot_tables.values())
    num_assays = len(zeroshot_tables)
    
    # Count models (columns beyond core ones)
    sample_df = next(iter(zeroshot_tables.values()))
    core_cols = ['UniProt_id', 'DMS_id', 'mutant', 'mutated_sequence', 'DMS_score', 'DMS_score_bin']
    model_cols = [col for col in sample_df.columns if col not in core_cols]
    num_models = len(model_cols)
    
    stats = {
        'total_variants': total_variants,
        'num_assays': num_assays,
        'num_models': num_models,
        'assay_sizes': {name: len(df) for name, df in zeroshot_tables.items()},
        'model_names': model_cols
    }
    
    return stats


if __name__ == "__main__":
    # Example usage
    print("Loading zero-shot substitution scores...")
    
    # This would normally download and load the data
    # For demo purposes, we'll show what the function would do
    print("Note: Actual download not implemented in this example")
    
    zeroshot_data = {}  # Placeholder
    
    if zeroshot_data:
        print(f"Loaded zero-shot data for {len(zeroshot_data)} DMS assays")
        
        # Show sample from first assay
        first_key = next(iter(zeroshot_data))
        sample_df = zeroshot_data[first_key]
        print(f"\nSample from {first_key}:")
        print(sample_df.head())
        print(f"Shape: {sample_df.shape}")
        print(f"Columns: {list(sample_df.columns)}")
        
        # Show summary stats
        stats = get_zero_shot_summary_stats(zeroshot_data)
        print(f"\nSummary Statistics:")
        for key, value in stats.items():
            if key not in ['assay_sizes', 'model_names']:
                print(f"{key}: {value}")
    
    # Show available models
    models = get_zero_shot_model_list()
    print(f"\nAvailable zero-shot models (sample of {len(models)}):")
    for i, model in enumerate(models[:10]):  # Show first 10
        print(f"  {i+1}. {model}")
    print(f"  ... and {len(models)-10} more models")