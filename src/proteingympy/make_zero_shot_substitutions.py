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
from .uniprot_utils import _query_uniprot_api


def get_zero_shot_substitution_data(cache_dir: str = ".cache") -> Dict[str, pd.DataFrame]:
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
    else:
        print(f"Zero-shot scores found in cache at {zip_path}")
    
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
    
    # Create basic mapping using UniProt API
    uniprot_mapping = _query_uniprot_api(entry_names)
    
    # Add UniProt_id to each DataFrame
    updated_tables = {}
    for i, (study_name, df) in enumerate(zeroshot_tables.items()):
        df_copy = df.copy()
        entry_name = entry_names[i]
        uniprot_id = uniprot_mapping.get(entry_name, None)
        df_copy['UniProt_id'] = uniprot_id
        updated_tables[study_name] = df_copy
    
    return updated_tables





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


def available_zero_shot_models() -> List[str]:
    """
    Get list of zero-shot models available in ProteinGym v1.2.
    
    Returns:
        List of model names (cleaned with underscores)
    """
    # Based on ProteinGym v1.2, there are 79 zero-shot models
    models = [
        "Site_Independent", "EVmutation",
        "DeepSequence_single", "DeepSequence_ensemble",
        "EVE_single", "EVE_ensemble",
        "Unirep", "Unirep_evotune",
        "MSA_Transformer_single", "MSA_Transformer_ensemble",
        "ESM1b", "ESM1v_single",
        "ESM1v_ensemble", "ESM2_8M",
        "ESM2_35M", "ESM2_150M",
        "ESM2_650M", "ESM2_3B",
        "ESM2_15B", "Wavenet",
        "RITA_s", "RITA_m",
        "RITA_l", "RITA_xl",
        "Progen2_small", "Progen2_medium",
        "Progen2_base", "Progen2_large",
        "Progen2_xlarge", "GEMME",
        "VESPA", "VESPAl",
        "VespaG", "ProtGPT2",
        "Tranception_S_no_retrieval", "Tranception_M_no_retrieval",
        "Tranception_L_no_retrieval", "Tranception_S",
        "Tranception_M", "Tranception_L",
        "TranceptEVE_S", "TranceptEVE_M",
        "TranceptEVE_L", "CARP_38M",
        "CARP_600K", "CARP_640M",
        "CARP_76M", "MIF",
        "MIFST", "ESM_IF1",
        "ProteinMPNN", "ProtSSN_k10_h512",
        "ProtSSN_k10_h768", "ProtSSN_k10_h1280",
        "ProtSSN_k20_h512", "ProtSSN_k20_h768",
        "ProtSSN_k20_h1280", "ProtSSN_k30_h512",
        "ProtSSN_k30_h768", "ProtSSN_k30_h1280",
        "ProtSSN_ensemble", "SaProt_650M_AF2",
        "SaProt_35M_AF2", "PoET",
        "MULAN_small", "ProSST_20",
        "ProSST_128", "ProSST_512",
        "ProSST_1024", "ProSST_2048",
        "ProSST_4096", "ESCOTT",
        "VenusREM", "RSALOR",
        "S2F", "S2F_MSA",
        "S3F", "S3F_MSA",
        "SiteRM"
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

# some random stats about the DataFrame chatbot generated -- dont need.
# Made function hidden
def _get_zero_shot_summary_stats(zeroshot_tables: Dict[str, pd.DataFrame]) -> Dict[str, Any]:
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
        stats = _get_zero_shot_summary_stats(zeroshot_data)
        print(f"\nSummary Statistics:")
        for key, value in stats.items():
            if key not in ['assay_sizes', 'model_names']:
                print(f"{key}: {value}")
    
    # Show available models
    models = available_zero_shot_models()
    print(f"\nAvailable zero-shot models (sample of {len(models)}):")
    for i, model in enumerate(models[:10]):  # Show first 10
        print(f"  {i+1}. {model}")
    print(f"  ... and {len(models)-10} more models")