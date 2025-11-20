"""
make_dms_substitutions.py - Python equivalent of make_DMS_substitutions.Rmd

Downloads and processes ProteinGym DMS substitution data.
Loads 217 DMS substitution assays with UniProt ID mapping.
"""

import os
import pandas as pd
import requests
from typing import Dict, List, Optional
import tempfile
import zipfile
from .data_import_funcs import _query_uniprot_api


def get_dms_substitution_data(cache_dir: str = ".cache", use_cache: bool = True) -> Dict[str, pd.DataFrame]:
    """
    Download and process ProteinGym DMS substitution data.
    
    Returns a dictionary of 217 DMS assays, each as a pandas DataFrame with columns:
    - UniProt_id: UniProt accession identifier  
    - DMS_id: DMS assay identifier
    - mutant: substitution description (e.g. A1P:D2N)
    - mutated_sequence: full amino acid sequence
    - DMS_score: experimental measurement (higher = more fit)
    - DMS_score_bin: binary fitness (1=fit, 0=not fit)
    
    Args:
        cache_dir: Directory to cache downloaded files
        use_cache: If True, use cached file if it exists. If False, force a fresh download.
        
    Returns:
        Dictionary mapping DMS study names to DataFrames
    """
    os.makedirs(cache_dir, exist_ok=True)
    zip_path = os.path.join(cache_dir, "DMS_ProteinGym_substitutions.zip")
    
    # Download if not cached or if use_cache is False
    if not use_cache or not os.path.exists(zip_path):
        if os.path.exists(zip_path):
            os.remove(zip_path)
        url = "https://zenodo.org/records/15293562/files/DMS_ProteinGym_substitutions.zip"
        print(f"Downloading {url} to {zip_path}...")
        response = requests.get(url, stream=True)
        response.raise_for_status()
        with open(zip_path, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
        print("Download complete.")
    else:
        print(f"Using cached file at {zip_path}.")
    
    # Extract and load data
    progym_tables = {}
    
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        file_list = [f for f in zip_ref.namelist() if f.endswith('.csv')]
        
        for csv_file in file_list:
            # Extract DMS study name (remove .csv extension)
            study_name = os.path.splitext(os.path.basename(csv_file))[0]
            
            # Read CSV from zip
            with zip_ref.open(csv_file) as f:
                df = pd.read_csv(f)
            
            # Add DMS_id column
            df['DMS_id'] = study_name
            
            # Convert DMS_score_bin to categorical
            df['DMS_score_bin'] = df['DMS_score_bin'].astype('category')
            
            progym_tables[study_name] = df
    
    # Add UniProt IDs
    progym_tables = _add_uniprot_ids(progym_tables)
    
    # Reorder columns
    cols = ['UniProt_id', 'DMS_id', 'mutant', 'mutated_sequence', 'DMS_score', 'DMS_score_bin']
    
    for study_name, df in progym_tables.items():
        # Select and reorder columns (keep any additional columns at end)
        available_cols = [col for col in cols if col in df.columns]
        other_cols = [col for col in df.columns if col not in cols]
        progym_tables[study_name] = df[available_cols + other_cols]
    
    return progym_tables


def _add_uniprot_ids(progym_tables: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
    """
    Add UniProt accession IDs to ProteinGym tables.
    
    Args:
        progym_tables: Dictionary of DMS DataFrames
        
    Returns:
        Updated dictionary with UniProt_id column added
    """
    # Extract entry names (first two elements split by underscore)
    study_names = list(progym_tables.keys())
    entry_names = []
    
    for name in study_names:
        parts = name.split('_')
        entry_name = f"{parts[0]}_{parts[1]}" if len(parts) >= 2 else parts[0]
        entry_names.append(entry_name)
    
    # Create mapping using UniProt API
    uniprot_mapping = _query_uniprot_api(entry_names)
    
    # Add UniProt_id to each DataFrame
    updated_tables = {}
    for i, (study_name, df) in enumerate(progym_tables.items()):
        df_copy = df.copy()
        entry_name = entry_names[i]
        uniprot_id = uniprot_mapping.get(entry_name, None)
        df_copy['UniProt_id'] = uniprot_id
        updated_tables[study_name] = df_copy
    
    return updated_tables


def get_dms_metadata(cache_dir: str = ".cache") -> pd.DataFrame:
    """
    Download and process DMS substitutions metadata/reference file.
    
    Args:
        cache_dir: Directory to cache downloaded files
        
    Returns:
        DataFrame with metadata for 217 DMS assays
    """
    os.makedirs(cache_dir, exist_ok=True)
    metadata_path = os.path.join(cache_dir, "DMS_substitutions.csv")
    
    if not os.path.exists(metadata_path):
        url = "https://zenodo.org/records/15293562/files/DMS_substitutions.csv"
        print(f"Downloading metadata from {url}...")
        response = requests.get(url, stream=True)
        response.raise_for_status()
        with open(metadata_path, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
        print("Metadata download complete.")
    
    # Load and process metadata
    df = pd.read_csv(metadata_path)
    
    # Convert categorical columns
    categorical_cols = [
        'taxon', 'source_organism', 'DMS_binarization_method', 
        'selection_type', 'selection_assay', 'raw_DMS_phenotype_name',
        'raw_DMS_directionality', 'raw_DMS_mutant_column', 
        'ProteinGym_version', 'coarse_selection_type'
    ]
    
    for col in categorical_cols:
        if col in df.columns:
            df[col] = df[col].astype('category')
    
    return df


if __name__ == "__main__":
    # Example usage
    print("Loading DMS substitution data...")
    dms_data = get_dms_substitution_data()
    print(f"Loaded {len(dms_data)} DMS assays")
    
    # Show sample from first assay
    first_key = next(iter(dms_data))
    sample_df = dms_data[first_key]
    print(f"\nSample from {first_key}:")
    print(sample_df.head())
    print(f"Shape: {sample_df.shape}")
    print(f"Columns: {list(sample_df.columns)}")