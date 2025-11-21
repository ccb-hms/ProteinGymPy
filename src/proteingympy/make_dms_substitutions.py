"""
make_dms_substitutions.py - Python equivalent of make_DMS_substitutions.Rmd

Downloads and processes ProteinGym DMS substitution data.
Loads 217 DMS substitution assays with UniProt ID mapping.
"""

import pandas as pd
import requests
from typing import Dict, List, Optional
from pathlib import Path
import zipfile

from .data_import_funcs import cached_download, get_cache_dir


def get_dms_substitution_data(cache_dir: str = None, use_cache: bool = True) -> Dict[str, pd.DataFrame]:
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
        cache_dir: Directory to cache downloaded files (uses default if None)
        use_cache: If True, use cached file if it exists. If False, force a fresh download.

    Returns:
        Dictionary mapping DMS study names to DataFrames
    """
    # Download using centralized caching utility
    zip_path = cached_download(
        url="https://zenodo.org/records/15293562/files/DMS_ProteinGym_substitutions.zip",
        filename="DMS_ProteinGym_substitutions.zip",
        cache_dir=cache_dir,
        use_cache=use_cache
    )
    
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
    
    # Create mapping - this is a simplified version
    # In a full implementation, you would use the UniProt API or a mapping service
    # For now, we'll create a basic mapping based on the R script's manual curation
    uniprot_mapping = _get_basic_uniprot_mapping(entry_names)
    
    # Add UniProt_id to each DataFrame
    updated_tables = {}
    for i, (study_name, df) in enumerate(progym_tables.items()):
        df_copy = df.copy()
        entry_name = entry_names[i]
        uniprot_id = uniprot_mapping.get(entry_name, None)
        df_copy['UniProt_id'] = uniprot_id
        updated_tables[study_name] = df_copy
    
    return updated_tables


def _get_basic_uniprot_mapping(entry_names: List[str]) -> Dict[str, Optional[str]]:
    """
    Map UniProt entry names to accession IDs using UniProt REST API.
    
    Args:
        entry_names: List of UniProt entry names (e.g., 'P53_HUMAN')
        
    Returns:
        Dictionary mapping entry name to UniProt accession ID
    """
    mapping = {}
    
    # Filter out special cases and duplicates
    unique_names = list(set(entry_names))
    names_to_query = []
    
    for name in unique_names:
        if name == "ANCSZ_Hobbs":
            mapping[name] = None
        else:
            names_to_query.append(name)
            
    if not names_to_query:
        return mapping
        
    # Batch queries to avoid URL length limits
    batch_size = 50
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    
    print(f"Querying UniProt API for {len(names_to_query)} entries...")
    
    for i in range(0, len(names_to_query), batch_size):
        batch = names_to_query[i:i+batch_size]
        
        # Construct query: id:NAME1 OR id:NAME2 ...
        query_parts = [f"id:{name}" for name in batch]
        query = " OR ".join(query_parts)
        
        params = {
            "query": query,
            "fields": "accession,id",
            "format": "json",
            "size": len(batch)
        }
        
        try:
            response = requests.get(base_url, params=params)
            response.raise_for_status()
            
            results = response.json().get("results", [])
            
            for result in results:
                # API returns 'primaryAccession' and 'uniProtkbId' (entry name)
                accession = result.get("primaryAccession")
                entry_name = result.get("uniProtkbId")
                
                if entry_name and accession:
                    mapping[entry_name] = accession
                    
        except Exception as e:
            print(f"Error querying UniProt API for batch {i//batch_size + 1}: {e}")
            
    # Ensure all requested names are in the mapping (None if not found)
    for name in entry_names:
        if name not in mapping:
            mapping[name] = None
            
    return mapping


def get_dms_metadata(cache_dir: str = None) -> pd.DataFrame:
    """
    Download and process DMS substitutions metadata/reference file.

    Args:
        cache_dir: Directory to cache downloaded files (uses default if None)

    Returns:
        DataFrame with metadata for 217 DMS assays
    """
    # Download using centralized caching utility
    metadata_path = cached_download(
        url="https://zenodo.org/records/15293562/files/DMS_substitutions.csv",
        filename="DMS_substitutions.csv",
        cache_dir=cache_dir,
        use_cache=True  # Metadata rarely changes, always use cache
    )

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