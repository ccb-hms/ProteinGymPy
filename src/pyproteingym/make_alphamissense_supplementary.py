"""
make_alphamissense_supplementary.py - Python equivalent of make_AM_supplementary.Rmd

Downloads and processes AlphaMissense pathogenicity scores for ProteinGym variants.
"""

import os
import pandas as pd
import requests
import zipfile
import tempfile
from typing import Optional


def get_alphamissense_proteingym_data(cache_dir: str = ".cache") -> pd.DataFrame:
    """
    Download and process AlphaMissense supplementary data for ProteinGym variants.
    
    This loads Table S8 from Cheng et al. 2023 containing AlphaMissense pathogenicity 
    scores for ~1.6M variants that match those in ProteinGym from 87 DMS experiments 
    across 72 proteins.
    
    Args:
        cache_dir: Directory to cache downloaded files
        
    Returns:
        DataFrame with columns:
        - DMS_id: DMS assay identifier  
        - Uniprot_ID: UniProt accession
        - variant_id: Variant identifier
        - AlphaMissense: Pathogenicity score (0-1, higher = more pathogenic)
    """
    os.makedirs(cache_dir, exist_ok=True)
    
    # File paths
    zip_path = os.path.join(cache_dir, "science.adg7492_data_s1_to_s9.zip")
    csv_path = os.path.join(cache_dir, "Supplementary_Data_S8_proteingym.csv")
    
    # Download supplementary data zip if not cached
    if not os.path.exists(csv_path):
        if not os.path.exists(zip_path):
            url = "https://www.science.org/doi/suppl/10.1126/science.adg7492/suppl_file/science.adg7492_data_s1_to_s9.zip"
            print(f"Downloading AlphaMissense supplementary data from {url}...")
            response = requests.get(url, stream=True)
            response.raise_for_status()
            with open(zip_path, "wb") as f:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)
            print("Download complete.")
        
        # Extract the specific file we need
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            # Look for the ProteinGym supplementary file
            target_file = None
            for file in zip_ref.namelist():
                if "proteingym" in file.lower() and file.endswith('.csv'):
                    target_file = file
                    break
            
            if target_file:
                print(f"Extracting {target_file}...")
                zip_ref.extract(target_file, cache_dir)
                # Move to standard name if needed
                extracted_path = os.path.join(cache_dir, target_file)
                if extracted_path != csv_path:
                    os.rename(extracted_path, csv_path)
            else:
                raise FileNotFoundError("Could not find ProteinGym supplementary CSV in zip file")
    
    # Load the data
    print("Loading AlphaMissense ProteinGym data...")
    df = pd.read_csv(csv_path)
    
    # Ensure AlphaMissense column is numeric
    df['AlphaMissense'] = pd.to_numeric(df['AlphaMissense'], errors='coerce')
    
    print(f"Loaded {len(df):,} AlphaMissense scores for ProteinGym variants")
    print(f"Data covers {df['DMS_id'].nunique()} DMS assays")
    print(f"Columns: {list(df.columns)}")
    
    return df


def get_alphamissense_summary_stats(df: pd.DataFrame) -> dict:
    """
    Generate summary statistics for AlphaMissense ProteinGym data.
    
    Args:
        df: AlphaMissense DataFrame from get_alphamissense_proteingym_data()
        
    Returns:
        Dictionary with summary statistics
    """
    stats = {
        'total_variants': len(df),
        'unique_dms_assays': df['DMS_id'].nunique(),
        'unique_proteins': df['Uniprot_ID'].nunique(),
        'alphamissense_mean': df['AlphaMissense'].mean(),
        'alphamissense_median': df['AlphaMissense'].median(),
        'alphamissense_std': df['AlphaMissense'].std(),
        'alphamissense_min': df['AlphaMissense'].min(),
        'alphamissense_max': df['AlphaMissense'].max(),
    }
    
    return stats


def filter_alphamissense_by_dms(df: pd.DataFrame, dms_ids: list) -> pd.DataFrame:
    """
    Filter AlphaMissense data for specific DMS assays.
    
    Args:
        df: AlphaMissense DataFrame
        dms_ids: List of DMS assay identifiers to include
        
    Returns:
        Filtered DataFrame
    """
    return df[df['DMS_id'].isin(dms_ids)].copy()


if __name__ == "__main__":
    # Example usage
    print("Loading AlphaMissense ProteinGym supplementary data...")
    am_data = get_alphamissense_proteingym_data()
    
    # Show summary statistics
    stats = get_alphamissense_summary_stats(am_data)
    print(f"\nSummary Statistics:")
    for key, value in stats.items():
        if isinstance(value, float):
            print(f"{key}: {value:.4f}")
        else:
            print(f"{key}: {value:,}")
    
    # Show sample data
    print(f"\nSample data:")
    print(am_data.head())
    
    # Show distribution across assays
    print(f"\nTop 10 assays by variant count:")
    assay_counts = am_data['DMS_id'].value_counts().head(10)
    print(assay_counts)