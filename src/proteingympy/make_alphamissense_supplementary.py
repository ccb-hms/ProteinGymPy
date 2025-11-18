"""
make_alphamissense_supplementary.py - Python equivalent of make_AM_supplementary.Rmd

Downloads and processes AlphaMissense pathogenicity scores for ProteinGym variants.
"""

import io
import json
import os
import re
import zipfile
from typing import Dict, Optional

import pandas as pd
import requests


def _is_proteingym_csv(member_name: str) -> bool:
    """Return True if the archive member looks like the ProteinGym CSV."""
    lower = member_name.lower()
    return lower.endswith('.csv') and 'proteingym' in lower and '__macosx' not in lower


def _extract_proteingym_csv(zip_ref: zipfile.ZipFile, dest_path: str) -> Optional[str]:
    """Extract the ProteinGym CSV from a (possibly nested) zip archive."""
    # First, check top-level members
    for member in zip_ref.namelist():
        if _is_proteingym_csv(member):
            os.makedirs(os.path.dirname(dest_path), exist_ok=True)
            with zip_ref.open(member) as src, open(dest_path, 'wb') as dst:
                dst.write(src.read())
            return member

    # Otherwise, recurse into any nested zip members
    for member in zip_ref.namelist():
        if member.lower().endswith('.zip'):
            with zip_ref.open(member) as nested_file:
                nested_bytes = nested_file.read()
            with zipfile.ZipFile(io.BytesIO(nested_bytes), 'r') as nested_zip:
                nested_result = _extract_proteingym_csv(nested_zip, dest_path)
                if nested_result:
                    return f"{member}->{nested_result}"

    return None


def _load_cached_uniprot_mapping(cache_path: str) -> Dict[str, Optional[str]]:
    if os.path.exists(cache_path):
        try:
            with open(cache_path, 'r', encoding='utf-8') as handle:
                data = json.load(handle)
            return {str(k): (v if v not in (None, "None") else None) for k, v in data.items()}
        except (json.JSONDecodeError, OSError):
            pass
    return {}


def _save_uniprot_mapping(cache_path: str, mapping: Dict[str, Optional[str]]) -> None:
    tmp_path = f"{cache_path}.tmp"
    with open(tmp_path, 'w', encoding='utf-8') as handle:
        json.dump(mapping, handle, indent=2, sort_keys=True)
    os.replace(tmp_path, cache_path)


def _query_uniprot_accessions(entry_names, session: Optional[requests.Session] = None, chunk_size: int = 50) -> Dict[str, Optional[str]]:
    """Query UniProt REST API to map UniProt entry names (SwissProt IDs) to accessions."""
    if not entry_names:
        return {}

    session = session or requests.Session()
    base_url = "https://rest.uniprot.org/uniprotkb/stream"
    params = {
        "format": "tsv",
        "fields": "accession,id",
    }

    results: Dict[str, Optional[str]] = {entry: None for entry in entry_names}

    for i in range(0, len(entry_names), chunk_size):
        chunk = entry_names[i : i + chunk_size]
        query = " OR ".join(f"id:{entry}" for entry in chunk)
        response = session.get(base_url, params={**params, "query": query}, timeout=30)
        response.raise_for_status()
        text = response.text.strip()
        if not text:
            continue
        lines = text.splitlines()
        # Expected header: Entry\tEntry Name
        for line in lines[1:]:
            parts = line.split('\t')
            if len(parts) < 2:
                continue
            accession, entry_name = parts[0].strip(), parts[1].strip()
            if entry_name:
                results[entry_name] = accession

    return results


def _add_uniprot_accessions(df: pd.DataFrame, cache_dir: str) -> pd.DataFrame:
    """Augment AlphaMissense data with UniProt accessions via the UniProt REST API."""
    entry_col = 'Uniprot_ID'
    if entry_col not in df.columns:
        return df

    cache_path = os.path.join(cache_dir, "alphamissense_uniprot_mapping.json")
    cached_mapping = _load_cached_uniprot_mapping(cache_path)

    entry_names = sorted(df[entry_col].dropna().unique())
    missing = [entry for entry in entry_names if cached_mapping.get(entry) in (None, '')]

    if missing:
        try:
            fetched = _query_uniprot_accessions(missing)
        except requests.RequestException as exc:
            print(
                "Warning: Unable to contact UniProt API to map SwissProt IDs to accessions. "
                f"Reason: {exc}. Falling back to original IDs."
            )
            df['SwissProt_ID'] = df[entry_col]
            return df

        # Merge fetched results into cache
        cached_mapping.update({entry: fetched.get(entry, cached_mapping.get(entry)) for entry in missing})
        _save_uniprot_mapping(cache_path, cached_mapping)

    df = df.copy()
    df['SwissProt_ID'] = df[entry_col]
    df[entry_col] = df['SwissProt_ID'].map(cached_mapping)

    missing_mask = df[entry_col].isna()
    if missing_mask.any():
        # Attempt to derive accession directly from SwissProt ID (prefix before underscore)
        derived_series = df.loc[missing_mask, 'SwissProt_ID'].apply(_derive_accession_from_entry)
        df.loc[missing_mask, entry_col] = derived_series

        # Cache any newly derived accessions
        derived_mapping = {
            swissprot: accession
            for swissprot, accession in zip(df.loc[missing_mask, 'SwissProt_ID'], derived_series)
            if accession
        }
        if derived_mapping:
            cached_mapping.update(derived_mapping)
            _save_uniprot_mapping(cache_path, cached_mapping)

        missing_mask = df[entry_col].isna()

    if missing_mask.any():
        unresolved = df.loc[missing_mask, 'SwissProt_ID'].unique()
        print(
            "Warning: UniProt accessions could not be resolved for the following SwissProt IDs: "
            + ", ".join(unresolved)
        )
        df[entry_col] = df[entry_col].fillna(df['SwissProt_ID'])

    return df


_ACCESSION_REGEX = re.compile(r'^[A-Z0-9]{6,10}$')


def _derive_accession_from_entry(entry_name: str) -> Optional[str]:
    """Best-effort fallback to derive an accession from a SwissProt entry name."""
    if not entry_name or '_' not in entry_name:
        return None
    candidate = entry_name.split('_', 1)[0]
    if _ACCESSION_REGEX.match(candidate):
        return candidate
    return None


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
        - Uniprot_ID: UniProt accession (resolved via UniProt API)
        - SwissProt_ID: Original AlphaMissense SwissProt entry name
    - variant_id: Variant identifier
    - AlphaMissense: Pathogenicity score (0-1, higher = more pathogenic)
    """
    os.makedirs(cache_dir, exist_ok=True)
    
    # File paths
    csv_path = os.path.join(cache_dir, "Supplementary_Data_S8_proteingym.csv")

    #url = "https://www.science.org/doi/suppl/10.1126/science.adg7492/suppl_file/science.adg7492_data_s1_to_s9.zip"
    # Science is blocking requests with TLS fingerprinting, so we rely on a local copy
    # Preferred zip path is the copy bundled with the package at src/
    repo_zip_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "science.adg7492_data_s1_to_s9.zip"))
    cache_zip_path = os.path.join(cache_dir, "science.adg7492_data_s1_to_s9.zip")

    # Prefer the repository zip file if present, otherwise fallback to cache
    if os.path.exists(repo_zip_path):
        zip_path = repo_zip_path
    else:
        zip_path = cache_zip_path

    # Extract CSV if not present
    if not os.path.exists(csv_path):
        if os.path.exists(zip_path):
            with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                extracted_name = _extract_proteingym_csv(zip_ref, csv_path)
                if extracted_name:
                    print(f"Extracted {extracted_name} from {zip_path} -> {csv_path}")
                else:
                    raise FileNotFoundError("Could not find ProteinGym supplementary CSV in zip file (including nested archives)")
        else:
            # Neither the CSV nor any local zip exists; we do not attempt to download
            raise FileNotFoundError(
                f"AlphaMissense supplementary CSV not found at {csv_path} and no local zip found at {repo_zip_path} or {cache_zip_path}."
            )
    
    # Load the data
    print("Loading AlphaMissense ProteinGym data...")
    
    df = pd.read_csv(csv_path)
    df = _add_uniprot_accessions(df, cache_dir)
    
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