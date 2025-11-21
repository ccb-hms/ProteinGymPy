import os
import requests
from typing import Dict, List, Optional

def get_dms_substitution_zip(cache_dir: str = ".cache/", use_cache: bool = True) -> str:
    """Download the DMS_ProteinGym_substitutions.zip file to the cache directory.
    
    Args:
        cache_dir: Directory to store the downloaded file.
        use_cache: If True, use cached file if it exists. If False, force a fresh download.
    
    Returns:
        Path to the downloaded zip file.
    """
    url = "https://zenodo.org/records/15293562/files/DMS_ProteinGym_substitutions.zip"
    os.makedirs(cache_dir, exist_ok=True)
    zip_path = os.path.join(cache_dir, "DMS_ProteinGym_substitutions.zip")
    
    if not use_cache or not os.path.exists(zip_path):
        if os.path.exists(zip_path):
            os.remove(zip_path)  
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
    return zip_path


def get_af2_structures_zip(cache_dir: str = ".cache/", use_cache: bool = True) -> str:
    """Download the ProteinGym_AF2_structures.zip file to the cache directory.

    Args:
        cache_dir: Directory to store the downloaded file.
        use_cache: If True, use cached file if it exists. If False, force a fresh download.

    Returns:
        Path to the downloaded zip file.
    """
    url = (
        "https://zenodo.org/records/15293562/files/ProteinGym_AF2_structures.zip?download=1"
    )
    os.makedirs(cache_dir, exist_ok=True)
    zip_path = os.path.join(cache_dir, "ProteinGym_AF2_structures.zip")

    if not use_cache or not os.path.exists(zip_path):
        if os.path.exists(zip_path):
            os.remove(zip_path)
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

    return zip_path


def _query_uniprot_api(entry_names: List[str]) -> Dict[str, Optional[str]]:
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