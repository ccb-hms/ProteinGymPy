import os
import requests

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