import sys
import os
from src.pyproteingym import get_dms_substitution_zip

def main():
    cache_dir = ".cache"
    zip_path = os.path.join(cache_dir, "DMS_ProteinGym_substitutions.zip")
    if not os.path.exists(zip_path):
        get_dms_substitution_zip(cache_dir)
    else:
        print(f"File already exists at {zip_path}.")


# if main.py is run in namespace -- then run these functions.
if __name__ == "__main__":
    main()

