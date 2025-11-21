"""
data_import_funcs.py - Base download utilities for ProteinGym data files.

This module provides centralized cache configuration and download utilities
with caching support, eliminating code duplication across data pipeline modules.
"""

import os
from functools import wraps
from pathlib import Path
from typing import Optional, Union, Callable
import requests


# ============================================================================
# Cache Configuration
# ============================================================================

# Default cache directory - can be overridden via environment variable
DEFAULT_CACHE_DIR = Path(os.getenv("PROTEINGYM_CACHE_DIR", ".cache"))


def get_cache_dir(cache_dir: Optional[Union[str, Path]] = None) -> Path:
    """
    Get the cache directory as a Path object.

    This function provides a centralized way to determine the cache directory,
    with the following precedence:
    1. Explicit cache_dir parameter (if provided)
    2. PROTEINGYM_CACHE_DIR environment variable (if set)
    3. Default ".cache" directory

    Args:
        cache_dir: Optional override for cache directory. Can be a string or Path object.

    Returns:
        Path object for the cache directory

    Examples:
        >>> get_cache_dir()
        PosixPath('.cache')

        >>> get_cache_dir("/tmp/my_cache")
        PosixPath('/tmp/my_cache')

        >>> os.environ["PROTEINGYM_CACHE_DIR"] = "/data/cache"
        >>> get_cache_dir()
        PosixPath('/data/cache')
    """
    if cache_dir is not None:
        return Path(cache_dir)
    return DEFAULT_CACHE_DIR


def set_default_cache_dir(cache_dir: Union[str, Path]) -> None:
    """
    Set the default cache directory globally.

    This updates the DEFAULT_CACHE_DIR module variable and also sets
    the PROTEINGYM_CACHE_DIR environment variable.

    Args:
        cache_dir: New default cache directory

    Example:
        >>> set_default_cache_dir("/data/proteingym_cache")
    """
    global DEFAULT_CACHE_DIR
    DEFAULT_CACHE_DIR = Path(cache_dir)
    os.environ["PROTEINGYM_CACHE_DIR"] = str(DEFAULT_CACHE_DIR)


# ============================================================================
# Download Utilities
# ============================================================================

def cached_download(
    url: str,
    filename: str,
    cache_dir: Optional[Union[str, Path]] = None,
    use_cache: bool = True,
    chunk_size: int = 8192
) -> Path:
    """
    Download a file with caching support.

    If the file already exists in the cache and use_cache is True, the cached
    version is used. Otherwise, the file is downloaded from the URL.

    Args:
        url: URL to download from
        filename: Name for the cached file
        cache_dir: Cache directory (uses default if None)
        use_cache: If True, use cached file if it exists. If False, force a fresh download.
        chunk_size: Download chunk size in bytes (default: 8192)

    Returns:
        Path to the cached file

    Raises:
        requests.HTTPError: If the download fails with an HTTP error
        requests.RequestException: If the download fails for other reasons

    Example:
        >>> zip_path = cached_download(
        ...     url="https://zenodo.org/records/15293562/files/data.zip",
        ...     filename="data.zip",
        ...     cache_dir=".cache"
        ... )
        Using cached file at .cache/data.zip
    """
    cache_path = get_cache_dir(cache_dir)
    cache_path.mkdir(parents=True, exist_ok=True)

    file_path = cache_path / filename

    if not use_cache or not file_path.exists():
        if file_path.exists():
            file_path.unlink()  # Remove existing file

        print(f"Downloading {url} to {file_path}...")
        response = requests.get(url, stream=True)
        response.raise_for_status()

        with open(file_path, "wb") as f:
            for chunk in response.iter_content(chunk_size=chunk_size):
                if chunk:
                    f.write(chunk)
        print("Download complete.")
    else:
        print(f"Using cached file at {file_path}")

    return file_path


def download_with_cache(filename: str, url_param: str = "url"):
    """
    Decorator to add caching support to download functions.

    The decorated function should return a URL string (or take a URL as a parameter).
    This decorator wraps the function to automatically handle downloading and caching.

    Args:
        filename: Name for the cached file
        url_param: Name of the URL parameter in the decorated function (default: "url")

    Returns:
        Decorator function

    Example:
        >>> @download_with_cache("data.zip")
        ... def get_data_url():
        ...     return "https://example.com/data.zip"
        ...
        >>> path = get_data_url(cache_dir=".cache", use_cache=True)
        Downloading https://example.com/data.zip to .cache/data.zip...
    """
    def decorator(func: Callable) -> Callable:
        @wraps(func)
        def wrapper(cache_dir: Optional[Union[str, Path]] = None, use_cache: bool = True, **kwargs):
            # Call the original function to get the URL
            # Check if the function expects a 'url' parameter
            import inspect
            sig = inspect.signature(func)

            if url_param in sig.parameters:
                # Function takes URL as parameter, pass through kwargs
                url = func(**kwargs)
            else:
                # Function returns URL
                url = func(**kwargs)

            cache_path = get_cache_dir(cache_dir)
            return cached_download(url, filename, cache_path, use_cache)
        return wrapper
    return decorator


def download_multiple(
    downloads: list[tuple[str, str]],
    cache_dir: Optional[Union[str, Path]] = None,
    use_cache: bool = True,
    chunk_size: int = 8192
) -> dict[str, Path]:
    """
    Download multiple files with caching support.

    Args:
        downloads: List of (url, filename) tuples to download
        cache_dir: Cache directory (uses default if None)
        use_cache: If True, use cached files if they exist
        chunk_size: Download chunk size in bytes

    Returns:
        Dictionary mapping filenames to their cached paths

    Example:
        >>> files = download_multiple([
        ...     ("https://example.com/data1.zip", "data1.zip"),
        ...     ("https://example.com/data2.zip", "data2.zip"),
        ... ])
        >>> files["data1.zip"]
        PosixPath('.cache/data1.zip')
    """
    results = {}
    for url, filename in downloads:
        results[filename] = cached_download(
            url=url,
            filename=filename,
            cache_dir=cache_dir,
            use_cache=use_cache,
            chunk_size=chunk_size
        )
    return results


# ============================================================================
# Specific ProteinGym Data Downloads
# ============================================================================

def get_dms_substitution_zip(cache_dir: Optional[Union[str, Path]] = None, use_cache: bool = True) -> Path:
    """
    Download the DMS_ProteinGym_substitutions.zip file to the cache directory.

    Args:
        cache_dir: Directory to store the downloaded file (uses default if None)
        use_cache: If True, use cached file if it exists. If False, force a fresh download.

    Returns:
        Path to the downloaded zip file.
    """
    return cached_download(
        url="https://zenodo.org/records/15293562/files/DMS_ProteinGym_substitutions.zip",
        filename="DMS_ProteinGym_substitutions.zip",
        cache_dir=cache_dir,
        use_cache=use_cache
    )


def get_af2_structures_zip(cache_dir: Optional[Union[str, Path]] = None, use_cache: bool = True) -> Path:
    """
    Download the ProteinGym_AF2_structures.zip file to the cache directory.

    Args:
        cache_dir: Directory to store the downloaded file (uses default if None)
        use_cache: If True, use cached file if it exists. If False, force a fresh download.

    Returns:
        Path to the downloaded zip file.
    """
    return cached_download(
        url="https://zenodo.org/records/15293562/files/ProteinGym_AF2_structures.zip?download=1",
        filename="ProteinGym_AF2_structures.zip",
        cache_dir=cache_dir,
        use_cache=use_cache
    )
