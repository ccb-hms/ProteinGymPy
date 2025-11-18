"""
Plot protein structures with DMS and model scores.

This module provides functionality to visualize DMS (Deep Mutational Scanning)
or model scores for amino acid substitutions on 3D protein structures using nglview.
"""

import warnings
import zipfile
from typing import Optional, Callable, List, Dict, Any, Union, Tuple, cast
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats
import nglview as nv
from Bio.PDB.PDBParser import PDBParser

from .data_import_funcs import get_af2_structures_zip

DEFAULT_CACHE_DIR = Path('.cache')


def filter_by_pos(
    df: pd.DataFrame,
    start_pos: Optional[int] = None,
    end_pos: Optional[int] = None
) -> pd.DataFrame:
    """
    Filter dataframe by position range.
    
    Parameters
    ----------
    df : pd.DataFrame
        Dataframe containing a 'pos' column with integer positions
    start_pos : int, optional
        First amino acid position to include. If None, uses minimum position.
    end_pos : int, optional
        Last amino acid position to include. If None, uses maximum position.
        
    Returns
    -------
    pd.DataFrame
        Filtered dataframe containing only rows within the specified position range
        
    Raises
    ------
    ValueError
        If 'pos' column is missing or not integer type, or if positions are out of range
    """
    # Check pos column exists
    if 'pos' not in df.columns:
        raise ValueError("The dataframe must contain a 'pos' column.")
    
    # Check pos column is integer
    if not pd.api.types.is_integer_dtype(df['pos']):
        raise ValueError("The 'pos' column must be an integer type.")
    
    # Get min/max positions
    min_pos = df['pos'].min()
    max_pos = df['pos'].max()
    
    # Validate user-provided positions
    if start_pos is not None and start_pos > max_pos:
        raise ValueError(
            f"start_pos ({start_pos}) is outside the assay range "
            f"({min_pos} to {max_pos})"
        )
    if end_pos is not None and end_pos < min_pos:
        raise ValueError(
            f"end_pos ({end_pos}) is outside the assay range "
            f"({min_pos} to {max_pos})"
        )
    
    # Set defaults if not provided
    if start_pos is None:
        start_pos = min_pos
    if end_pos is None:
        end_pos = max_pos
    
    # Filter the dataframe
    filtered_df = df[(df['pos'] >= start_pos) & (df['pos'] <= end_pos)].copy()
    
    return filtered_df


def get_prot_ids(names: Union[str, List[str]]) -> Union[str, List[str]]:
    """
    Extract protein ID from assay names.
    
    Extracts the first two underscore-separated parts of assay names.
    For example, "C6KNH7_9INFA_Lee_2018" returns "C6KNH7_9INFA".
    
    Parameters
    ----------
    names : str or list of str
        Assay name(s) to extract protein IDs from
        
    Returns
    -------
    str or list of str
        Protein ID(s) extracted from assay names
    """
    def extract_id(name: str) -> str:
        parts = name.split('_')
        return '_'.join(parts[:2]) if len(parts) >= 2 else name
    
    if isinstance(names, str):
        return extract_id(names)
    else:
        return [extract_id(name) for name in names]


def quantile_normalize_scores(scores: Union[np.ndarray, Any]) -> np.ndarray:
    """
    Normalize scores using rank-based normal quantile transformation.
    
    This transformation converts values into z-scores using a rank-based 
    normal quantile approach:
    1. Compute empirical CDF to get percentile ranks
    2. Apply inverse normal CDF (quantile function) to convert to z-scores
    
    The resulting values preserve rank order but are approximately normally
    distributed (mean=0, SD=1). Typical values range between -3 and 3.
    
    Parameters
    ----------
    scores : np.ndarray or array-like
        Array of scores to normalize
        
    Returns
    -------
    np.ndarray
        Quantile-normalized scores (z-scores)
        
    Raises
    ------
    ImportError
        If scipy is not installed
    """
    if stats is None:
        raise ImportError("scipy is required for quantile normalization")
    
    # Convert to numpy array if needed
    scores_array = np.asarray(scores)
    
    # Calculate ranks (handling ties with average method)
    ranks = stats.rankdata(scores_array, method='average')
    
    # Convert ranks to probabilities [0, 1]
    n = len(scores_array)
    probabilities = ranks / (n + 1)  # Use (n+1) to avoid 0 and 1
    
    # Apply inverse normal CDF
    quantile_scores = stats.norm.ppf(probabilities)
    
    return quantile_scores


def color_line(
    df: pd.DataFrame,
    quant_norm: bool = True,
    color_palette: Optional[List[str]] = None,
    n: int = 200
) -> pd.DataFrame:
    """
    Map aggregate scores to colors using quantile normalization.
    
    Parameters
    ----------
    df : pd.DataFrame
        Dataframe with 'aggregate_score' column
    quant_norm : bool, default=True
        Whether to apply quantile normalization
    color_palette : list of str, optional
        List of hex color codes. If None, uses a blue-cyan-yellow palette.
    n : int, default=200
        Number of colors in the palette
        
    Returns
    -------
    pd.DataFrame
        Input dataframe with added 'quant_score' and 'color' columns
    """
    df = df.copy()
    
    # Create default color palette if not provided (parula-like)
    if color_palette is None:
        # Simple approximation of MATLAB's parula colormap
        color_palette = _create_parula_palette(n)
    
    if quant_norm:
        # Apply quantile normalization
        df['quant_score'] = quantile_normalize_scores(df['aggregate_score'].values)
        
        # Clamp scores between -3 and 3
        df['quant_clamped'] = np.clip(df['quant_score'], -3, 3)
        
        # Map -3 to 0 and 3 to (n-1)
        color_indices = ((df['quant_clamped'] + 3) / 6 * (n - 1)).astype(int)
        
        # Assign colors
        df['color'] = [color_palette[idx] for idx in color_indices]
    
    return df


def _create_parula_palette(n: int = 200) -> List[str]:
    """
    Create a parula-like color palette.
    
    Approximates MATLAB's parula colormap with blue-cyan-yellow gradient.
    
    Parameters
    ----------
    n : int, default=200
        Number of colors to generate
        
    Returns
    -------
    list of str
        List of hex color codes
    """
    try:
        from matplotlib import cm
        from matplotlib.colors import rgb2hex
        
        # Use viridis as approximation (similar to parula)
        cmap = cm.get_cmap('viridis', n)
        colors = [rgb2hex(cmap(i)) for i in range(n)]
        return colors
    except ImportError:
        warnings.warn("matplotlib not available, using simple gradient")
        # Fallback: simple blue to yellow gradient
        colors = []
        for i in range(n):
            ratio = i / (n - 1)
            r = int(ratio * 255)
            g = int(ratio * 255)
            b = int((1 - ratio) * 255)
            colors.append(f'#{r:02x}{g:02x}{b:02x}')
        return colors


def get_color_function(
    color_scheme: Optional[str],
    values: List[float]
) -> Callable[[float], str]:
    """
    Create a color mapping function for score values.
    
    Parameters
    ----------
    color_scheme : str, optional
        Color scheme name. Options: 'EVE' for EVE-style coloring,
        or None for red-white-blue scheme.
    values : list of float
        Breakpoint values for color mapping (e.g., [min, 0, max])
        
    Returns
    -------
    callable
        Function that maps a numeric value to a hex color code
    """
    if color_scheme == "EVE":
        # EVE color scheme: black -> purple -> cyan -> yellow
        colors = ['#000000', '#9440e8', '#00CED1', '#fde662']
    else:
        # Default: red -> white -> blue
        colors = ['#ff0000', '#ffffff', '#0000ff']
    
    # Create interpolation function
    def interpolate_color(value: float) -> str:
        """Map a value to a color using linear interpolation."""
        # Clamp value to range
        value = np.clip(value, values[0], values[-1])
        
        # Find which segment the value falls in
        for i in range(len(values) - 1):
            if values[i] <= value <= values[i + 1]:
                # Linear interpolation between colors
                t = (value - values[i]) / (values[i + 1] - values[i])
                
                # Interpolate RGB components
                color1 = _hex_to_rgb(colors[i])
                color2 = _hex_to_rgb(colors[i + 1])
                
                r = int(color1[0] + t * (color2[0] - color1[0]))
                g = int(color1[1] + t * (color2[1] - color1[1]))
                b = int(color1[2] + t * (color2[2] - color1[2]))
                
                return f'#{r:02x}{g:02x}{b:02x}'
        
        return colors[-1]
    
    return interpolate_color


def _hex_to_rgb(hex_color: str) -> Tuple[int, int, int]:
    """Convert hex color to RGB tuple."""
    hex_color = hex_color.lstrip('#')
    r = int(hex_color[0:2], 16)
    g = int(hex_color[2:4], 16)
    b = int(hex_color[4:6], 16)
    return (r, g, b)


def _rgb_to_hex(rgb: Tuple[int, int, int]) -> str:
    """Convert RGB tuple to hex color."""
    return f'#{rgb[0]:02x}{rgb[1]:02x}{rgb[2]:02x}'


def _ensure_af2_structures(cache_dir: Union[str, Path] = DEFAULT_CACHE_DIR) -> Path:
    """Ensure AF2 structure archive is downloaded and extracted.

    Parameters
    ----------
    cache_dir : str or Path, optional
        Directory where the archive and extracted files should live.

    Returns
    -------
    Path
        Path to the directory containing extracted PDB files.
    """
    cache_dir = Path(cache_dir)
    structures_root = cache_dir / 'ProteinGym_AF2_structures'
    cache_dir.mkdir(parents=True, exist_ok=True)

    zip_path = Path(
        get_af2_structures_zip(cache_dir=str(cache_dir), use_cache=True)
    )

    needs_extract = True
    if structures_root.exists():
        # Check if there are any PDB files already extracted
        needs_extract = not any(structures_root.rglob('*.pdb'))

    if needs_extract:
        structures_root.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall(structures_root)

    return structures_root


def _find_pdb_for_prot_id(prot_id: str, structures_root: Path) -> Path:
    """Locate the PDB file corresponding to a ProteinGym protein ID."""
    if not prot_id:
        raise ValueError("Protein ID is empty; cannot locate PDB file.")

    pdb_paths = list(structures_root.rglob('*.pdb'))
    if not pdb_paths:
        raise FileNotFoundError(
            f"No PDB files found in extracted directory: {structures_root}"
        )

    exact_matches = [p for p in pdb_paths if p.stem == prot_id]
    if len(exact_matches) == 1:
        return exact_matches[0]
    if len(exact_matches) > 1:
        raise ValueError(
            f"Multiple PDB files exactly matching '{prot_id}' found: {exact_matches}"
        )

    fuzzy_matches = [p for p in pdb_paths if prot_id in p.stem]
    if not fuzzy_matches:
        raise FileNotFoundError(
            f"Could not find a PDB file for protein ID '{prot_id}' in {structures_root}"
        )
    if len(fuzzy_matches) > 1:
        raise ValueError(
            f"Multiple candidate PDB files found for protein ID '{prot_id}': {fuzzy_matches}"
        )

    return fuzzy_matches[0]


def plot_structure(
    assay_name: str,
    pdb_file: Optional[Union[str, Path]] = None,
    data_scores: str = "DMS",
    dms_data: Optional[Dict[str, pd.DataFrame]] = None,
    start_pos: Optional[int] = None,
    end_pos: Optional[int] = None,
    aggregate_fun: Callable = np.mean,
    color_scheme: Optional[str] = None
) -> Any:
    """
    Visualize DMS and model scores on 3D protein structures.
    
    Plots DMS or model scores for amino acid substitutions on a 3D protein
    structure for a chosen assay using nglview.
    
    Parameters
    ----------
    assay_name : str
        Valid DMS assay name (e.g., "C6KNH7_9INFA_Lee_2018")
    pdb_file : str or Path, optional
        Path to PDB file. If None, attempts to load from standard location.
    data_scores : str, default="DMS"
        Data source for scores. Options:
        - "DMS" for experimental DMS scores
        - Model name for zero-shot predictions
        - Supervised model name for semi-supervised predictions
    dms_data : dict, optional
        Dictionary mapping assay names to DataFrames with mutation data.
        If None, loads from standard location.
    start_pos : int, optional
        First amino acid position to plot. If None, uses minimum position.
    end_pos : int, optional
        Last amino acid position to plot. If None, uses maximum position.
    aggregate_fun : callable, default=np.mean
        Function to aggregate scores per position (e.g., np.mean, np.max, np.min)
    color_scheme : str, optional
        Color scheme for visualization. Options:
        - None: blue-white-red gradient
        - "EVE": EVE-style black-purple-cyan-yellow gradient
        
    Returns
    -------
    tuple
        (nglview.NGLWidget, matplotlib.figure.Figure)
        Interactive 3D protein structure viewer with colored residues and colorbar figure
        
    Raises
    ------
    ValueError
        If invalid assay_name or data_scores is provided
    FileNotFoundError
        If PDB file cannot be found
        
    Notes
    -----
    For model scores, a rank-based normal quantile transformation is applied
    to normalize predictions across different models. This preserves rank order
    while standardizing the distribution (mean=0, SD=1).
    
    Required columns in dms_data DataFrames:
    - 'mutant': Mutation identifier (e.g., "A1P:D2N")
    - 'DMS_score': Experimental fitness measurement
    
    Examples
    --------
    >>> from proteingympy.plot_structure import plot_structure
    >>> 
    >>> # Plot DMS scores for a specific region
    >>> view, fig = plot_structure(
    ...     assay_name="C6KNH7_9INFA_Lee_2018",
    ...     start_pos=20,
    ...     end_pos=50,
    ...     aggregate_fun=np.max
    ... )
    >>> 
    >>> # Plot zero-shot model predictions
    >>> view, fig = plot_structure(
    ...     assay_name="C6KNH7_9INFA_Lee_2018",
    ...     data_scores="GEMME",
    ...     start_pos=20,
    ...     end_pos=50
    ... )
    >>> 
    >>> # Plot with EVE color scheme
    >>> view, fig = plot_structure(
    ...     assay_name="ACE2_HUMAN_Chan_2020",
    ...     data_scores="DMS",
    ...     color_scheme="EVE"
    ... )
    """
    # Import data loading functions 
    try:
        from .make_dms_substitutions import get_dms_substitution_data
        from .make_zero_shot_substitutions import (
            get_zero_shot_scores_data,
            get_zero_shot_model_list
        )
        from .make_supervised_scores import (
            get_supervised_scores_data,
            get_supervised_model_list
        )
    except ImportError as exc:
        raise ImportError(
            "Required ProteinGym data pipeline functions not found. "
            "Please ensure the make_* modules are available."
        ) from exc
    
    # Validate data_scores
    zero_shot_models: List[str] = []
    supervised_models: List[str] = []
    try:
        zero_shot_models = get_zero_shot_model_list()
    except Exception:
        warnings.warn("Could not load zero-shot model list")
    try:
        supervised_models = get_supervised_model_list()
    except Exception:
        warnings.warn("Could not load supervised model list")

    valid_scores = ['DMS'] + zero_shot_models + supervised_models
    
    if data_scores not in valid_scores:
        raise ValueError(
            f"Invalid data_scores '{data_scores}'. "
            f"Must be 'DMS' or a valid model name from ProteinGym."
        )
    
    # Load appropriate data based on data_scores
    if data_scores == "DMS":
        if dms_data is None:
            print("'dms_data' not provided, loading with get_dms_substitution_data()")
            dms_data = get_dms_substitution_data()
        
        if dms_data is not None and assay_name not in dms_data:
            raise ValueError(f"Assay '{assay_name}' not found in dms_data")
        
        if dms_data is None:
            raise ValueError("Could not load DMS data")
        
        df = dms_data[assay_name].copy()
        df = df.rename(columns={'DMS_score': 'pg_scores'})
        
    elif data_scores in zero_shot_models:
        print(f"Using zero-shot model scores: {data_scores}")
        data = get_zero_shot_scores_data()
        if assay_name not in data:
            raise ValueError(f"Assay '{assay_name}' not found in zero-shot data")
        if data_scores not in data[assay_name].columns:
            raise ValueError(
                f"Model '{data_scores}' not available for assay '{assay_name}'"
            )
        df = data[assay_name][['mutant', data_scores]].copy()
        df = df.rename(columns={data_scores: 'pg_scores'})
        
    else:  # Supervised model
        print(f"Using semi-supervised model scores: {data_scores}")
        supervised_data, _ = get_supervised_scores_data()
        if not supervised_data:
            raise ValueError("Could not load supervised model data")
        if assay_name not in supervised_data:
            raise ValueError(
                f"Assay '{assay_name}' not found in supervised model data"
            )
        df = supervised_data[assay_name]
        if data_scores not in df.columns:
            raise ValueError(
                f"Model '{data_scores}' not available for assay '{assay_name}'"
            )
        df = df[['mutant', data_scores]].copy()
        df = df.rename(columns={data_scores: 'pg_scores'})
    
    # Load PDB file
    if pdb_file is None:
        prot_id = cast(str, get_prot_ids(assay_name))
        structures_root = _ensure_af2_structures(DEFAULT_CACHE_DIR)
        pdb_file = _find_pdb_for_prot_id(prot_id, structures_root)
    
    pdb_path = Path(pdb_file)
    if not pdb_path.exists():
        raise FileNotFoundError(f"PDB file not found: {pdb_file}")
    
    # Parse PDB structure
    if PDBParser is None:
        raise ImportError("BioPython is required for PDB parsing")
    
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', str(pdb_path))
    
    # Extract position and amino acids from mutant strings
    df['ref'] = df['mutant'].str[0]
    df['pos'] = df['mutant'].str.extract(r'(\d+)').astype(int)
    df['alt'] = df['mutant'].str[-1]
    
    # Aggregate scores by position
    df_agg = df.groupby('pos').agg(
        aggregate_score=('pg_scores', aggregate_fun)
    ).reset_index()
    
    # Filter by position range
    filtered_df = filter_by_pos(df_agg, start_pos, end_pos)
    start_pos = int(filtered_df['pos'].min())
    end_pos = int(filtered_df['pos'].max())
    
    # Prepare color mapping
    if data_scores == "DMS":
        if color_scheme == "EVE":
            min_val = filtered_df['aggregate_score'].min()
            max_val = filtered_df['aggregate_score'].max()
            mid1 = min_val + (max_val - min_val) / 3
            mid2 = min_val + (max_val - min_val) * 2 / 3
            values = [min_val, mid1, mid2, max_val]
        else:
            values = [
                filtered_df['aggregate_score'].min(),
                0,
                filtered_df['aggregate_score'].max()
            ]
        
        col_fun = get_color_function(color_scheme, values)
        filtered_df['color'] = filtered_df['aggregate_score'].apply(col_fun)
        
    else:  # Model scores - use quantile normalization
        if color_scheme == "EVE":
            col_palette = None  # Will use default in color_line
        else:
            col_palette = _create_parula_palette(n=200)
        
        filtered_df = color_line(
            filtered_df,
            quant_norm=True,
            color_palette=col_palette
        )
        
        # Apply EVE coloring if specified
        if color_scheme == "EVE":
            min_val = filtered_df['quant_clamped'].min()
            max_val = filtered_df['quant_clamped'].max()
            mid1 = min_val + (max_val - min_val) / 3
            mid2 = min_val + (max_val - min_val) * 2 / 3
            values = [min_val, mid1, mid2, max_val]
            col_fun = get_color_function("EVE", values)
            filtered_df['color'] = filtered_df['quant_clamped'].apply(col_fun)
    
    # Create nglview widget
    if nv is None:
        raise ImportError("nglview is required for 3D visualization")
    
    from nglview.color import ColormakerRegistry
    
    view = nv.show_file(str(pdb_path), default=False)
    view.stage.set_parameters(**{
        "clipNear": 0, 
        "clipFar": 100, 
        "clipDist": 10,
        "fogNear": 0, 
        "fogFar": 1000,
        "backgroundColor": "white",
    })
    
    # Build color scheme as list of [color, selection] pairs
    color_scheme_list = []
    for _, row in filtered_df.iterrows():
        pos = int(row['pos'])
        color = row['color']
        color_scheme_list.append([color, str(pos)])
    
    # Register the custom color scheme
    scheme_id = ColormakerRegistry.add_selection_scheme(
        "custom_colors", 
        color_scheme_list
    )
    
    # Add cartoon with the custom color scheme
    view.add_cartoon(selection="protein", color="custom_colors")
    
    # Center view on selected region
    view.center()
    
    # Create colorbar
    try:
        import matplotlib.pyplot as plt
        from matplotlib.colors import LinearSegmentedColormap, Normalize
        
        fig, ax = plt.subplots(figsize=(6, 0.6))
        fig.subplots_adjust(bottom=0.5)
        
        if data_scores == "DMS":
            # Use the actual color scheme
            if color_scheme == "EVE":
                # EVE colors: black -> purple -> cyan -> yellow
                colors_list = ['#000000', '#9440e8', '#00CED1', '#fde662']
                n_bins = 100
                cmap = LinearSegmentedColormap.from_list('eve', colors_list, N=n_bins)
                vmin = filtered_df['aggregate_score'].min()
                vmax = filtered_df['aggregate_score'].max()
                label = 'DMS Score'
            else:
                # Default: red -> white -> blue
                colors_list = ['#ff0000', '#ffffff', '#0000ff']
                n_bins = 100
                cmap = LinearSegmentedColormap.from_list('default', colors_list, N=n_bins)
                vmin = filtered_df['aggregate_score'].min()
                vmax = filtered_df['aggregate_score'].max()
                label = 'DMS Score'
        else:
            # Model scores use quantile normalization
            if color_scheme == "EVE":
                colors_list = ['#000000', '#9440e8', '#00CED1', '#fde662']
                n_bins = 100
                cmap = LinearSegmentedColormap.from_list('eve', colors_list, N=n_bins)
                vmin = filtered_df['quant_clamped'].min()
                vmax = filtered_df['quant_clamped'].max()
                label = f'{data_scores} Score (Quantile Normalized)'
            else:
                cmap = plt.get_cmap('viridis')
                vmin = filtered_df['quant_clamped'].min()
                vmax = filtered_df['quant_clamped'].max()
                label = f'{data_scores} Score (Quantile Normalized)'
        
        norm = Normalize(vmin=vmin, vmax=vmax)
        fig.colorbar(
            plt.cm.ScalarMappable(norm=norm, cmap=cmap),
            cax=ax, 
            orientation='horizontal', 
            label=label
        )
        
    except ImportError:
        warnings.warn("matplotlib not available, skipping colorbar generation")
        fig = None
    
    return view, fig
