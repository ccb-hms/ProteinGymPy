"""
plot_dms_heatmap.py - Python equivalent of plot_dms_heatmap.R

Visualize DMS scores along a protein as a heatmap.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap

from typing import Optional, Dict, Tuple, Any


def _filter_by_pos(
    df: pd.DataFrame, 
    start_pos: Optional[int] = None, 
    end_pos: Optional[int] = None
) -> pd.DataFrame:
    """
    Filter dataframe by position range.
    
    Args:
        df: DataFrame with 'pos' column
        start_pos: First amino acid position (default: min position)
        end_pos: Last amino acid position (default: max position)
        
    Returns:
        Filtered DataFrame
        
    Raises:
        ValueError: If pos column is missing or invalid, or positions out of range
    """
    # Check pos column
    if 'pos' not in df.columns:
        raise ValueError("The dataframe must contain a 'pos' column.")
    
    if not pd.api.types.is_integer_dtype(df['pos']):
        raise ValueError("The 'pos' column must be an integer type.")
    
    # Grab minimum and maximum values of the pos column
    min_pos = df['pos'].min()
    max_pos = df['pos'].max()
    
    # Check if user-provided start_pos or end_pos is within the range
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
    
    # If start or end is None, default to min or max "pos"
    if start_pos is None:
        start_pos = min_pos
    if end_pos is None:
        end_pos = max_pos
    
    # Filter the dataframe based on the specified positions
    filtered_df = df[(df['pos'] >= start_pos) & (df['pos'] <= end_pos)].copy()
    
    return filtered_df


def _filter_exact_coord(
    assay_pos: pd.DataFrame,
    start_pos: Optional[int] = None,
    end_pos: Optional[int] = None,
    exact_coord: Optional[bool] = None
) -> pd.DataFrame:
    """
    Optionally fill in missing positions with NA values.
    
    Args:
        assay_pos: DataFrame with position data
        start_pos: First position
        end_pos: Last position
        exact_coord: If True, create consecutive positions; if False/None, use only available positions
        
    Returns:
        DataFrame with optional filled positions
    """
    if exact_coord is None:
        print(
            "'exact_coord' not provided, "
            "using only positions available in assay."
        )
        return assay_pos
    
    elif exact_coord is False:
        return assay_pos
    
    elif exact_coord is True:
        if start_pos is None:
            start_pos = int(assay_pos['pos'].min())
        if end_pos is None:
            end_pos = int(assay_pos['pos'].max())
        
        # Create a sequence of consecutive positions
        all_pos = pd.DataFrame({'pos': range(start_pos, end_pos + 1)})
        
        # Merge with full sequence and fill missing values with NA
        assay_pos = all_pos.merge(assay_pos, on='pos', how='left')
        
        return assay_pos
    
    else:
        return assay_pos


def _make_colormap_dms(mat: np.ndarray, color_scheme: str = "default") -> Tuple[Any, float, float]:
    """
    Create a colormap for the heatmap.
    
    Args:
        mat: Matrix of DMS scores
        color_scheme: Either "default" (red-white-blue) or "EVE" (popEVE portal colors)
        
    Returns:
        Tuple of (colormap, vmin, vmax)
    """
    vmin = np.nanmin(mat)
    vmax = np.nanmax(mat)
    
    if color_scheme == "EVE":
        # EVE color scheme: black -> purple -> cyan -> yellow
        halfpt = vmin / 2
        colors = ["#000000", "#9440e8", "#00CED1", "#fde662"]
        positions = [0, (halfpt - vmin) / (vmax - vmin), (0 - vmin) / (vmax - vmin), 1]
        cmap = LinearSegmentedColormap.from_list(
            "eve", list(zip(positions, colors))
        )
    else:
        # Default color scheme: red-white-blue
        cmap = LinearSegmentedColormap.from_list(
            "default", ["red", "white", "blue"]
        )
    
    return cmap, vmin, vmax


def plot_dms_heatmap(
    assay_name: str,
    dms_data: Optional[Dict[str, pd.DataFrame]] = None,
    start_pos: Optional[int] = None,
    end_pos: Optional[int] = None,
    exact_coord: bool = False,
    cluster_rows: bool = False,
    cluster_columns: bool = False,
    color_scheme: str = "default",
    figsize: Tuple[float, float] = (12, 8),
    **kwargs
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Visualize DMS scores along a protein.
    
    Plots DMS scores for amino acid substitutions along a protein in a defined 
    DMS assay. The x-axis shows amino acid positions where DMS mutations exist, 
    and the y-axis represents possible amino acid residues, ordered by default 
    based on physiochemical groupings.
    
    Args:
        assay_name: Valid DMS assay name (key in dms_data dictionary)
        dms_data: Dictionary of DMS assays with DataFrames. If None, attempts to 
            load from proteingympy.make_dms_substitutions.get_dms_substitution_data()
        start_pos: First amino acid position to plot (default: first available position)
        end_pos: Last amino acid position to plot (default: last available position)
        exact_coord: If True, plot precise start_pos and end_pos coordinates, filling 
            missing positions with NaN. If False, plot only positions with available data.
        cluster_rows: If True, cluster amino acid rows (requires scipy)
        cluster_columns: If True, cluster position columns (requires scipy)
        color_scheme: "default" for red-white-blue or "EVE" for popEVE portal colors
        figsize: Figure size as (width, height) in inches
        **kwargs: Additional arguments passed to seaborn.heatmap()
        
    Returns:
        Tuple of (figure, axes) objects
        
    Raises:
        ValueError: If assay contains only multiple amino acid sites, or invalid parameters
        ImportError: If required packages are missing
        
    Examples:
        >>> from proteingympy.make_dms_substitutions import get_dms_substitution_data
        >>> dms_data = get_dms_substitution_data()
        >>> fig, ax = plot_dms_heatmap(
        ...     assay_name="A0A192B1T2_9HIV1_Haddox_2018",
        ...     dms_data=dms_data,
        ...     start_pos=10,
        ...     end_pos=80
        ... )
        >>> plt.show()
        
        >>> # With EVE color scheme and exact coordinates
        >>> fig, ax = plot_dms_heatmap(
        ...     assay_name="A0A192B1T2_9HIV1_Haddox_2018",
        ...     dms_data=dms_data,
        ...     start_pos=10,
        ...     end_pos=80,
        ...     exact_coord=True,
        ...     color_scheme="EVE"
        ... )
        >>> plt.show()
    """
    # Check dependencies for clustering
    if cluster_rows or cluster_columns:
        try:
            from scipy.cluster import hierarchy
            from scipy.spatial.distance import pdist
        except ImportError:
            raise ImportError(
                "Clustering requires scipy. Install with: pip install scipy"
            )
    
    # Load dms_data if not provided
    if dms_data is None:
        print(
            "'dms_data' not provided, "
            "using DMS data loaded with get_dms_substitution_data()"
        )
        try:
            from .make_dms_substitutions import get_dms_substitution_data
            dms_data = get_dms_substitution_data()
        except ImportError:
            raise ImportError(
                "Could not import get_dms_substitution_data. "
                "Please provide dms_data parameter."
            )
    
    # Extract the specified assay
    if assay_name not in dms_data:
        raise ValueError(
            f"Assay '{assay_name}' not found in dms_data. "
            f"Available assays: {list(dms_data.keys())[:5]}..."
        )
    
    assay_df = dms_data[assay_name].copy()
    
    # Filter out multiple aa sites (those containing ':')
    assay_df = assay_df[~assay_df['mutant'].str.contains(':', na=False)]
    
    # Stop if all rows are multiple sites
    if len(assay_df) == 0:
        raise ValueError(
            f"Unable to plot DMS substitution heatmap; "
            f"assay '{assay_name}' contains only multiple amino acid sites."
        )
    
    # Wrangle the data: extract ref, pos, alt from mutant string
    # Mutant format: A1P (ref=A, pos=1, alt=P)
    assay_df['ref'] = assay_df['mutant'].str[0]
    assay_df['pos'] = assay_df['mutant'].str.extract(r'(\d+)')[0].astype(int)
    assay_df['alt'] = assay_df['mutant'].str[-1]
    
    # Select relevant columns
    assay_df = assay_df[['ref', 'pos', 'alt', 'DMS_score']]
    
    # Reshape to wide format
    assay_wide = assay_df.pivot_table(
        index='pos',
        columns='alt',
        values='DMS_score',
        aggfunc='first'  # In case of duplicates, take first
    ).reset_index()
    
    # Also keep ref column
    ref_by_pos = assay_df[['pos', 'ref']].drop_duplicates().set_index('pos')
    assay_wide = assay_wide.join(ref_by_pos, on='pos')
    
    # Subset to start_pos and end_pos
    if start_pos is None:
        print("'start_pos' not provided, using the first position in the protein.")
    
    if end_pos is None:
        print("'end_pos' not provided, using the last position in the protein.")
    
    assay_pos = _filter_by_pos(
        df=assay_wide,
        start_pos=start_pos,
        end_pos=end_pos
    )
    
    # Apply exact_coord filtering
    assay_pos = _filter_exact_coord(
        assay_pos,
        start_pos=start_pos,
        end_pos=end_pos,
        exact_coord=exact_coord
    )
    
    # Get column annotation (reference amino acids)
    column_annotation = assay_pos[['ref', 'pos']].drop_duplicates()
    
    # Check for clustering with NA values
    if column_annotation['ref'].isna().any() and cluster_columns:
        raise ValueError(
            "Protein range includes missing values, preventing clustering of columns. "
            "Try setting exact_coord argument to False."
        )
    
    # Fill NaN in annotations with space
    column_annotation['ref'] = column_annotation['ref'].fillna(' ')
    
    # Convert to matrix
    pos_values = assay_pos['pos'].values
    
    # Get amino acid columns (exclude 'pos' and 'ref')
    aa_cols = [col for col in assay_pos.columns if col not in ['pos', 'ref']]
    
    heatmap_matrix = assay_pos[aa_cols].values.T
    
    # Reorder rows based on physiochemical properties
    physiochem_order = "DEKRHNQSTPGAVILMCFYW"
    aa_order = list(physiochem_order)
    
    # Reorder matrix rows
    row_labels = aa_cols
    reordered_indices = []
    reordered_labels = []
    
    for aa in aa_order:
        if aa in row_labels:
            idx = row_labels.index(aa)
            reordered_indices.append(idx)
            reordered_labels.append(aa)
    
    reordered_matrix = heatmap_matrix[reordered_indices, :]
    
    # Apply clustering if requested
    if cluster_rows:
        # Remove rows that are all NaN
        valid_rows = ~np.all(np.isnan(reordered_matrix), axis=1)
        if valid_rows.sum() > 1:
            from scipy.cluster import hierarchy
            from scipy.spatial.distance import pdist
            
            valid_matrix = reordered_matrix[valid_rows]
            # Replace NaN with 0 for distance calculation
            filled_matrix = np.nan_to_num(valid_matrix, nan=0)
            
            if filled_matrix.shape[0] > 1:
                row_linkage = hierarchy.linkage(pdist(filled_matrix), method='average')
                row_dendro = hierarchy.dendrogram(row_linkage, no_plot=True)
                row_order = row_dendro['leaves']
                
                # Apply ordering to valid rows only
                temp_matrix = valid_matrix[row_order]
                temp_labels = [reordered_labels[i] for i, v in enumerate(valid_rows) if v]
                temp_labels = [temp_labels[i] for i in row_order]
                
                reordered_matrix = temp_matrix
                reordered_labels = temp_labels
    
    if cluster_columns:
        # Remove columns that are all NaN
        valid_cols = ~np.all(np.isnan(reordered_matrix), axis=0)
        if valid_cols.sum() > 1:
            from scipy.cluster import hierarchy
            from scipy.spatial.distance import pdist
            
            valid_matrix = reordered_matrix[:, valid_cols]
            # Replace NaN with 0 for distance calculation
            filled_matrix = np.nan_to_num(valid_matrix, nan=0)
            
            if filled_matrix.shape[1] > 1:
                col_linkage = hierarchy.linkage(pdist(filled_matrix.T), method='average')
                col_dendro = hierarchy.dendrogram(col_linkage, no_plot=True)
                col_order = col_dendro['leaves']
                
                # Apply ordering to valid columns only
                temp_matrix = valid_matrix[:, col_order]
                temp_pos = pos_values[valid_cols][col_order]
                temp_ref = column_annotation['ref'].values[valid_cols][col_order]
                
                reordered_matrix = temp_matrix
                pos_values = temp_pos
                column_annotation = pd.DataFrame({'ref': temp_ref, 'pos': temp_pos})
    
    # Create colormap
    cmap, vmin, vmax = _make_colormap_dms(reordered_matrix, color_scheme)
    
    # Create the heatmap
    fig, ax = plt.subplots(figsize=figsize)
    
    # Set up heatmap arguments
    heatmap_kwargs = {
        'cmap': cmap,
        'vmin': vmin,
        'vmax': vmax,
        'cbar_kws': {'label': 'DMS Score'},
        'xticklabels': pos_values,
        'yticklabels': reordered_labels,
        'linewidths': 0,
        'square': False,
    }
    
    # Update with user-provided kwargs
    heatmap_kwargs.update(kwargs)
    
    # Create heatmap
    sns.heatmap(
        reordered_matrix,
        ax=ax,
        **heatmap_kwargs
    )
    
    # Add reference amino acid annotations on top
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(np.arange(len(pos_values)) + 0.5)
    ax2.set_xticklabels(column_annotation['ref'].values, fontsize=10)
    ax2.tick_params(length=0)
    
    # Labels
    ax.set_xlabel('Position', fontsize=12)
    ax.set_ylabel('Amino Acid', fontsize=12)
    ax2.set_xlabel('Reference AA', fontsize=10, labelpad=10)
    
    plt.tight_layout()
    
    return fig, ax
