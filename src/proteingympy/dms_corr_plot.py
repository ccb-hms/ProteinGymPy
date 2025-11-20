"""
dms_corr_plot.py - Python equivalent of dms_corr_plot.R

Compare DMS experimental scores with model predictions using Spearman correlation
and visualization with matplotlib.
"""

from typing import Optional, Tuple, Dict
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec

try:
    import seaborn as sns
except ImportError:
    sns = None

try:
    from scipy import stats
except ImportError:
    stats = None

from .make_dms_substitutions import get_dms_substitution_data
from .make_alphamissense_supplementary import get_alphamissense_proteingym_data
from .make_zero_shot_substitutions import (
    get_zero_shot_substitution_data,
    available_zero_shot_models
)
from .make_supervised_scores import (
    get_supervised_substitution_data,
    available_supervised_models
)


def _filter_alphamissense_table(
    am_table: Optional[pd.DataFrame],
    uniprot_id: str,
    cache_dir: str = ".cache"
) -> pd.DataFrame:
    """
    Filter AlphaMissense table for a specific UniProt ID.

    Args:
        am_table: AlphaMissense DataFrame (if None, will load default)
        uniprot_id: UniProt accession identifier
        cache_dir: Cache directory for downloading data

    Returns:
        Filtered DataFrame for the specified protein
    """
    if am_table is None:
        # Load default AlphaMissense data
        am_table = get_alphamissense_proteingym_data(cache_dir=cache_dir)

        # Rename columns to match model_tables convention
        am_table = am_table.rename(columns={
            'Uniprot_ID': 'UniProt_id',
            'variant_id': 'mutant'
        })

    # Filter for uniprot_id
    filtered_table = am_table[am_table['UniProt_id'] == uniprot_id].copy()

    # Check if table is empty after filtering
    if len(filtered_table) == 0:
        raise ValueError(
            f"No AlphaMissense information found for protein accession '{uniprot_id}'. "
            f"Check that the UniProt ID is correct."
        )

    return filtered_table


def _filter_model_table(
    model_table: Dict[str, pd.DataFrame],
    uniprot_id: str
) -> pd.DataFrame:
    """
    Filter model table (zero-shot or supervised) for a specific UniProt ID.

    Args:
        model_table: Dictionary of DataFrames (one per DMS assay)
        uniprot_id: UniProt accession identifier

    Returns:
        Combined DataFrame for the specified protein
    """
    # Extract assays containing uniprot_id
    model_list = []
    for study_name, df in model_table.items():
        if df is not None and 'UniProt_id' in df.columns:
            if (df['UniProt_id'] == uniprot_id).any():
                model_list.append(df)

    # Error if no assay found for uniprot_id
    if len(model_list) == 0:
        raise ValueError(
            f"No ProteinGym assay found for protein accession '{uniprot_id}'. "
            f"Check that the UniProt ID is correct."
        )

    # Combine into one dataframe
    combined_table = pd.concat(model_list, ignore_index=True)

    return combined_table


def _filter_dms_uniprot(
    dms_table: Optional[Dict[str, pd.DataFrame]],
    uniprot_id: str,
    cache_dir: str = ".cache"
) -> pd.DataFrame:
    """
    Filter DMS table for a specific UniProt ID.

    Args:
        dms_table: Dictionary of DMS DataFrames (if None, will load default)
        uniprot_id: UniProt accession identifier
        cache_dir: Cache directory for downloading data

    Returns:
        Combined DataFrame for the specified protein with DMS scores
    """
    if dms_table is None:
        print("'dms_table' not provided, using default table from get_dms_substitution_data()")
        dms_table = get_dms_substitution_data(cache_dir=cache_dir)

    # Filter for uniprot_id - remove any tables without UniProt_id or with NaN
    dms_list = []
    for study_name, df in dms_table.items():
        if df is not None and 'UniProt_id' in df.columns:
            # Check if UniProt_id column has the target ID
            if (df['UniProt_id'] == uniprot_id).any():
                dms_list.append(df)

    # Check if table is empty after filtering
    if len(dms_list) == 0:
        raise ValueError(
            f"No DMS substitution information found for protein accession '{uniprot_id}'. "
            f"Check that the UniProt ID is correct."
        )

    # Combine into one dataframe
    combined_table = pd.concat(dms_list, ignore_index=True)

    return combined_table


def _get_model_dataframe(
    model: str,
    uniprot_id: str,
    cache_dir: str = ".cache"
) -> pd.DataFrame:
    """
    Get the appropriate model dataframe based on model name.

    Args:
        model: Model name (e.g., "AlphaMissense", "GEMME", "ESM_1v")
        uniprot_id: UniProt accession identifier
        cache_dir: Cache directory for downloading data

    Returns:
        DataFrame with columns: DMS_id, UniProt_id, mutant, [model_name]
    """
    zero_shot_models = available_zero_shot_models()
    supervised_models = available_supervised_models()

    if model == "AlphaMissense":
        print(f"Using AlphaMissense model from get_alphamissense_proteingym_data()")
        model_table = _filter_alphamissense_table(
            am_table=None,
            uniprot_id=uniprot_id,
            cache_dir=cache_dir
        )
        # Select relevant columns
        cols = ['DMS_id', 'UniProt_id', 'mutant', 'AlphaMissense']
        available_cols = [col for col in cols if col in model_table.columns]
        model_table = model_table[available_cols]

    elif model in zero_shot_models:
        # Load zero-shot data
        zeroshot_data = get_zero_shot_substitution_data(cache_dir=cache_dir)
        model_table = _filter_model_table(
            model_table=zeroshot_data,
            uniprot_id=uniprot_id
        )
        # Select relevant columns
        cols = ['DMS_id', 'UniProt_id', 'mutant', model]
        available_cols = [col for col in cols if col in model_table.columns]
        model_table = model_table[available_cols]

    elif model in supervised_models:
        # Load supervised data (default to random_5)
        supervised_data, _ = get_supervised_substitution_data(
            fold_type="random_5",
            cache_dir=cache_dir
        )
        model_table = _filter_model_table(
            model_table=supervised_data,
            uniprot_id=uniprot_id
        )
        # Select relevant columns
        cols = ['DMS_id', 'UniProt_id', 'mutant', model]
        available_cols = [col for col in cols if col in model_table.columns]
        model_table = model_table[available_cols]

    else:
        raise ValueError(f"Model '{model}' not recognized.")

    return model_table


def _merge_model_dms_tables(
    model_df: pd.DataFrame,
    dms_df: pd.DataFrame,
    model: str
) -> pd.DataFrame:
    """
    Merge model and DMS dataframes by UniProt_id and mutant.

    Args:
        model_df: Model DataFrame
        dms_df: DMS DataFrame
        model: Model name

    Returns:
        Merged DataFrame with averaged scores per mutant
    """
    # Check that UniProt IDs are the same
    uniprot_model = model_df['UniProt_id'].unique()
    uniprot_dms = dms_df['UniProt_id'].unique()

    if len(uniprot_model) != 1 or len(uniprot_dms) != 1 or uniprot_model[0] != uniprot_dms[0]:
        raise ValueError("UniProt IDs must be the same across model and DMS tables")

    # Merge tables
    merged_table = pd.merge(
        model_df,
        dms_df,
        on=['UniProt_id', 'mutant'],
        how='left',
        suffixes=('_model', '_dms')
    )

    # Handle DMS_id columns
    if 'DMS_id_model' in merged_table.columns and 'DMS_id_dms' in merged_table.columns:
        merged_table['DMS_id'] = merged_table['DMS_id_dms']
        merged_table = merged_table.drop(columns=['DMS_id_model', 'DMS_id_dms'])

    # Handle DMS_score columns (keep the DMS one if both exist)
    if 'DMS_score_model' in merged_table.columns and 'DMS_score_dms' in merged_table.columns:
        merged_table['DMS_score'] = merged_table['DMS_score_dms']
        merged_table = merged_table.drop(columns=['DMS_score_model', 'DMS_score_dms'])
    elif 'DMS_score_dms' in merged_table.columns:
        merged_table['DMS_score'] = merged_table['DMS_score_dms']
        merged_table = merged_table.drop(columns=['DMS_score_dms'])

    # Select relevant columns and remove NaN
    cols = ['UniProt_id', 'mutant', model, 'DMS_score']
    available_cols = [col for col in cols if col in merged_table.columns]
    merged_table = merged_table[available_cols].dropna()

    # Average model and DMS scores across multiple studies per protein
    merged_table = merged_table.groupby(['UniProt_id', 'mutant'], as_index=False).agg({
        model: 'mean',
        'DMS_score': 'mean'
    })

    merged_table = merged_table.rename(columns={
        model: 'mean_model',
        'DMS_score': 'mean_dms'
    })

    return merged_table


def _calculate_spearman_correlation(
    merged_table: pd.DataFrame
) -> Tuple[float, float]:
    """
    Calculate Spearman correlation between model scores and DMS scores.

    Args:
        merged_table: Merged DataFrame with mean_model and mean_dms columns

    Returns:
        Tuple of (correlation coefficient, p-value)
    """
    if stats is None:
        raise ImportError("scipy is required for correlation analysis. Install with: pip install scipy")

    correlation, pvalue = stats.spearmanr(
        merged_table['mean_model'],
        merged_table['mean_dms']
    )

    return correlation, pvalue


def dms_corr_plot(
    uniprot_id: str,
    model: str = "AlphaMissense",
    dms_table: Optional[Dict[str, pd.DataFrame]] = None,
    cache_dir: str = ".cache",
    figsize: Tuple[float, float] = (10, 8),
    bins: int = 60,
    cmap: str = 'viridis'
) -> Figure:
    """
    Compare DMS experimental scores with model predictions using Spearman correlation.

    Creates a 2D density plot (hexbin) with marginal distributions showing the
    relationship between experimental DMS scores and predicted model scores.

    Args:
        uniprot_id: Valid UniProt accession identifier
        model: Model to plot (default: "AlphaMissense")
        dms_table: Optional dict of DMS DataFrames (if None, loads default)
        cache_dir: Directory to cache downloaded files
        figsize: Figure size as (width, height)
        bins: Number of bins for hexbin plot
        cmap: Colormap for hexbin plot

    Returns:
        matplotlib Figure object

    Examples:
        >>> # Use defaults (AlphaMissense vs DMS scores)
        >>> fig = dms_corr_plot(uniprot_id="Q9NV35")
        >>> plt.show()

        >>> # Compare specific model with DMS
        >>> fig = dms_corr_plot(
        ...     uniprot_id="P04637",
        ...     model="Kermut"
        ... )
        >>> plt.show()

    Notes:
        - Requires scipy for correlation analysis
        - Requires matplotlib for visualization
        - seaborn is optional but recommended for better styling
        - Model names can be from zero-shot, supervised, or AlphaMissense
        - Use available_zero_shot_models() and available_supervised_models()
          to see available models

    """
    # Check dependencies
    if stats is None:
        raise ImportError(
            "scipy is required for correlation analysis. "
            "Install with: pip install scipy"
        )

    # Validate uniprot_id
    if not isinstance(uniprot_id, str):
        raise ValueError("uniprot_id must be a string")

    # Validate model
    valid_models = (
        available_zero_shot_models() +
        available_supervised_models() +
        ["AlphaMissense"]
    )

    if model not in valid_models:
        raise ValueError(f"Invalid model specified: {model}")

    # Load model data for uniprot_id
    print(f"Loading data for {model}...")
    model_df = _get_model_dataframe(model, uniprot_id, cache_dir)

    # Load DMS data for uniprot_id
    print(f"Loading DMS scores for {uniprot_id}...")
    dms_df = _filter_dms_uniprot(dms_table, uniprot_id, cache_dir)

    # Merge tables
    print("Merging model predictions with DMS scores...")
    merged_table = _merge_model_dms_tables(
        model_df=model_df,
        dms_df=dms_df,
        model=model
    )

    # Check if merged table is empty
    if len(merged_table) == 0:
        raise ValueError(
            f"No common mutants between {model} and DMS scores for accession '{uniprot_id}'"
        )

    print(f"Found {len(merged_table)} common mutants")

    # Calculate correlation
    correlation, pvalue = _calculate_spearman_correlation(merged_table)

    print(f"Spearman r = {correlation:.3f}, p-value = {pvalue:.2e}")

    # Create figure with GridSpec for layout
    fig = plt.figure(figsize=figsize)
    gs = GridSpec(
        3, 3,
        figure=fig,
        hspace=0.05,
        wspace=0.05,
        width_ratios=[1, 4, 0.3],
        height_ratios=[1, 4, 0.5]
    )

    # Main plot (center)
    ax_main = fig.add_subplot(gs[1, 1])

    # Marginal plots
    ax_top = fig.add_subplot(gs[0, 1], sharex=ax_main)
    ax_right = fig.add_subplot(gs[1, 2], sharey=ax_main)

    # Colorbar axis
    ax_cbar = fig.add_subplot(gs[2, 1])

    # Create hexbin plot (x=DMS, y=model - matching R version)
    x = merged_table['mean_dms'].values
    y = merged_table['mean_model'].values

    hexbin = ax_main.hexbin(
        x, y,
        gridsize=bins,
        cmap=cmap,
        mincnt=1,
        edgecolors='none'
    )

    # Style main plot
    ax_main.set_xlabel('DMS score', fontsize=16)
    ax_main.set_ylabel(f'{model} score', fontsize=16)
    ax_main.tick_params(labelsize=16)

    # Add colorbar
    cbar = plt.colorbar(hexbin, cax=ax_cbar, orientation='horizontal')
    cbar.set_label('Count', fontsize=16)
    cbar.ax.tick_params(labelsize=16)

    # Top marginal (histogram + KDE for DMS scores)
    ax_top.hist(x, bins=30, color='#B0C4DE', edgecolor='black', alpha=0.7, density=True)
    if sns is not None:
        try:
            from scipy.stats import gaussian_kde
            kde = gaussian_kde(x)
            x_range = np.linspace(x.min(), x.max(), 100)
            ax_top.plot(x_range, kde(x_range), 'k-', linewidth=1.5)
        except:
            pass  # Skip KDE if it fails
    ax_top.set_ylabel('Density', fontsize=16)
    ax_top.tick_params(labelbottom=False, labelsize=16)
    ax_top.spines['top'].set_visible(False)
    ax_top.spines['right'].set_visible(False)

    # Right marginal (histogram + KDE for model scores)
    ax_right.hist(y, bins=30, color='#B0C4DE', edgecolor='black',
                  alpha=0.7, orientation='horizontal', density=True)
    if sns is not None:
        try:
            from scipy.stats import gaussian_kde
            kde = gaussian_kde(y)
            y_range = np.linspace(y.min(), y.max(), 100)
            ax_right.plot(kde(y_range), y_range, 'k-', linewidth=1.5)
        except:
            pass  # Skip KDE if it fails
    ax_right.set_xlabel('Density', fontsize=16)
    ax_right.tick_params(labelleft=False, labelsize=16)
    ax_right.spines['top'].set_visible(False)
    ax_right.spines['right'].set_visible(False)

    # Add title with correlation info
    title = (
        f"UniProt ID: {uniprot_id}\n"
        f"Spearman r = {correlation:.2f}, p-value = {pvalue:.2e}"
    )
    fig.suptitle(title, fontsize=14, y=0.98)

    return fig


def get_available_models() -> Dict[str, list]:
    """
    Get dictionary of all available models in ProteinGym.

    Returns:
        Dictionary with keys 'zero_shot', 'supervised', and 'other'
    """
    return {
        'zero_shot': available_zero_shot_models(),
        'supervised': available_supervised_models(),
        'other': ['AlphaMissense']
    }


if __name__ == "__main__":
    # Example usage
    print("DMS Correlation Plot Example")
    print("=" * 50)

    # Show available models
    models = get_available_models()
    print("\nAvailable models:")
    print(f"  Zero-shot ({len(models['zero_shot'])}): {models['zero_shot'][:5]}...")
    print(f"  Supervised ({len(models['supervised'])}): {models['supervised']}")
    print(f"  Other: {models['other']}")

    # Example usage (commented out to avoid actual execution)
    print("\nExample usage:")
    print("  fig = dms_corr_plot(uniprot_id='Q9NV35')")
    print("  plt.show()")

