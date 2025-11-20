"""
benchmark_models.py - Python equivalent of benchmark_models.R

Benchmark Variant Effect Prediction Models. `benchmark_models()` plots one of the five model performance metrics ("AUC", "MCC", "NDCG", "Spearman", "Top_recall") for up to 5 user-specified variant effect prediction tools listed in `available_models()`. 
"""

from typing import List, Dict, Optional, Any
import os
import pickle
import warnings

try:
    import numpy as np
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
except Exception as e:
    raise ImportError(
        "Required packages not found. Please install pandas, numpy, seaborn, and matplotlib."
        " Example: pip install pandas numpy seaborn matplotlib"
    ) from e


def available_models():
    """
    Returns a list of available model names for users to choose from.
    Equivalent to the R function `available_models()`.
    """
    models = [
       'Site_Independent', 'EVmutation', 'DeepSequence_single',
       'DeepSequence_ensemble', 'EVE_single', 'EVE_ensemble', 'Unirep',
       'Unirep_evotuned', 'MSA_Transformer_single', 'MSA_Transformer_ensemble',
       'ESM_1b', 'ESM_1v_single', 'ESM_1v_ensemble', 'ESM2_8M', 'ESM2_35M',
       'ESM2_150M', 'ESM2_650M', 'ESM2_3B', 'ESM2_15B', 'Wavenet', 'RITA_S',
       'RITA_M', 'RITA_L', 'RITA_XL', 'Progen2_S', 'Progen2_M', 'Progen2_Base',
       'Progen2_L', 'Progen2_XL', 'GEMME', 'VESPA', 'VESPAl', 'VespaG',
       'ProtGPT2', 'Tranception_S_no_retrieval', 'Tranception_M_no_retrieval',
       'Tranception_L_no_retrieval', 'Tranception_S', 'Tranception_M',
       'Tranception_L', 'TranceptEVE_S', 'TranceptEVE_M', 'TranceptEVE_L',
       'CARP_38M', 'CARP_600K', 'CARP_640M', 'CARP_76M', 'MIF', 'MIF_ST',
       'ESM_IF1', 'ProteinMPNN', 'ProtSSN_k_10_h_512', 'ProtSSN_k_10_h_768',
       'ProtSSN_k_10_h_1280', 'ProtSSN_k_20_h_512', 'ProtSSN_k_20_h_768',
       'ProtSSN_k_20_h_1280', 'ProtSSN_k_30_h_512', 'ProtSSN_k_30_h_768',
       'ProtSSN_k_30_h_1280', 'ProtSSN_ensemble', 'SaProt_650M', 'SaProt_35M',
       'PoET_200M', 'MULAN', 'ProSST_K_20', 'ProSST_K_128', 'ProSST_K_512',
       'ProSST_K_1024', 'ProSST_K_2048', 'ProSST_K_4096', 'ESCOTT', 'VenusREM',
       'RSALOR', 'S2F', 'S2F_MSA', 'S3F', 'S3F_MSA', 'SiteRM'
    ]
    return models


from typing import List, Union

def normalize_model_names(models: Union[str, List[str]]) -> List[str]:
    """Normalize model names to match DataFrame / available_models format."""
    if isinstance(models, str):
        models = [models]
    return [m.strip().replace(" ", "_").lower() for m in models]

from typing import Iterable, List, Union

def check_metric_argument(user_metric: str) -> None:
    """
    Validate the metric argument.
    """
    valid_metrics = ["AUC", "MCC", "NDCG", "Spearman", "Top_recall"]

    # Ensure scalar string
    if not isinstance(user_metric, str):
        raise TypeError(
            f"Metric must be a single string, not {type(user_metric).__name__}: {user_metric}"
        )

    if user_metric not in valid_metrics:
        raise ValueError(
            f"Invalid metric specified: {user_metric}. "
            f"Valid metrics are: {', '.join(valid_metrics)}"
        )



from typing import List, Union

def check_model_argument(models: Union[str, List[str]]) -> None:

    models = normalize_model_names(models)

    # Accept either a string or a list
    if isinstance(models, str):
        models = [models]

    # Combine both sets of valid models
    valid_models = available_models()
    valid_models = normalize_model_names(valid_models)
    
    # Check validity
    invalid_models = [m for m in models if m not in valid_models]
    if invalid_models:
        raise ValueError(f"Invalid model(s) specified: {invalid_models}")

    # Check limit
    if len(models) > 5:
        raise ValueError("Select up to 5 models for comparison")


from proteingympy.make_zeroshot_dms_benchmarks import get_zero_shot_metrics

def benchmark_models(
    metric: Optional[str] = None,
    models: Optional[List[str]] = None,
    metric_tables: Optional[Dict[str, pd.DataFrame]] = None,
) -> Any:
    """
    Plot model benchmark scores similar to the R implementation.

    Parameters
    ----------
    metric : str, optional
        One of "AUC", "MCC", "NDCG", "Spearman", "Top_recall".
        If not provided, defaults to "Spearman" with a printed message.
    models : list of str, optional
        A list (up to 5) of model names to compare. If None, must provide at least one.
    metric_tables : dict, optional
        A dict mapping metric name -> pandas.DataFrame. If not provided, the
        function will call get_zero_shot_metrics().

    Returns
    -------
    matplotlib.figure.Figure
        The figure containing the plot.
    """

    # Default metric
    if metric is None:
        print("No metric specified. Using default Spearman correlation")
        metric = "Spearman"
    else:
        check_metric_argument(metric)

    # Ensure models argument is provided
    if models is None:
        raise ValueError("Select at least one model from available_models()")

    # Normalize models input
    if isinstance(models, str):
        models = [models]
    models_normalized = [m.strip().replace(" ", "_").lower() for m in models]

    # Check models are valid
    check_model_argument(models_normalized)

    # Load metric tables
    if metric_tables is None:
        metric_tables = get_zero_shot_metrics()

    if metric not in metric_tables:
        raise KeyError(
            f"Metric '{metric}' not found in metric_tables. Available metrics: {list(metric_tables.keys())}"
        )

    selected_table = metric_tables[metric]

    # Ensure DataFrame
    if not isinstance(selected_table, pd.DataFrame):
        raise ValueError("Each metric in metric_tables must be a pandas.DataFrame")

    # Normalize DataFrame columns for matching
    selected_table.columns = (
        selected_table.columns
        .str.strip()
        .str.replace(r"\s+", "_", regex=True)
        .str.lower()
    )

    # Check missing columns
    missing_cols = [m for m in models_normalized if m not in selected_table.columns]
    if missing_cols:
        raise KeyError(f"Selected model(s) not present in metric table: {missing_cols}")

    # Select columns
    selected_table = selected_table.loc[:, models_normalized]

    # If Spearman, take absolute value
    if metric == "Spearman":
        res = selected_table.abs()
    else:
        res = selected_table.copy()

    # Long form
    res_long = res.reset_index(drop=True).melt(var_name="model", value_name="score")

    # Compute mean per model and order descending
    model_means = res_long.groupby("model", observed=True)["score"].mean()
    ordered_models = list(model_means.sort_values(ascending=False).index)

    # Map normalized model names back to original names
    normalized_to_original = {
        m.strip().replace(" ", "_").lower(): m for m in available_models()
    }
    res_long["model"] = res_long["model"].map(normalized_to_original)
    ordered_models_original = [normalized_to_original[m] for m in ordered_models]

    # Set categorical order
    res_long["model"] = pd.Categorical(
        res_long["model"], categories=ordered_models_original, ordered=True
    )

    # Plotting
    sns.set(style="whitegrid")
    n_models = len(ordered_models_original)
    fig_width = max(8, n_models * 1.2)
    fig, ax = plt.subplots(figsize=(fig_width, 6))

    ggplot2_colors = ["#F8766D", "#E6D922", "#10C876", "#029FD3", "#C77CFF"]
    palette = ggplot2_colors[:n_models]

    sns.violinplot(
        x="model",
        y="score",
        data=res_long,
        order=ordered_models_original,
        ax=ax,
        inner=None,
        cut=0,
        palette=palette,
        linewidth=0.8,
        hue="model",
        legend=False
    )

    sns.boxplot(
        x="model",
        y="score",
        data=res_long,
        order=ordered_models_original,
        ax=ax,
        width=0.2,
        linewidth=1,
        showcaps=True,
        boxprops={"zorder": 2, "facecolor": "white"},
        showfliers=False,
        orient="vertical"
    )

    sns.stripplot(
        x="model",
        y="score",
        data=res_long,
        order=ordered_models_original,
        ax=ax,
        color="k",
        size=3,
        jitter=0.15,
        alpha=0.5,
        zorder=3
    )

    ax.set_xlabel("")
    ax.set_ylabel(f"{metric} score", fontsize=15)
    ax.tick_params(axis="x", labelrotation=15, labelsize=12)
    ax.tick_params(axis="y", labelsize=12)

    # Remove legend
    ax.get_legend() and ax.get_legend().remove()
    plt.tight_layout()

    return fig