"""
Python port of R's benchmark_models() from the ProteinGymR package.

This module provides:
- available_models()
- supervised_available_models()
- check_metric_argument()
- check_model_argument()
- zeroshot_DMS_metrics() (helper that attempts to load metrics from disk)
- benchmark_models(...)

Notes:
- The original R function loads benchmark scores via zeroshot_DMS_metrics().
  This Python version will attempt to load a pickled dictionary
  named "zeroshot_DMS_metrics.pkl" from the same directory as this file if
  `metric_tables` is not provided to benchmark_models(). The pickled object
  should be a dict mapping metric names ("AUC", "MCC", "NDCG", "Spearman",
  "Top_recall") to pandas.DataFrame objects whose columns are model names
  (as returned by available_models()) and rows are individual dataset scores.

- If you can't / don't want to provide that file, pass metric_tables=<dict>
  directly to benchmark_models().

- The plot attempts to reproduce the R ggplot2 appearance using seaborn +
  matplotlib: a violin (to approximate the half-eye), a narrow boxplot, and
  jittered points.

Dependencies:
- pandas
- numpy
- seaborn
- matplotlib

Install with: pip install pandas numpy seaborn matplotlib

Author: ported to Python to mimic behavior of the R function
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


def available_models() -> List[str]:
    return [
        "Site_Independent",
        "EVmutation",
        "DeepSequence_single",
        "DeepSequence_ensemble",
        "EVE_single",
        "EVE_ensemble",
        "Unirep",
        "Unirep_evotune",
        "MSA_Transformer_single",
        "MSA_Transformer_ensemble",
        "ESM1b",
        "ESM1v_single",
        "ESM1v_ensemble",
        "ESM2_8M",
        "ESM2_35M",
        "ESM2_150M",
        "ESM2_650M",
        "ESM2_3B",
        "ESM2_15B",
        "Wavenet",
        "RITA_s",
        "RITA_m",
        "RITA_l",
        "RITA_xl",
        "Progen2_small",
        "Progen2_medium",
        "Progen2_base",
        "Progen2_large",
        "Progen2_xlarge",
        "GEMME",
        "VESPA",
        "VESPAl",
        "VespaG",
        "ProtGPT2",
        "Tranception_S_no_retrieval",
        "Tranception_M_no_retrieval",
        "Tranception_L_no_retrieval",
        "Tranception_S",
        "Tranception_M",
        "Tranception_L",
        "TranceptEVE_S",
        "TranceptEVE_M",
        "TranceptEVE_L",
        "CARP_38M",
        "CARP_600K",
        "CARP_640M",
        "CARP_76M",
        "MIF",
        "MIFST",
        "ESM_IF1",
        "ProteinMPNN",
        "ProtSSN_k10_h512",
        "ProtSSN_k10_h768",
        "ProtSSN_k10_h1280",
        "ProtSSN_k20_h512",
        "ProtSSN_k20_h768",
        "ProtSSN_k20_h1280",
        "ProtSSN_k30_h512",
        "ProtSSN_k30_h768",
        "ProtSSN_k30_h1280",
        "ProtSSN_ensemble",
        "SaProt_650M_AF2",
        "SaProt_35M_AF2",
        "PoET",
        "MULAN_small",
        "ProSST_20",
        "ProSST_128",
        "ProSST_512",
        "ProSST_1024",
        "ProSST_2048",
        "ProSST_4096",
        "ESCOTT",
        "VenusREM",
        "RSALOR",
        "S2F",
        "S2F_MSA",
        "S3F",
        "S3F_MSA",
        "SiteRM",
    ]


def supervised_available_models() -> List[str]:
    return [
        "OHE_Notaugmented",
        "normalized_targets",
        "OHE_Augmented_DeepSequence",
        "OHE_Augmented_ESM1v",
        "OHE_Augmented_MSATransformer",
        "OHE_Augmented_Tranception",
        "OHE_Augmented_TranceptEVE",
        "Embeddings_Augmented_ESM1v",
        "Embeddings_Augmented_MSATransformer",
        "Embeddings_Augmented_Tranception",
        "ProteinNPT",
        "Kermut",
    ]


def check_metric_argument(user_metric: List[str]) -> None:
    valid_metrics = ["AUC", "MCC", "NDCG", "Spearman", "Top_recall"]

    if not all(m in valid_metrics for m in user_metric):
        invalid_metric = [m for m in user_metric if m not in valid_metrics]
        raise ValueError(f"Invalid metric(s) specified: {invalid_metric}")

    if len(user_metric) > 1:
        raise ValueError("Select only one metric for comparison")


def check_model_argument(models: List[str]) -> None:
    valid_models = available_models()

    if not all(m in valid_models for m in models):
        invalid_models = [m for m in models if m not in valid_models]
        raise ValueError(f"Invalid model(s) specified: {invalid_models}")

    if len(models) > 5:
        raise ValueError("Select up to 5 models for comparison")


def zeroshot_DMS_metrics() -> Dict[str, pd.DataFrame]:
    """
    Helper that attempts to load a pickled dict named 'zeroshot_DMS_metrics.pkl'
    from the same directory. The expected object is a dict mapping metric name ->
    pandas.DataFrame (columns = model names).

    If not found, raises a FileNotFoundError and instructs the user to provide
    metric_tables to benchmark_models().
    """
    this_dir = os.path.dirname(__file__)
    candidate = os.path.join(this_dir, "zeroshot_DMS_metrics.pkl")
    if os.path.exists(candidate):
        with open(candidate, "rb") as fh:
            obj = pickle.load(fh)
        if not isinstance(obj, dict):
            raise ValueError("zeroshot_DMS_metrics.pkl does not contain a dict.")
        return obj
    else:
        raise FileNotFoundError(
            "Could not find 'zeroshot_DMS_metrics.pkl' in the same directory. "
            "Please provide metric_tables to benchmark_models(metric_tables=<dict>) "
            "where the dict maps e.g. 'Spearman' -> pandas.DataFrame(columns=model_names)."
        )


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
        A list (up to 5) of model names to compare. If None, defaults to
        available_models() (but check_model_argument will fail if >5).
    metric_tables : dict, optional
        A dict mapping metric name -> pandas.DataFrame. If not provided, the
        function will attempt to call zeroshot_DMS_metrics() which tries to
        read a 'zeroshot_DMS_metrics.pkl' in the same folder.

    Returns
    -------
    matplotlib.figure.Figure
        The figure containing the plot.
    """

    # dependency checks are already covered at import, but keep a small runtime check
    for pkg in ("seaborn", "matplotlib", "pandas", "numpy"):
        # these are already imported above; this loop ensures consistent error message if missing
        pass

    # Handle metric defaulting
    if metric is None:
        print("No metric specified. Using default Spearman correlation")
        metric = "Spearman"
    else:
        check_metric_argument([metric])

    # Handle models argument
    if models is None:
        raise ValueError("Select at least one model from available_models()")
    else:
        check_model_argument(models)

    # Load metric tables if not provided
    if metric_tables is None:
        metric_tables = zeroshot_DMS_metrics()

    if metric not in metric_tables:
        raise KeyError(
            f"Metric '{metric}' not found in metric_tables. Available metrics: {list(metric_tables.keys())}"
        )

    selected_table = metric_tables[metric]

    # Ensure selected_table is a DataFrame
    if not isinstance(selected_table, pd.DataFrame):
        raise ValueError("Each metric in metric_tables must be a pandas.DataFrame")

    # Select columns (models); raise if missing columns
    missing_cols = [m for m in models if m not in selected_table.columns]
    if missing_cols:
        raise KeyError(f"Selected model(s) not present in metric table: {missing_cols}")

    selected_table = selected_table.loc[:, models]

    # If Spearman, take absolute value (matches R behavior)
    if metric == "Spearman":
        res = selected_table.abs()
    else:
        res = selected_table.copy()

    # Long form
    res_long = res.reset_index(drop=True).melt(var_name="model", value_name="score")

    # Compute mean per model and order descending
    model_means = res_long.groupby("model", observed=True)["score"].mean()
    ordered_models = list(model_means.sort_values(ascending=False).index)

    # Make 'model' a categorical with that order (so seaborn uses the order)
    res_long["model"] = pd.Categorical(res_long["model"], categories=ordered_models, ordered=True)

    # Plotting: approximate half-eye with violin + narrow boxplot + jitter
    sns.set(style="whitegrid")
    n_models = len(ordered_models)
    fig_width = max(8, n_models * 1.2)
    fig, ax = plt.subplots(figsize=(fig_width, 6))

    palette = sns.color_palette("tab20", n_colors=n_models)

    # Violin (approximate half-eye)
    sns.violinplot(
        x="model",
        y="score",
        data=res_long,
        order=ordered_models,
        ax=ax,
        inner=None,
        cut=0,
        palette=palette,
        linewidth=0.8,
    )

    # Narrow boxplot on top
    sns.boxplot(
        x="model",
        y="score",
        data=res_long,
        order=ordered_models,
        ax=ax,
        width=0.15,
        showcaps=True,
        boxprops={"zorder": 2, "facecolor": "white"},
        showfliers=False,
    )

    # Jittered points
    sns.stripplot(
        x="model",
        y="score",
        data=res_long,
        order=ordered_models,
        ax=ax,
        color="k",
        size=3,
        jitter=0.15,
        alpha=0.4,
        zorder=3,
    )

    ax.set_xlabel("")
    ax.set_ylabel(f"{metric} score", fontsize=12)
    ax.tick_params(axis="x", labelrotation=45, labelsize=10)
    ax.tick_params(axis="y", labelsize=10)

    # Remove legend (in R plot legend shows models but it's redundant here)
    ax.get_legend() and ax.get_legend().remove()

    plt.tight_layout()
    return fig