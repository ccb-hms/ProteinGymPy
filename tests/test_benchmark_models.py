#!/usr/bin/env python3
"""
Tests for benchmark_models().

"""
import pytest
import sys
import os
import pandas as pd
import numpy as np

# Add src to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))


from proteingympy.benchmark_models import (
    check_metric_argument,
    check_model_argument,
    benchmark_models,
    available_models,
)
from proteingympy.make_zeroshot_dms_benchmarks import get_zero_shot_metrics


# ---------------------------------------------------------
# check_metric_argument() TESTS
# ---------------------------------------------------------

def test_check_metric_argument_invalid():
    """Invalid metric should raise an error"""
    with pytest.raises(ValueError, match="Invalid metric specified: Pearson"):
        check_metric_argument("Pearson")


def test_check_metric_argument_multiple():
    """Passing a list should raise TypeError because only single string allowed"""
    with pytest.raises(TypeError, match="Metric must be a single string"):
        check_metric_argument(["AUC", "MCC"])


# ---------------------------------------------------------
# check_model_argument() TESTS
# ---------------------------------------------------------

def test_check_model_argument_invalid():
    with pytest.raises(ValueError, match="Invalid model\\(s\\) specified: \\['Wrong_model'\\]"):
        check_model_argument("Wrong_model")


def test_check_model_argument_too_many():
    too_many = ["Site_Independent", "EVE_single", "GEMME", "VESPA", "PoET", "CARP_640M"]
    with pytest.raises(ValueError, match="Select up to 5 models for comparison"):
        check_model_argument(too_many)


# ---------------------------------------------------------
# benchmark_models() TESTS
# ---------------------------------------------------------
import pytest
import pandas as pd
from proteingympy.benchmark_models import benchmark_models, available_models

# 1. Helper: Create a fake metric table for the tests
@pytest.fixture
def example_metric_tables():
    # Let's use a subset of available models for a small table
    model_names = available_models()[:4]
    data = {
        model: [0.8 + i*0.01 for i in range(5)] for model in model_names
    }
    # Slight difference for each model
    metric_tables = {
        "Spearman": pd.DataFrame(data),
        "AUC": pd.DataFrame(data)
    }
    return metric_tables

# 2. Test: Successful plot returns a matplotlib Figure
def test_benchmark_models_returns_figure(example_metric_tables):
    import matplotlib
    models = available_models()[:3]
    fig = benchmark_models(
        metric="Spearman",
        models=models,
        metric_tables=example_metric_tables
    )
    assert isinstance(fig, matplotlib.figure.Figure)
    # Optionally: check number of axes, or other simple plot properties

# 3. Test: Raises ValueError for >5 models
def test_benchmark_models_raises_too_many_models(example_metric_tables):
    # Use 6 models, which should raise
    models = available_models()[:6]
    with pytest.raises(ValueError, match="up to 5 models"):
        benchmark_models(
            metric="Spearman",
            models=models,
            metric_tables=example_metric_tables
        )

# 4. Test: Raises KeyError for invalid metric
def test_benchmark_models_raises_invalid_metric(example_metric_tables):
    models = available_models()[:3]
    with pytest.raises(KeyError, match="not found in metric_tables"):
        benchmark_models(
            metric="Top_recall",  # Not present in metric_tables
            models=models,
            metric_tables=example_metric_tables
        )

# 5. Test: Raises KeyError for invalid model column
def test_benchmark_models_raises_invalid_model_column(example_metric_tables):
    models = available_models()[:2] + ["NOT_A_REAL_MODEL"]
    with pytest.raises(ValueError, match="Invalid model.*specified"):
        benchmark_models(
            metric="Spearman",
            models=models,
            metric_tables=example_metric_tables
        )