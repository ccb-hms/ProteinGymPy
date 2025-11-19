from .data_import_funcs import get_dms_substitution_zip

# Import all new data pipeline functions
from .make_dms_substitutions import get_dms_substitution_data, get_dms_metadata
from .make_alphamissense_supplementary import get_alphamissense_proteingym_data, get_alphamissense_summary_stats
from .make_metadata import create_complete_metadata_table, save_metadata_csv
from .make_supervised_scores import get_supervised_scores_data, get_supervised_model_list
from .make_zero_shot_substitutions import get_zero_shot_scores_data, get_zero_shot_model_list
from .make_zeroshot_dms_benchmarks import get_zero_shot_benchmark_data, get_benchmark_summary_stats, get_top_models_by_metric
from .model_corr_plot import model_corr_plot, get_available_models
from .dms_corr_plot import dms_corr_plot
from .benchmark_models import benchmark_models
from .benchmark_models import available_models


# For convenience, expose the most commonly used functions
__all__ = [
    'get_dms_substitution_zip',  # Original function
    'get_dms_substitution_data',  # New comprehensive DMS data loader
    'get_dms_metadata',
    'get_alphamissense_proteingym_data',
    'get_alphamissense_summary_stats',
    'get_supervised_scores_data',
    'get_supervised_model_list',
    'get_zero_shot_scores_data',
    'get_zero_shot_model_list',
    'get_zero_shot_benchmark_data',
    'get_benchmark_summary_stats',
    'get_top_models_by_metric',
    'create_complete_metadata_table',
    'save_metadata_csv',
    'model_corr_plot',
    'get_available_models',
    'benchmark_models',
    'available_models',
    'main'
]


def main() -> None:
    print("Hello from proteingympy!")
    print("Available data pipeline functions:")
    print("  - get_dms_substitution_data(): Load 217 DMS substitution assays")
    print("  - get_alphamissense_proteingym_data(): Load AlphaMissense pathogenicity scores")
    print("  - get_supervised_scores_data(): Load supervised model predictions")
    print("  - get_zero_shot_scores_data(): Load zero-shot model predictions")
    print("  - get_zero_shot_benchmark_data(): Load benchmarking metrics")
    print("  - get_dms_metadata(): Load DMS assay metadata")
    print("  - create_complete_metadata_table(): Generate dataset metadata")
    print("  - model_corr_plot(): Compare prediction scores between models")
    print("  - get_available_models(): Get list of all available models")
    print("  - available_models(): Get list of all available zero-shot models for benchmarking")
    print("  - benchmark_models(): Benchmark zero-shot models on DMS data")
