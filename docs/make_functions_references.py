import os

functions = [
    'get_dms_substitution_data',
    'get_dms_metadata',
    'get_alphamissense_proteingym_data',
    'get_supervised_substitution_data',
    'available_supervised_models',
    'get_zero_shot_substitution_data',
    'available_zero_shot_models',
    'get_zero_shot_metrics',
    'create_complete_metadata_table',
    'benchmark_models',
    'dms_corr_plot',
    'model_corr_plot',
    'plot_dms_heatmap',
    'plot_structure'
]

os.makedirs('docs/reference', exist_ok=True)

for func in functions:
    with open(f'docs/reference/{func}.md', 'w') as f:
        f.write(f'# {func}\n\n::: proteingympy.{func}\n')