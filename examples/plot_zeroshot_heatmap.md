# plot_zeroshot_heatmap

## Overview

`plot_zeroshot_heatmap()` creates heatmap visualizations of zero-shot model predictions for amino acid substitutions along a protein sequence. This function is the Python equivalent of the R function `plot_zeroshot_heatmap()` from ProteinGymR.

## Function Signature

```python
def plot_zeroshot_heatmap(
    assay_name: str,
    model: str,
    model_data: Optional[Dict[str, pd.DataFrame]] = None,
    start_pos: Optional[int] = None,
    end_pos: Optional[int] = None,
    exact_coord: bool = False,
    cluster_rows: bool = False,
    cluster_columns: bool = False,
    color_scheme: str = "default",
    figsize: Tuple[float, float] = (12, 8),
    **kwargs
) -> Tuple[Figure, Axes]
```

## Parameters

- **assay_name** (`str`): Valid assay name that exists as a key in the `model_data` dictionary. To see available assays, check the keys of your loaded data.

- **model** (`str`): Name of the zero-shot model to plot. Must be a column in the assay dataframe. Common models include:
  - ESM1v
  - EVE
  - AlphaMissense
  - ProteinNPT
  - Kermut
  - And 70+ other zero-shot models

- **model_data** (`Dict[str, pd.DataFrame]`, optional): Dictionary of zero-shot assays where keys are assay names and values are DataFrames containing model predictions. If `None`, the function attempts to load data using `get_zero_shot_substitution_data()`.

- **start_pos** (`int`, optional): First amino acid position to plot. Defaults to the first available position in the assay.

- **end_pos** (`int`, optional): Last amino acid position to plot. Defaults to the last available position in the assay.

- **exact_coord** (`bool`, default=`False`): If `True`, plots precise start_pos and end_pos coordinates, filling missing positions with NaN (shown as grey in the heatmap). If `False`, only plots positions with available data.

- **cluster_rows** (`bool`, default=`False`): If `True`, performs hierarchical clustering on amino acid rows. Requires scipy.

- **cluster_columns** (`bool`, default=`False`): If `True`, performs hierarchical clustering on position columns. Requires scipy.

- **color_scheme** (`str`, default=`"default"`): Color scheme for the heatmap:
  - `"default"`: Parula-like colormap (blue → cyan → yellow)
  - `"EVE"`: popEVE portal color scheme (black → purple → cyan → yellow)

- **figsize** (`Tuple[float, float]`, default=`(12, 8)`): Figure size as (width, height) in inches.

- **kwargs**: Additional keyword arguments passed to `seaborn.heatmap()`.

## Returns

- **Tuple[Figure, Axes]**: A tuple containing the matplotlib Figure and Axes objects.

## Data Format

The input data should be structured as a dictionary where:
- Keys are assay names (strings)
- Values are pandas DataFrames with at least these columns:
  - `mutant`: Mutation identifier (e.g., "A1P" means position 1 mutated from A to P)
  - Model columns: One or more columns containing model prediction scores

Example DataFrame structure:
```
   mutant    ESM1v      EVE    AlphaMissense
0    A1C    -2.34     1.23         0.45
1    A1D     0.12    -0.56         0.67
2    A1E     1.45     2.11         0.23
...
```

## Examples

### Basic Usage

```python
from proteingympy.plot_zeroshot_heatmap import plot_zeroshot_heatmap
from proteingympy import get_zero_shot_substitution_data
import matplotlib.pyplot as plt

# Load zero-shot model data
model_data = get_zero_shot_substitution_data()

# Create a basic heatmap
fig, ax = plot_zeroshot_heatmap(
    assay_name="ACE2_HUMAN_Chan_2020",
    model="ESM1v",
    model_data=model_data,
    start_pos=1,
    end_pos=100
)

plt.show()
```

### Using EVE Color Scheme

```python
# Create heatmap with popEVE portal colors
fig, ax = plot_zeroshot_heatmap(
    assay_name="ACE2_HUMAN_Chan_2020",
    model="EVE",
    model_data=model_data,
    start_pos=50,
    end_pos=150,
    color_scheme="EVE"
)

plt.savefig("eve_heatmap.png", dpi=300, bbox_inches='tight')
```

### With Clustering

```python
# Cluster amino acids by similar prediction patterns
fig, ax = plot_zeroshot_heatmap(
    assay_name="SHOC2_HUMAN_Newby_2022",
    model="AlphaMissense",
    model_data=model_data,
    start_pos=100,
    end_pos=200,
    cluster_rows=True,
    figsize=(14, 10)
)

plt.show()
```

### Exact Coordinates (Show Gaps)

```python
# Show all positions including those without data (as grey/NaN)
fig, ax = plot_zeroshot_heatmap(
    assay_name="P53_HUMAN_Giacomelli_2018",
    model="ProteinNPT",
    model_data=model_data,
    start_pos=100,
    end_pos=200,
    exact_coord=True
)

plt.show()
```

## Interpretation

The heatmap displays:
- **X-axis**: Amino acid positions along the protein sequence
- **Y-axis**: Possible amino acid residues, ordered by physiochemical properties (DEKRHNQSTPGAVILMCFYW)
- **Top annotation**: Reference amino acid at each position
- **Color intensity**: Model prediction score
  - Default scheme: Blue = low scores, Yellow = high scores
  - EVE scheme: Black/Purple = low scores, Cyan/Yellow = high scores
- **Grey cells**: Missing data (NaN values)

Higher model scores typically indicate predicted higher fitness or less pathogenic variants, while lower scores suggest deleterious effects (exact interpretation depends on the specific model).

## Physiochemical Ordering

Amino acids are ordered by default according to their physiochemical properties:
1. **DE** - Acidic (negatively charged)
2. **KRH** - Basic (positively charged)
3. **NQ** - Polar uncharged (amides)
4. **ST** - Polar uncharged (hydroxyls)
5. **PGAVIL** - Small/hydrophobic
6. **MC** - Sulfur-containing
7. **FYW** - Aromatic

This ordering helps identify patterns in how different types of amino acid substitutions affect the protein.

## Error Handling

The function raises:
- `ValueError`: If the assay_name is not found, model column doesn't exist, or assay contains only multiple amino acid sites
- `ImportError`: If required packages (scipy for clustering) are missing

## Notes

1. The function automatically filters out multiple amino acid substitutions (those containing ':' in the mutant identifier).
2. When `exact_coord=True` and `cluster_columns=True` are both used, missing values can prevent clustering. Set `exact_coord=False` to avoid this.
3. The colorbar shows three tick marks: minimum, midpoint, and maximum values with 2 decimal places.
4. Clustering uses hierarchical clustering with average linkage and Euclidean distance (after replacing NaN with 0).

## Comparison with R Version

This Python implementation maintains feature parity with the R version from ProteinGymR:

| Feature | R Version | Python Version |
|---------|-----------|----------------|
| Basic heatmap | ✓ (ComplexHeatmap) | ✓ (seaborn) |
| Position filtering | ✓ | ✓ |
| Exact coordinates | ✓ | ✓ |
| Row/column clustering | ✓ | ✓ (requires scipy) |
| EVE color scheme | ✓ | ✓ |
| Default color scheme | ✓ (parula) | ✓ (parula-like) |
| Reference AA annotation | ✓ | ✓ |
| Physiochemical ordering | ✓ | ✓ |

Key differences:
- R uses ComplexHeatmap package; Python uses matplotlib/seaborn
- R integrates with Bioconductor ExperimentHub; Python uses local/remote data loading
- Python returns Figure and Axes objects for further customization

## See Also

- `plot_dms_heatmap()`: Similar function for plotting experimental DMS scores (rather than model predictions)
- `plot_structure()`: Visualize predictions on 3D protein structures
- `model_corr_plot()`: Compare predictions between two models
- `get_zero_shot_substitution_data()`: Load zero-shot model prediction data
- `available_zero_shot_models()`: List all available zero-shot models
