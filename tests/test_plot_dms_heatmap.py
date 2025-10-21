"""
Tests for plot_dms_heatmap module.
"""

import pytest
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from proteingympy.plot_dms_heatmap import (
    plot_dms_heatmap,
    _filter_by_pos,
    _filter_exact_coord,
    _make_colormap_dms
)


@pytest.fixture
def sample_dms_data():
    """Create sample DMS data for testing."""
    # Create a simple dataset with mutants
    data = []
    positions = [1, 2, 3, 4, 5, 10, 11, 12]
    amino_acids = ['A', 'C', 'D', 'E', 'F', 'G']
    
    for pos in positions:
        for aa in amino_acids:
            mutant = f"A{pos}{aa}"  # Reference is always A
            score = np.random.randn()  # Random DMS score
            data.append({
                'mutant': mutant,
                'DMS_score': score,
                'UniProt_id': 'TEST_PROTEIN',
                'DMS_id': 'TEST_ASSAY'
            })
    
    df = pd.DataFrame(data)
    
    return {'TEST_ASSAY': df}


def test_filter_by_pos():
    """Test position filtering function."""
    df = pd.DataFrame({
        'pos': [1, 2, 3, 4, 5, 10, 11, 12],
        'value': range(8)
    })
    
    # Test with both start and end
    filtered = _filter_by_pos(df, start_pos=3, end_pos=10)
    assert len(filtered) == 5
    assert filtered['pos'].min() == 3
    assert filtered['pos'].max() == 10
    
    # Test with only start
    filtered = _filter_by_pos(df, start_pos=5)
    assert filtered['pos'].min() == 5
    assert filtered['pos'].max() == 12
    
    # Test with only end
    filtered = _filter_by_pos(df, end_pos=5)
    assert filtered['pos'].min() == 1
    assert filtered['pos'].max() == 5


def test_filter_by_pos_errors():
    """Test error handling in position filtering."""
    df = pd.DataFrame({
        'pos': [1, 2, 3, 4, 5],
        'value': range(5)
    })
    
    # Test with start_pos outside range
    with pytest.raises(ValueError, match="start_pos.*outside the assay range"):
        _filter_by_pos(df, start_pos=10)
    
    # Test with end_pos outside range
    with pytest.raises(ValueError, match="end_pos.*outside the assay range"):
        _filter_by_pos(df, end_pos=0)
    
    # Test missing pos column
    df_no_pos = pd.DataFrame({'value': range(5)})
    with pytest.raises(ValueError, match="must contain a 'pos' column"):
        _filter_by_pos(df_no_pos)


def test_filter_exact_coord():
    """Test exact coordinate filtering."""
    df = pd.DataFrame({
        'pos': [1, 3, 5],  # Missing positions 2 and 4
        'value': ['a', 'b', 'c']
    })
    
    # Test with exact_coord=True
    result = _filter_exact_coord(df, start_pos=1, end_pos=5, exact_coord=True)
    assert len(result) == 5
    assert list(result['pos']) == [1, 2, 3, 4, 5]
    assert pd.isna(result[result['pos'] == 2]['value'].iloc[0])
    assert pd.isna(result[result['pos'] == 4]['value'].iloc[0])
    
    # Test with exact_coord=False
    result = _filter_exact_coord(df, exact_coord=False)
    assert len(result) == 3
    assert list(result['pos']) == [1, 3, 5]


def test_make_colormap_dms():
    """Test colormap creation."""
    mat = np.array([[-1, 0, 1], [-2, 0, 2]])
    
    # Test default colormap
    cmap, vmin, vmax = _make_colormap_dms(mat, color_scheme="default")
    assert vmin == -2
    assert vmax == 2
    assert cmap is not None
    
    # Test EVE colormap
    cmap, vmin, vmax = _make_colormap_dms(mat, color_scheme="EVE")
    assert vmin == -2
    assert vmax == 2
    assert cmap is not None


def test_plot_dms_heatmap_basic(sample_dms_data):
    """Test basic heatmap plotting."""
    fig, ax = plot_dms_heatmap(
        assay_name='TEST_ASSAY',
        dms_data=sample_dms_data,
        start_pos=1,
        end_pos=5,
        figsize=(8, 6)
    )
    
    assert fig is not None
    assert ax is not None
    assert isinstance(fig, plt.Figure)
    
    plt.close(fig)


def test_plot_dms_heatmap_eve_colors(sample_dms_data):
    """Test heatmap with EVE color scheme."""
    fig, ax = plot_dms_heatmap(
        assay_name='TEST_ASSAY',
        dms_data=sample_dms_data,
        start_pos=1,
        end_pos=5,
        color_scheme='EVE',
        figsize=(8, 6)
    )
    
    assert fig is not None
    plt.close(fig)


def test_plot_dms_heatmap_exact_coord(sample_dms_data):
    """Test heatmap with exact coordinates."""
    fig, ax = plot_dms_heatmap(
        assay_name='TEST_ASSAY',
        dms_data=sample_dms_data,
        start_pos=1,
        end_pos=12,
        exact_coord=True,
        figsize=(10, 6)
    )
    
    assert fig is not None
    plt.close(fig)


def test_plot_dms_heatmap_invalid_assay(sample_dms_data):
    """Test error handling for invalid assay name."""
    with pytest.raises(ValueError, match="not found in dms_data"):
        plot_dms_heatmap(
            assay_name='NONEXISTENT_ASSAY',
            dms_data=sample_dms_data
        )


def test_plot_dms_heatmap_only_multiple_sites():
    """Test error handling when assay has only multiple amino acid sites."""
    # Create data with only multi-site mutants
    data = {
        'MULTI_SITE_ASSAY': pd.DataFrame({
            'mutant': ['A1P:D2N', 'A1G:D2E'],  # All contain ':'
            'DMS_score': [0.5, -0.5],
            'UniProt_id': 'TEST',
            'DMS_id': 'MULTI_SITE_ASSAY'
        })
    }
    
    with pytest.raises(ValueError, match="contains only multiple amino acid sites"):
        plot_dms_heatmap(
            assay_name='MULTI_SITE_ASSAY',
            dms_data=data
        )


def test_plot_dms_heatmap_clustering(sample_dms_data):
    """Test clustering functionality."""
    try:
        import scipy
        
        # Test row clustering
        fig, ax = plot_dms_heatmap(
            assay_name='TEST_ASSAY',
            dms_data=sample_dms_data,
            start_pos=1,
            end_pos=5,
            cluster_rows=True,
            figsize=(8, 6)
        )
        assert fig is not None
        plt.close(fig)
        
        # Test column clustering
        fig, ax = plot_dms_heatmap(
            assay_name='TEST_ASSAY',
            dms_data=sample_dms_data,
            start_pos=1,
            end_pos=5,
            cluster_columns=True,
            figsize=(8, 6)
        )
        assert fig is not None
        plt.close(fig)
        
    except ImportError:
        pytest.skip("scipy not installed")


def test_plot_dms_heatmap_clustering_with_na():
    """Test that clustering with NA values raises appropriate error."""
    data = {
        'TEST': pd.DataFrame({
            'mutant': ['A1P', 'A2C'],
            'DMS_score': [0.5, -0.5],
            'UniProt_id': 'TEST',
            'DMS_id': 'TEST'
        })
    }
    
    try:
        import scipy
        
        # This should raise an error when exact_coord=True creates NaN values
        with pytest.raises(ValueError, match="missing values, preventing clustering"):
            plot_dms_heatmap(
                assay_name='TEST',
                dms_data=data,
                start_pos=1,
                end_pos=10,  # Will create NaN for positions 3-10
                exact_coord=True,
                cluster_columns=True
            )
    except ImportError:
        pytest.skip("scipy not installed")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
