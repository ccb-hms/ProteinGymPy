"""
Example usage of plot_zeroshot_heatmap function

This script demonstrates how to create heatmaps of zero-shot model predictions
for amino acid substitutions along a protein sequence.
"""

import sys
import os

# Add src to path for local development
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from proteingympy.plot_zeroshot_heatmap import plot_zeroshot_heatmap
import matplotlib.pyplot as plt
import pandas as pd


def create_example_data():
    """
    Create example zero-shot model data for demonstration.
    
    This simulates the structure of data that would be loaded from
    get_zero_shot_substitution_data() function.
    """
    import numpy as np
    
    # Create synthetic data for a protein with positions 1-50
    AMINO_ACIDS = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                   'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    mutations = [
        {
            'mutant': f'{ref_aa}{pos}{alt_aa}',
            'ESM1v': np.random.randn() * 2,
            'EVE': np.random.randn() * 2,
            'AlphaMissense': np.random.uniform(0, 1)
        }
        for pos in range(1, 51)
        for ref_aa in AMINO_ACIDS
        for alt_aa in AMINO_ACIDS
        if ref_aa != alt_aa
    ]
    
    df = pd.DataFrame(mutations)
    
    # Return as a dictionary with assay name as key
    return {'EXAMPLE_PROTEIN_Test_2024': df}


def example_basic_heatmap():
    """Example 1: Basic heatmap with default settings"""
    print("Example 1: Basic heatmap")
    
    # Load or create example data
    model_data = create_example_data()
    
    # Create heatmap for ESM1v model
    fig, ax = plot_zeroshot_heatmap(
        assay_name="EXAMPLE_PROTEIN_Test_2024",
        model="ESM1v",
        model_data=model_data,
        start_pos=1,
        end_pos=30
    )
    
    plt.savefig('example_zeroshot_heatmap_basic.png', dpi=150, bbox_inches='tight')
    print("Saved: example_zeroshot_heatmap_basic.png")
    plt.close()


def example_eve_colorscheme():
    """Example 2: Using EVE color scheme"""
    print("\nExample 2: EVE color scheme")
    
    model_data = create_example_data()
    
    # Create heatmap with EVE colors
    fig, ax = plot_zeroshot_heatmap(
        assay_name="EXAMPLE_PROTEIN_Test_2024",
        model="EVE",
        model_data=model_data,
        start_pos=10,
        end_pos=40,
        color_scheme="EVE"
    )
    
    plt.savefig('example_zeroshot_heatmap_eve.png', dpi=150, bbox_inches='tight')
    print("Saved: example_zeroshot_heatmap_eve.png")
    plt.close()


def example_with_clustering():
    """Example 3: With row clustering"""
    print("\nExample 3: With clustering")
    
    model_data = create_example_data()
    
    # Create heatmap with clustering
    fig, ax = plot_zeroshot_heatmap(
        assay_name="EXAMPLE_PROTEIN_Test_2024",
        model="AlphaMissense",
        model_data=model_data,
        start_pos=20,
        end_pos=50,
        cluster_rows=True,
        figsize=(14, 8)
    )
    
    plt.savefig('example_zeroshot_heatmap_clustered.png', dpi=150, bbox_inches='tight')
    print("Saved: example_zeroshot_heatmap_clustered.png")
    plt.close()


def example_exact_coordinates():
    """Example 4: Using exact coordinates to show gaps"""
    print("\nExample 4: Exact coordinates with gaps")
    
    model_data = create_example_data()
    
    # Remove some positions to create gaps
    df = model_data['EXAMPLE_PROTEIN_Test_2024']
    # Keep only odd positions between 10 and 30
    df_filtered = df[df['mutant'].str.extract(r'(\d+)')[0].astype(int) % 2 == 1]
    df_filtered = df_filtered[
        (df_filtered['mutant'].str.extract(r'(\d+)')[0].astype(int) >= 10) &
        (df_filtered['mutant'].str.extract(r'(\d+)')[0].astype(int) <= 30)
    ]
    model_data_gaps = {'EXAMPLE_PROTEIN_Test_2024': df_filtered}
    
    # Create heatmap with exact coordinates (will show NaN for missing positions)
    fig, ax = plot_zeroshot_heatmap(
        assay_name="EXAMPLE_PROTEIN_Test_2024",
        model="ESM1v",
        model_data=model_data_gaps,
        start_pos=10,
        end_pos=30,
        exact_coord=True
    )
    
    plt.savefig('example_zeroshot_heatmap_gaps.png', dpi=150, bbox_inches='tight')
    print("Saved: example_zeroshot_heatmap_gaps.png")
    plt.close()


if __name__ == "__main__":
    print("Generating zero-shot heatmap examples...\n")
    print("=" * 60)
    
    # Run examples
    example_basic_heatmap()
    example_eve_colorscheme()
    
    try:
        example_with_clustering()
    except ImportError as e:
        print(f"\nSkipping clustering example (requires scipy): {e}")
    
    example_exact_coordinates()
    
    print("\n" + "=" * 60)
    print("\nAll examples completed!")
    print("\nNote: These examples use synthetic data for demonstration.")
    print("For real data, use get_zero_shot_substitution_data() from proteingympy")
