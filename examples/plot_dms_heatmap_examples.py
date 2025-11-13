"""
Examples of using plot_dms_heatmap() to visualize DMS scores.

This script demonstrates various ways to create DMS heatmap visualizations
using the proteingympy package.
"""

import matplotlib.pyplot as plt
from proteingympy import get_dms_substitution_data, plot_dms_heatmap


def example_basic_heatmap():
    """Basic heatmap with default settings."""
    print("Example 1: Basic heatmap")
    print("-" * 50)
    
    # Load DMS data
    dms_data = get_dms_substitution_data()
    
    # Create basic heatmap
    fig, ax = plot_dms_heatmap(
        assay_name="A0A192B1T2_9HIV1_Haddox_2018",
        dms_data=dms_data,
        start_pos=10,
        end_pos=80
    )
    
    plt.savefig("heatmap_basic.png", dpi=300, bbox_inches='tight')
    print("Saved: heatmap_basic.png\n")
    plt.close()


def example_eve_colorscheme():
    """Heatmap with EVE color scheme and exact coordinates."""
    print("Example 2: EVE color scheme with exact coordinates")
    print("-" * 50)
    
    dms_data = get_dms_substitution_data()
    
    fig, ax = plot_dms_heatmap(
        assay_name="A0A192B1T2_9HIV1_Haddox_2018",
        dms_data=dms_data,
        start_pos=10,
        end_pos=80,
        exact_coord=True,
        color_scheme="EVE",
        figsize=(14, 8)
    )
    
    plt.savefig("heatmap_eve.png", dpi=300, bbox_inches='tight')
    print("Saved: heatmap_eve.png\n")
    plt.close()


def example_clustered_heatmap():
    """Heatmap with row clustering."""
    print("Example 3: Clustered heatmap")
    print("-" * 50)
    
    dms_data = get_dms_substitution_data()
    
    fig, ax = plot_dms_heatmap(
        assay_name="A0A192B1T2_9HIV1_Haddox_2018",
        dms_data=dms_data,
        start_pos=50,
        end_pos=100,
        cluster_rows=True,
        figsize=(12, 8)
    )
    
    plt.savefig("heatmap_clustered.png", dpi=300, bbox_inches='tight')
    print("Saved: heatmap_clustered.png\n")
    plt.close()


def example_custom_styling():
    """Heatmap with custom styling options."""
    print("Example 4: Custom styling")
    print("-" * 50)
    
    dms_data = get_dms_substitution_data()
    
    fig, ax = plot_dms_heatmap(
        assay_name="A0A192B1T2_9HIV1_Haddox_2018",
        dms_data=dms_data,
        start_pos=20,
        end_pos=60,
        color_scheme="default",
        figsize=(15, 6),
        linewidths=0.5,
        linecolor='lightgray',
        cbar_kws={'label': 'Fitness Score', 'shrink': 0.8}
    )
    
    ax.set_title("DMS Heatmap: A0A192B1T2_9HIV1_Haddox_2018 (positions 20-60)", 
                 fontsize=14, pad=20)
    
    plt.savefig("heatmap_custom.png", dpi=300, bbox_inches='tight')
    print("Saved: heatmap_custom.png\n")
    plt.close()


def example_multiple_assays():
    """Create heatmaps for multiple assays side by side."""
    print("Example 5: Multiple assays comparison")
    print("-" * 50)
    
    dms_data = get_dms_substitution_data()
    
    # Get list of available assays (use first 2 for demo)
    assay_names = list(dms_data.keys())[:2]
    
    fig, axes = plt.subplots(1, 2, figsize=(20, 8))
    
    for i, assay in enumerate(assay_names):
        print(f"Processing assay {i+1}: {assay}")
        
        # Get position range for this assay
        df = dms_data[assay]
        df = df[~df['mutant'].str.contains(':', na=False)]
        if len(df) > 0:
            df['pos'] = df['mutant'].str.extract(r'(\d+)')[0].astype(int)
            start = int(df['pos'].min())
            end = min(start + 50, int(df['pos'].max()))
            
            # Create subplot
            plt.sca(axes[i])
            _, _ = plot_dms_heatmap(
                assay_name=assay,
                dms_data=dms_data,
                start_pos=start,
                end_pos=end,
                figsize=(10, 8)  # Will be ignored due to existing axes
            )
            axes[i].set_title(assay, fontsize=10)
    
    plt.tight_layout()
    plt.savefig("heatmap_multiple.png", dpi=300, bbox_inches='tight')
    print("Saved: heatmap_multiple.png\n")
    plt.close()


def main():
    """Run all examples."""
    print("=" * 50)
    print("DMS Heatmap Visualization Examples")
    print("=" * 50)
    print()
    
    try:
        example_basic_heatmap()
        example_eve_colorscheme()
        example_clustered_heatmap()
        example_custom_styling()
        # example_multiple_assays()  # Commented out as it takes longer
        
        print("=" * 50)
        print("All examples completed successfully!")
        print("=" * 50)
        
    except Exception as e:
        print(f"Error running examples: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
