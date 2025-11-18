#!/usr/bin/env python3
"""
Example script demonstrating dms_corr_plot functionality.

This script shows how to compare DMS experimental scores with model predictions
using Spearman correlation and visualization.
"""

import os
import sys

# Add src to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import matplotlib.pyplot as plt
from proteingympy.dms_corr_plot import (
    dms_corr_plot,
    get_available_models
)


def example_basic_usage():
    """
    Example 1: Basic usage with default model (AlphaMissense vs DMS)
    """
    print("=" * 70)
    print("Example 1: Basic Usage - AlphaMissense vs DMS scores")
    print("=" * 70)

    # Use a known UniProt ID from the example
    uniprot_id = "Q9NV35"

    print(f"\nComparing AlphaMissense predictions with DMS scores for {uniprot_id}...")

    try:
        fig = dms_corr_plot(
            uniprot_id=uniprot_id,
            model="AlphaMissense"
        )

        # Save figure
        output_path = f"dms_corr_{uniprot_id}_AlphaMissense.png"
        fig.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"\n✓ Figure saved to {output_path}")

        # Show interactively (comment out if running headless)
        # plt.show()

        plt.close(fig)

    except Exception as e:
        print(f"\n✗ Error: {e}")
        print("This is expected if the data is not yet downloaded.")


def example_zero_shot_model():
    """
    Example 2: Compare DMS scores with a zero-shot model
    """
    print("\n" + "=" * 70)
    print("Example 2: Zero-shot Model - GEMME vs DMS scores")
    print("=" * 70)

    # Use P04637 (TP53) as another example
    uniprot_id = "P04637"

    print(f"\nComparing GEMME predictions with DMS scores for {uniprot_id}...")

    try:
        fig = dms_corr_plot(
            uniprot_id=uniprot_id,
            model="GEMME",
            bins=40,  # Fewer bins for different visualization
            cmap='plasma'  # Different colormap
        )

        # Save figure
        output_path = f"dms_corr_{uniprot_id}_GEMME.png"
        fig.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"\n✓ Figure saved to {output_path}")

        plt.close(fig)

    except Exception as e:
        print(f"\n✗ Error: {e}")
        print("This is expected if the data is not yet downloaded.")


def example_large_language_model():
    """
    Example 3: Compare DMS scores with a large language model
    """
    print("\n" + "=" * 70)
    print("Example 3: Large Language Model - ESM vs DMS scores")
    print("=" * 70)

    uniprot_id = "Q9NV35"

    print(f"\nComparing ESM predictions with DMS scores for {uniprot_id}...")

    try:
        fig = dms_corr_plot(
            uniprot_id=uniprot_id,
            model="ESM_2_t48_15B_UR50D",
            figsize=(12, 10),  # Larger figure
            cmap='magma'
        )

        # Save figure
        output_path = f"dms_corr_{uniprot_id}_ESM.png"
        fig.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"\n✓ Figure saved to {output_path}")

        plt.close(fig)

    except Exception as e:
        print(f"\n✗ Error: {e}")
        print("This is expected if the data is not yet downloaded.")


def example_supervised_model():
    """
    Example 4: Compare DMS scores with a supervised model
    """
    print("\n" + "=" * 70)
    print("Example 4: Supervised Model - EVE vs DMS scores")
    print("=" * 70)

    uniprot_id = "Q9NV35"

    print(f"\nComparing EVE predictions with DMS scores for {uniprot_id}...")

    try:
        fig = dms_corr_plot(
            uniprot_id=uniprot_id,
            model="EVE",
            cmap='viridis'
        )

        # Save figure
        output_path = f"dms_corr_{uniprot_id}_EVE.png"
        fig.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"\n✓ Figure saved to {output_path}")

        plt.close(fig)

    except Exception as e:
        print(f"\n✗ Error: {e}")
        print("This is expected if the data is not yet downloaded.")


def show_available_models():
    """
    Show all available models for comparison
    """
    print("\n" + "=" * 70)
    print("Available Models in ProteinGym")
    print("=" * 70)

    models = get_available_models()

    print(f"\nZero-shot models ({len(models['zero_shot'])}):")
    for i, model in enumerate(models['zero_shot'][:10], 1):  # Show first 10
        print(f"  {i:2d}. {model}")
    if len(models['zero_shot']) > 10:
        print(f"  ... and {len(models['zero_shot']) - 10} more")

    print(f"\nSupervised models ({len(models['supervised'])}):")
    for i, model in enumerate(models['supervised'], 1):
        print(f"  {i:2d}. {model}")

    print(f"\nOther models ({len(models['other'])}):")
    for i, model in enumerate(models['other'], 1):
        print(f"  {i:2d}. {model}")

    print(f"\nTotal: {sum(len(v) for v in models.values())} models available")


def main():
    """
    Run all examples
    """
    print("ProteinGym DMS Correlation Plot Examples")
    print("=" * 70)
    print("\nThese examples demonstrate how to compare DMS experimental scores")
    print("with model predictions using Spearman correlation.")
    print("\nNote: Examples require downloaded data. First run will download")
    print("data from Zenodo, which may take some time.")
    print("=" * 70)

    # Show available models
    show_available_models()

    # Run examples
    example_basic_usage()
    example_zero_shot_model()
    example_large_language_model()
    example_supervised_model()

    print("\n" + "=" * 70)
    print("Examples completed!")
    print("=" * 70)
    print("\nGenerated files:")
    print("  - dms_corr_Q9NV35_AlphaMissense.png")
    print("  - dms_corr_P04637_GEMME.png")
    print("  - dms_corr_Q9NV35_ESM.png")
    print("  - dms_corr_Q9NV35_EVE.png")
    print("\nYou can use these functions in your own analysis:")
    print("  from proteingympy.dms_corr_plot import dms_corr_plot")
    print("  fig = dms_corr_plot(uniprot_id='YOUR_ID', model='MODEL_NAME')")
    print("  plt.show()")


if __name__ == "__main__":
    main()
