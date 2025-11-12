#!/usr/bin/env python3
"""
Example script demonstrating model_corr_plot functionality.

This script shows how to compare prediction scores between different
ProteinGym models using Spearman correlation and visualization.
"""

import os
import sys

# Add src to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import matplotlib.pyplot as plt
from proteingympy.model_corr_plot import (
    model_corr_plot,
    get_available_models
)


def example_basic_usage():
    """
    Example 1: Basic usage with default models (AlphaMissense vs GEMME)
    """
    print("=" * 70)
    print("Example 1: Basic Usage - AlphaMissense vs GEMME")
    print("=" * 70)

    # Use a known UniProt ID from the example
    uniprot_id = "Q9NV35"

    print(f"\nComparing AlphaMissense and GEMME predictions for {uniprot_id}...")

    try:
        fig = model_corr_plot(
            uniprot_id=uniprot_id,
            model1="AlphaMissense",
            model2="GEMME"
        )

        # Save figure
        output_path = f"model_corr_{uniprot_id}_AlphaMissense_vs_GEMME.png"
        fig.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"\n✓ Figure saved to {output_path}")

        # Show interactively (comment out if running headless)
        # plt.show()

        plt.close(fig)

    except Exception as e:
        print(f"\n✗ Error: {e}")
        print("This is expected if the data is not yet downloaded.")


def example_custom_models():
    """
    Example 2: Compare two zero-shot models
    """
    print("\n" + "=" * 70)
    print("Example 2: Custom Models - ESM vs Tranception")
    print("=" * 70)

    # Use P04637 (TP53) as another example
    uniprot_id = "P04637"

    print(f"\nComparing ESM and Tranception predictions for {uniprot_id}...")

    try:
        fig = model_corr_plot(
            uniprot_id=uniprot_id,
            model1="ESM_2_t48_15B_UR50D",
            model2="Tranception_L",
            bins=40,  # Fewer bins for different visualization
            cmap='plasma'  # Different colormap
        )

        # Save figure
        output_path = f"model_corr_{uniprot_id}_ESM_vs_Tranception.png"
        fig.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"\n✓ Figure saved to {output_path}")

        plt.close(fig)

    except Exception as e:
        print(f"\n✗ Error: {e}")
        print("This is expected if the data is not yet downloaded.")


def example_supervised_models():
    """
    Example 3: Compare supervised models
    """
    print("\n" + "=" * 70)
    print("Example 3: Supervised Models - EVE vs DeepSequence")
    print("=" * 70)

    uniprot_id = "Q9NV35"

    print(f"\nComparing EVE and DeepSequence predictions for {uniprot_id}...")

    try:
        fig = model_corr_plot(
            uniprot_id=uniprot_id,
            model1="EVE",
            model2="DeepSequence",
            figsize=(12, 10),  # Larger figure
            cmap='magma'
        )

        # Save figure
        output_path = f"model_corr_{uniprot_id}_EVE_vs_DeepSequence.png"
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
    for i, model in enumerate(models['zero_shot'], 1):
        print(f"  {i:2d}. {model}")

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
    print("ProteinGym Model Correlation Plot Examples")
    print("=" * 70)
    print("\nThese examples demonstrate how to compare prediction scores")
    print("between different ProteinGym models using Spearman correlation.")
    print("\nNote: Examples require downloaded data. First run will download")
    print("data from Zenodo, which may take some time.")
    print("=" * 70)

    # Show available models
    show_available_models()

    # Run examples
    example_basic_usage()
    example_custom_models()
    example_supervised_models()

    print("\n" + "=" * 70)
    print("Examples completed!")
    print("=" * 70)
    print("\nGenerated files:")
    print("  - model_corr_Q9NV35_AlphaMissense_vs_GEMME.png")
    print("  - model_corr_P04637_ESM_vs_Tranception.png")
    print("  - model_corr_Q9NV35_EVE_vs_DeepSequence.png")
    print("\nYou can use these functions in your own analysis:")
    print("  from proteingympy.model_corr_plot import model_corr_plot")
    print("  fig = model_corr_plot(uniprot_id='YOUR_ID', model1='MODEL1', model2='MODEL2')")
    print("  plt.show()")


if __name__ == "__main__":
    main()
