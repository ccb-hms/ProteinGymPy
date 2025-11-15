#!/usr/bin/env python3
"""
Example script demonstrating benchmark_models functionality.

This script shows how to visualize one of the five performance metrics 
across different models available in ProteinGym.
"""

import os
import sys

# Add src to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import matplotlib.pyplot as plt
from proteingympy.benchmark_models import (
    benchmark_models,
    available_models
)


def example_basic_usage():
    """
    Example 1: Basic usage with default models (EVmutation vs GEMME)
    """
    print("=" * 70)
    print("Example 1: Basic Usage - EVmutation vs GEMME")
    print("=" * 70)

    metric = "AUC"

    print(f"\nComparing EVmutation and GEMME performance using {metric}...")

    try:
        fig = benchmark_models(
            metric=metric,
            models=["EVmutation", "GEMME"]
        )

        output_path = f"benchmark_{metric}_EVmut_vs_GEMME.png"
        fig.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"\n✓ Figure saved to {output_path}")

        # plt.show()  # Uncomment for interactive use

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

    models = available_models()

    print(f"\nTotal: {len(models)} models available")
    print("\nAvailable models:")
    for name in models:
        print(f"  {name}")


def main():
    """
    Run all examples
    """
    print("ProteinGym Model Correlation Plot Examples")
    print("=" * 70)
    print("\nThese examples demonstrate how to compare prediction scores")
    print("between different ProteinGym models using performance metrics (e.g. Spearman correlation).")
    print("\nNote: Examples require downloaded data. The first run may download")
    print("data from Zenodo, which can take several minutes.")
    print("=" * 70)

    # Show available models
    show_available_models()

    # Run examples
    example_basic_usage()

    print("\n" + "=" * 70)
    print("Examples completed!")
    print("=" * 70)
    print("\nGenerated files:")
    print("  - benchmark_AUC_EVmut_vs_GEMME.png")
    print("\nYou can use these functions in your own analysis:")
    print("  from proteingympy.benchmark_models import benchmark_models")
    print("  fig = benchmark_models(metric='ONE_METRIC', models=['LIST', 'OF', 'MODELS'])")
    print("  plt.show()")

if __name__ == "__main__":
    main()