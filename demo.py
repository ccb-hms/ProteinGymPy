#!/usr/bin/env python3
"""
Simple demonstration of proteinGymPy functionality.

This shows the basic functions without downloading large datasets.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from proteingympy import (
    available_supervised_models,
    available_zero_shot_models,
    create_complete_metadata_table
)

def main():
    print("🧬 ProteinGymPy Demo")
    print("=" * 50)
    
    # Show available models
    print("\n📋 Available Models:")
    
    supervised = available_supervised_models()
    print(f"\nSupervised models ({len(supervised)}):")
    for i, model in enumerate(supervised[:5], 1):
        print(f"  {i}. {model}")
    if len(supervised) > 5:
        print(f"  ... and {len(supervised)-5} more")
    
    zeroshot = available_zero_shot_models()  
    print(f"\nZero-shot models ({len(zeroshot)}):")
    for i, model in enumerate(zeroshot[:5], 1):
        print(f"  {i}. {model}")
    if len(zeroshot) > 5:
        print(f"  ... and {len(zeroshot)-5} more")
    
    # Show metadata capabilities
    print(f"\n📊 Dataset Metadata:")
    metadata_df = create_complete_metadata_table()
    print(f"  Generated metadata for {len(metadata_df)} datasets")
    print(f"  Available data types: {', '.join(metadata_df['SourceType'].unique())}")
    
    print(f"\n📖 Sample metadata entries:")
    sample = metadata_df[['Title', 'SourceType', 'DataProvider']].head(3)
    for _, row in sample.iterrows():
        print(f"  • {row['Title'][:50]}...")
        print(f"    Type: {row['SourceType']}, Provider: {row['DataProvider']}")
    
    print(f"\n✅ All functions working correctly!")
    print(f"\nTo run full examples with data download:")
    print(f"  python examples/proteingym_pipeline_examples.py")
    print(f"\nTo run tests:")
    print(f"  python -m pytest tests/ -v")

if __name__ == "__main__":
    main()