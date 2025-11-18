#!/usr/bin/env python
"""Test script for export_structure_html function."""

from proteingympy.plot_structure import export_structure_html

# Test 1: Export DMS scores with full structure
print("Test 1: Exporting ACE2_HUMAN_Chan_2020 with DMS scores...")
export_structure_html(
    assay_name="ACE2_HUMAN_Chan_2020",
    output_path="docs/assets/ace2_dms_structure.html",
    full_structure=True,
    title="ACE2 Human - DMS Scores"
)

print("\nTest 1 completed successfully!")
