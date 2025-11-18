"""Test script for plot_structure with full_structure=True"""

from proteingympy.plot_structure import plot_structure

# Test the function
print("Testing plot_structure with full_structure=True...")
view = plot_structure(
    assay_name="ACE2_HUMAN_Chan_2020",
    data_scores="DMS",
    color_scheme="EVE",
    full_structure=True
)

print(f"\nView type: {type(view)}")
print(f"Number of representations: {len(view._representations)}")

if view._representations:
    print("\nRepresentations added:")
    for i, rep in enumerate(view._representations):
        print(f"  {i}: {rep}")
else:
    print("\nWARNING: No representations were added!")

print("\nIf running in a Jupyter environment, the view object should display the 3D structure.")
print("Grey residues should show the full structure, with colored residues overlaid.")
