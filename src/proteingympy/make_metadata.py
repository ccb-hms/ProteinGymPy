"""
make_metadata.py - Python equivalent of make_metadata.R

Creates metadata for ProteinGym datasets including descriptions, sources, and data types.
"""

import os
import pandas as pd
from typing import Dict, List, Optional
from dataclasses import dataclass
import glob


@dataclass
class DatasetMetadata:
    """Metadata for a ProteinGym dataset."""
    title: str
    description: str
    bioc_version: str
    genome: Optional[str]
    source_type: str
    source_url: str
    source_version: Optional[str]
    species: Optional[str]
    taxonomy_id: Optional[str]
    coordinate_1_based: bool
    data_provider: str
    maintainer: str
    data_class: str
    dispatch_class: str
    data_path: str


def create_alphamissense_metadata() -> DatasetMetadata:
    """Create metadata for AlphaMissense supplementary data."""
    return DatasetMetadata(
        title="AlphaMissense pathogenicity scores for variants in ProteinGym",
        description="Supplementary table from Cheng et al. 2023 containing AlphaMissense pathogenicity scores for mutations found in ProteinGym DMS substitution data",
        bioc_version="3.20",
        genome=None,
        source_type="CSV",
        source_url="https://www.science.org/doi/10.1126/science.adg7492",
        source_version=None,
        species=None,
        taxonomy_id=None,
        coordinate_1_based=True,
        data_provider="Cheng et al. 2023",
        maintainer="ProteinGymPy Package <noreply@example.com>",
        data_class="DataFrame",
        dispatch_class="CSV",
        data_path="ProteinGymPy/cheng_proteingym_variants.csv"
    )


def create_dms_substitutions_metadata() -> DatasetMetadata:
    """Create metadata for DMS substitution scores."""
    return DatasetMetadata(
        title="ProteinGym deep mutational scanning (DMS) assays for substitutions",
        description="ProteinGym DMS information for 217 assays from Notin et al. 2023",
        bioc_version="3.20",
        genome=None,
        source_type="ZIP",
        source_url="https://proteingym.org/",
        source_version="1.1",
        species=None,
        taxonomy_id=None,
        coordinate_1_based=True,
        data_provider="Marks Lab at Harvard Medical School",
        maintainer="ProteinGymPy Package <noreply@example.com>",
        data_class="Dict",
        dispatch_class="Pickle",
        data_path="ProteinGymPy/progym217_dms_subs_v1_1.pkl"
    )


def create_dms_reference_metadata() -> DatasetMetadata:
    """Create metadata for DMS substitution reference file."""
    return DatasetMetadata(
        title="ProteinGym metadata for 217 DMS substitution assays",
        description="Reference file for ProteinGym v1.1 217 DMS assays from Notin et al. 2023",
        bioc_version="3.20",
        genome=None,
        source_type="CSV",
        source_url="https://proteingym.org/",
        source_version="1.1",
        species=None,
        taxonomy_id=None,
        coordinate_1_based=True,
        data_provider="Marks Lab at Harvard Medical School",
        maintainer="ProteinGymPy Package <noreply@example.com>",
        data_class="DataFrame",
        dispatch_class="CSV",
        data_path="ProteinGymPy/ref_file_217_dms_subs_v1_1.csv"
    )


def create_zeroshot_scores_metadata() -> DatasetMetadata:
    """Create metadata for zero-shot DMS substitution benchmarks."""
    return DatasetMetadata(
        title="ProteinGym zero-shot DMS substitution benchmarks",
        description="Zero-shot DMS substitution benchmarks from Notin et al. 2023 using Spearman, NDCG, AUC, MCC, and Top-K recall metrics",
        bioc_version="3.20",
        genome=None,
        source_type="CSV",
        source_url="https://proteingym.org/",
        source_version="1.1",
        species=None,
        taxonomy_id=None,
        coordinate_1_based=True,
        data_provider="Marks Lab at Harvard Medical School",
        maintainer="ProteinGymPy Package <noreply@example.com>",
        data_class="Dict",
        dispatch_class="Pickle",
        data_path="ProteinGymPy/zeroshot_dms_subs_v1_1.pkl"
    )


def create_zeroshot_summary_v12_metadata() -> DatasetMetadata:
    """Create metadata for zero-shot summary scores v1.2."""
    return DatasetMetadata(
        title="ProteinGym zero-shot DMS substitution benchmarks v1.2",
        description="Zero-shot DMS substitution benchmarks from Notin et al. 2023 using Spearman, NDCG, AUC, MCC, and Top-K recall metrics for 79 models",
        bioc_version="3.21",
        genome=None,
        source_type="CSV",
        source_url="https://zenodo.org/records/14997691",
        source_version="1.2",
        species=None,
        taxonomy_id=None,
        coordinate_1_based=True,
        data_provider="Marks Lab at Harvard Medical School",
        maintainer="ProteinGymPy Package <noreply@example.com>",
        data_class="Dict",
        dispatch_class="Pickle",
        data_path="ProteinGymPy/zeroshot_summary_scores_v1_2.pkl"
    )


def create_zeroshot_model_scores_metadata() -> DatasetMetadata:
    """Create metadata for zero-shot model scores."""
    return DatasetMetadata(
        title="ProteinGym zero-shot DMS substitution scores for 79 models",
        description="Zero-shot substitution scores for 79 models across 217 DMS assays from Notin et al. 2023",
        bioc_version="3.21",
        genome=None,
        source_type="CSV",
        source_url="https://zenodo.org/records/14997691",
        source_version="1.2",
        species=None,
        taxonomy_id=None,
        coordinate_1_based=True,
        data_provider="Marks Lab at Harvard Medical School",
        maintainer="ProteinGymPy Package <noreply@example.com>",
        data_class="Dict",
        dispatch_class="Pickle",
        data_path="ProteinGymPy/zeroshot_scores_v1_2.pkl"
    )


def create_supervised_metadata(fold_type: str) -> DatasetMetadata:
    """Create metadata for supervised model scores."""
    return DatasetMetadata(
        title=f"ProteinGym semi-supervised model prediction scores for 12 models ({fold_type})",
        description=f"Semi-supervised prediction scores for 12 models across 217 DMS assays with {fold_type} 5 variable folds from Notin et al. 2023",
        bioc_version="3.21",
        genome=None,
        source_type="CSV",
        source_url="https://zenodo.org/records/14997691",
        source_version="1.2",
        species=None,
        taxonomy_id=None,
        coordinate_1_based=True,
        data_provider="Marks Lab at Harvard Medical School",
        maintainer="ProteinGymPy Package <noreply@example.com>",
        data_class="Dict",
        dispatch_class="Pickle",
        data_path=f"ProteinGymPy/supervised_{fold_type}5_scores_v1_2.pkl"
    )


def create_supervised_summary_metadata() -> DatasetMetadata:
    """Create metadata for supervised summary metrics."""
    return DatasetMetadata(
        title="ProteinGym summary metrics of semi-supervised scores for 12 models",
        description="Performance metrics of 12 semi-supervised models across 217 DMS assays using 5 variable fold for contiguous, modulo, and random settings from Notin et al. 2023",
        bioc_version="3.21",
        genome=None,
        source_type="CSV",
        source_url="https://zenodo.org/records/14997691",
        source_version="1.2",
        species=None,
        taxonomy_id=None,
        coordinate_1_based=True,
        data_provider="Marks Lab at Harvard Medical School",
        maintainer="ProteinGymPy Package <noreply@example.com>",
        data_class="DataFrame",
        dispatch_class="CSV",
        data_path="ProteinGymPy/supervised_summary_scores_v1_2.csv"
    )


def generate_pdb_metadata(pdb_directory: str) -> List[DatasetMetadata]:
    """
    Generate metadata for PDB structure files.
    
    Args:
        pdb_directory: Directory containing PDB files
        
    Returns:
        List of DatasetMetadata objects for each PDB file
    """
    metadata_list = []
    
    if os.path.exists(pdb_directory):
        pdb_files = glob.glob(os.path.join(pdb_directory, "*.pdb"))
        
        for pdb_file in pdb_files:
            filename = os.path.basename(pdb_file)
            pdb_id = os.path.splitext(filename)[0]
            
            metadata = DatasetMetadata(
                title=f"Protein structure for {pdb_id}",
                description=f"AlphaFold2 predicted protein structure for {pdb_id} from ProteinGym v1.2 curated by Notin et al. 2023",
                bioc_version="3.21",
                genome=None,
                source_type="PDB",
                source_url="https://zenodo.org/records/14997691",
                source_version="1.2",
                species=None,
                taxonomy_id=None,
                coordinate_1_based=True,
                data_provider="Marks Lab at Harvard Medical School",
                maintainer="ProteinGymPy Package <noreply@example.com>",
                data_class="str",
                dispatch_class="FilePath",
                data_path=f"ProteinGymPy/{filename}"
            )
            
            metadata_list.append(metadata)
    
    return metadata_list


def create_complete_metadata_table(pdb_directory: Optional[str] = None) -> pd.DataFrame:
    """
    Create complete metadata table for all ProteinGym datasets.
    
    Args:
        pdb_directory: Optional directory containing PDB files
        
    Returns:
        DataFrame with all metadata entries
    """
    metadata_entries = []
    
    # Core datasets
    metadata_entries.append(create_alphamissense_metadata())
    metadata_entries.append(create_dms_substitutions_metadata())
    metadata_entries.append(create_dms_reference_metadata())
    metadata_entries.append(create_zeroshot_scores_metadata())
    metadata_entries.append(create_zeroshot_summary_v12_metadata())
    metadata_entries.append(create_zeroshot_model_scores_metadata())
    
    # Supervised datasets
    for fold_type in ["contiguous", "modulo", "random"]:
        metadata_entries.append(create_supervised_metadata(fold_type))
    
    metadata_entries.append(create_supervised_summary_metadata())
    
    # PDB files if directory provided
    if pdb_directory:
        pdb_metadata = generate_pdb_metadata(pdb_directory)
        metadata_entries.extend(pdb_metadata)
    
    # Convert to DataFrame
    rows = []
    for metadata in metadata_entries:
        row = {
            'Title': metadata.title,
            'Description': metadata.description,
            'BiocVersion': metadata.bioc_version,
            'Genome': metadata.genome,
            'SourceType': metadata.source_type,
            'SourceUrl': metadata.source_url,
            'SourceVersion': metadata.source_version,
            'Species': metadata.species,
            'TaxonomyId': metadata.taxonomy_id,
            'Coordinate_1_based': metadata.coordinate_1_based,
            'DataProvider': metadata.data_provider,
            'Maintainer': metadata.maintainer,
            'DataClass': metadata.data_class,
            'DispatchClass': metadata.dispatch_class,
            'DataPath': metadata.data_path
        }
        rows.append(row)
    
    return pd.DataFrame(rows)


def save_metadata_csv(output_path: str, pdb_directory: Optional[str] = None):
    """
    Create and save complete metadata CSV file.
    
    Args:
        output_path: Path where to save the metadata CSV
        pdb_directory: Optional directory containing PDB files
    """
    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Create metadata table
    metadata_df = create_complete_metadata_table(pdb_directory)
    
    # Save to CSV
    metadata_df.to_csv(output_path, index=False)
    print(f"Metadata saved to {output_path}")
    print(f"Total datasets: {len(metadata_df)}")


if __name__ == "__main__":
    # Example usage
    print("Creating ProteinGym metadata...")
    
    # Create metadata table
    metadata_df = create_complete_metadata_table()
    
    print(f"Created metadata for {len(metadata_df)} datasets")
    print("\nSample entries:")
    print(metadata_df[['Title', 'SourceType', 'DataProvider']].head(10))
    
    # Optionally save to CSV
    output_path = "metadata/proteingym_metadata.csv"
    save_metadata_csv(output_path)