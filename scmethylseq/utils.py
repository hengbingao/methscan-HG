"""
Utility functions for scMethylSeq
"""

import os
import logging
from pathlib import Path
from typing import List, Dict, Optional, Tuple
import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


def validate_file_exists(filepath: str, description: str = "File") -> Path:
    """
    Validate that a file exists
    
    Parameters:
    -----------
    filepath : str
        Path to file
    description : str
        Description for error message
        
    Returns:
    --------
    Path
        Validated Path object
        
    Raises:
    -------
    FileNotFoundError
        If file doesn't exist
    """
    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"{description} not found: {filepath}")
    return path


def create_output_directory(dirpath: str, exist_ok: bool = True) -> Path:
    """
    Create output directory
    
    Parameters:
    -----------
    dirpath : str
        Directory path
    exist_ok : bool
        Whether to allow existing directory
        
    Returns:
    --------
    Path
        Created directory path
    """
    path = Path(dirpath)
    path.mkdir(parents=True, exist_ok=exist_ok)
    logger.info(f"Created output directory: {path}")
    return path


def count_files_in_directory(dirpath: str, pattern: str = "*") -> int:
    """
    Count files matching pattern in directory
    
    Parameters:
    -----------
    dirpath : str
        Directory path
    pattern : str
        File pattern (glob)
        
    Returns:
    --------
    int
        Number of files
    """
    path = Path(dirpath)
    if not path.exists():
        return 0
    return len(list(path.glob(pattern)))


def get_file_size_gb(filepath: str) -> float:
    """
    Get file size in GB
    
    Parameters:
    -----------
    filepath : str
        File path
        
    Returns:
    --------
    float
        File size in GB
    """
    path = Path(filepath)
    if not path.exists():
        return 0.0
    return path.stat().st_size / (1024 ** 3)


def check_dependencies() -> Dict[str, bool]:
    """
    Check if required system dependencies are available
    
    Returns:
    --------
    Dict[str, bool]
        Dictionary of dependency name to availability
    """
    import subprocess
    import shutil
    
    dependencies = {
        'bedtools': shutil.which('bedtools') is not None,
        'tabix': shutil.which('tabix') is not None,
        'allcools': shutil.which('allcools') is not None,
    }
    
    return dependencies


def print_dependency_status():
    """
    Print status of system dependencies
    """
    deps = check_dependencies()
    
    print("\n" + "=" * 60)
    print("System Dependencies Status")
    print("=" * 60)
    
    for dep, available in deps.items():
        status = "✓ Available" if available else "✗ Not found"
        print(f"{dep:20s}: {status}")
    
    if not all(deps.values()):
        print("\nWarning: Some dependencies are missing.")
        print("Please install missing dependencies before running the pipeline.")
    
    print("=" * 60 + "\n")


def estimate_memory_usage(n_cells: int, n_genes: int) -> Dict[str, float]:
    """
    Estimate memory usage for analysis
    
    Parameters:
    -----------
    n_cells : int
        Number of cells
    n_genes : int
        Number of genes
        
    Returns:
    --------
    Dict[str, float]
        Memory estimates in GB
    """
    # Rough estimates based on typical data
    matrix_gb = (n_cells * n_genes * 8) / (1024 ** 3)  # float64
    
    estimates = {
        'methylation_matrix': matrix_gb,
        'pca_computation': matrix_gb * 1.5,
        'umap_computation': (n_cells * 100 * 8) / (1024 ** 3),  # 100 neighbors
        'total_recommended': matrix_gb * 3,
    }
    
    return estimates


def print_memory_estimate(n_cells: int, n_genes: int):
    """
    Print memory usage estimates
    
    Parameters:
    -----------
    n_cells : int
        Number of cells
    n_genes : int
        Number of genes
    """
    estimates = estimate_memory_usage(n_cells, n_genes)
    
    print("\n" + "=" * 60)
    print("Memory Usage Estimates")
    print("=" * 60)
    print(f"Dataset: {n_cells:,} cells × {n_genes:,} genes")
    print()
    print(f"Methylation matrix:     {estimates['methylation_matrix']:.2f} GB")
    print(f"PCA computation:        {estimates['pca_computation']:.2f} GB")
    print(f"UMAP computation:       {estimates['umap_computation']:.2f} GB")
    print()
    print(f"Recommended RAM:        {estimates['total_recommended']:.2f} GB")
    print("=" * 60 + "\n")


def merge_metadata_files(file_list: List[str], output_file: str):
    """
    Merge multiple metadata files
    
    Parameters:
    -----------
    file_list : List[str]
        List of metadata file paths
    output_file : str
        Output merged file
    """
    logger.info(f"Merging {len(file_list)} metadata files...")
    
    dfs = []
    for f in file_list:
        df = pd.read_csv(f)
        dfs.append(df)
    
    # Merge on 'cell' column
    merged = dfs[0]
    for df in dfs[1:]:
        merged = pd.merge(merged, df, on='cell', how='outer')
    
    merged.to_csv(output_file, index=False)
    logger.info(f"Merged metadata saved to {output_file}")


def filter_matrix_by_cells(matrix_file: str,
                           cell_list: List[str],
                           output_file: str):
    """
    Filter methylation matrix to keep only specified cells
    
    Parameters:
    -----------
    matrix_file : str
        Input matrix file
    cell_list : List[str]
        List of cells to keep
    output_file : str
        Output filtered matrix
    """
    logger.info(f"Filtering matrix to {len(cell_list)} cells...")
    
    # Load matrix
    df = pd.read_csv(matrix_file, sep='\t', header=None)
    
    # Get cell names from first row
    cell_names = df.iloc[0, 1:].values
    
    # Find indices of cells to keep
    keep_indices = [0]  # Always keep first column (gene names)
    for i, cell in enumerate(cell_names):
        if cell in cell_list:
            keep_indices.append(i + 1)
    
    # Filter columns
    filtered = df.iloc[:, keep_indices]
    
    # Save
    filtered.to_csv(output_file, sep='\t', header=False, index=False)
    logger.info(f"Filtered matrix saved to {output_file}")


def compute_methylation_statistics(adata, 
                                   context: str = 'CG') -> pd.DataFrame:
    """
    Compute methylation statistics per cell
    
    Parameters:
    -----------
    adata : AnnData
        AnnData object
    context : str
        Methylation context ('CG' or 'CH')
        
    Returns:
    --------
    pd.DataFrame
        Statistics per cell
    """
    from anndata import AnnData
    
    if context == 'CG':
        key = 'X_cg'
    else:
        key = 'X_ch'
    
    if key not in adata.obsm:
        raise ValueError(f"{key} not found in adata.obsm")
    
    data = adata.obsm[key]
    
    # Exclude missing values (2)
    valid_data = data.copy()
    valid_data[valid_data == 2] = np.nan
    
    stats = pd.DataFrame({
        'cell': adata.obs_names,
        f'mean_{context}': np.nanmean(valid_data, axis=1),
        f'median_{context}': np.nanmedian(valid_data, axis=1),
        f'std_{context}': np.nanstd(valid_data, axis=1),
        f'n_valid_genes': np.sum(~np.isnan(valid_data), axis=1),
    })
    
    return stats


def export_cluster_markers(adata,
                          cluster_key: str = 'leiden',
                          output_dir: str = 'markers',
                          n_genes: int = 100):
    """
    Export top marker genes for each cluster
    
    Parameters:
    -----------
    adata : AnnData
        AnnData object with clustering
    cluster_key : str
        Key in obs for cluster assignments
    output_dir : str
        Output directory
    n_genes : int
        Number of top genes per cluster
    """
    try:
        import scanpy as sc
    except ImportError:
        logger.error("scanpy required for marker gene analysis")
        return
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    logger.info(f"Computing marker genes for {cluster_key}...")
    
    # Compute differential genes
    sc.tl.rank_genes_groups(adata, groupby=cluster_key, method='wilcoxon')
    
    # Export for each cluster
    clusters = adata.obs[cluster_key].unique()
    
    for cluster in clusters:
        marker_df = sc.get.rank_genes_groups_df(adata, group=str(cluster))
        marker_df = marker_df.head(n_genes)
        
        output_file = output_dir / f"cluster_{cluster}_markers.csv"
        marker_df.to_csv(output_file, index=False)
        
        logger.info(f"Cluster {cluster}: {len(marker_df)} markers saved")


def convert_formats(input_file: str,
                   output_file: str,
                   input_format: str = 'csv',
                   output_format: str = 'h5ad'):
    """
    Convert between different file formats
    
    Parameters:
    -----------
    input_file : str
        Input file path
    output_file : str
        Output file path
    input_format : str
        Input format ('csv', 'tsv', 'h5ad')
    output_format : str
        Output format ('csv', 'tsv', 'h5ad')
    """
    if input_format == output_format:
        logger.warning("Input and output formats are the same")
        return
    
    logger.info(f"Converting {input_format} to {output_format}...")
    
    # Implement conversion logic
    # This is a placeholder for format conversion
    raise NotImplementedError("Format conversion not yet implemented")


class ProgressTracker:
    """
    Track pipeline progress
    """
    
    def __init__(self, output_file: str = 'pipeline_progress.txt'):
        self.output_file = output_file
        self.steps = {}
        self.load_progress()
    
    def load_progress(self):
        """Load existing progress"""
        if Path(self.output_file).exists():
            with open(self.output_file, 'r') as f:
                for line in f:
                    if ':' in line:
                        step, status = line.strip().split(':', 1)
                        self.steps[step] = status.strip()
    
    def mark_complete(self, step: str):
        """Mark step as complete"""
        self.steps[step] = 'COMPLETE'
        self.save_progress()
        logger.info(f"✓ Completed: {step}")
    
    def mark_failed(self, step: str, error: str = ''):
        """Mark step as failed"""
        self.steps[step] = f'FAILED: {error}'
        self.save_progress()
        logger.error(f"✗ Failed: {step}")
    
    def is_complete(self, step: str) -> bool:
        """Check if step is complete"""
        return self.steps.get(step) == 'COMPLETE'
    
    def save_progress(self):
        """Save progress to file"""
        with open(self.output_file, 'w') as f:
            for step, status in self.steps.items():
                f.write(f"{step}: {status}\n")
    
    def print_summary(self):
        """Print progress summary"""
        print("\n" + "=" * 60)
        print("Pipeline Progress Summary")
        print("=" * 60)
        
        for step, status in self.steps.items():
            print(f"{step:30s}: {status}")
        
        print("=" * 60 + "\n")
