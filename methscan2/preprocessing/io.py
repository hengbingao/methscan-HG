"""
Data I/O Functions for MethSCAn2
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Optional, Union, List
from glob import glob
from joblib import Parallel, delayed
from tqdm.auto import tqdm
import warnings

from ..core.methylation_data import MethylationData


def read_cov_files(
    file_paths: Union[str, List[str]],
    sample_names: Optional[List[str]] = None,
    genome: str = 'hg38',
    min_coverage: int = 1,
    n_jobs: int = -1,
    verbose: bool = True
) -> MethylationData:
    """
    Read Bismark coverage files (.cov format)
    
    Format: chr start end methylation_percentage count_methylated count_unmethylated
    
    Parameters:
        file_paths: Path pattern (with wildcards) or list of file paths
        sample_names: Sample names for each file. If None, use filenames
        genome: Genome assembly name (e.g., 'hg38', 'mm10')
        min_coverage: Minimum read coverage to include a site
        n_jobs: Number of parallel jobs. -1 uses all CPUs
        verbose: Show progress bar
        
    Returns:
        MethylationData object
    """
    # Parse file paths
    if isinstance(file_paths, str):
        files = sorted(glob(file_paths))
    else:
        files = file_paths
    
    if len(files) == 0:
        raise ValueError(f"No files found matching pattern: {file_paths}")
    
    if verbose:
        print(f"Found {len(files)} files to read")
    
    # Generate sample names if not provided
    if sample_names is None:
        sample_names = [Path(f).stem for f in files]
    
    if len(sample_names) != len(files):
        raise ValueError(f"Number of sample names ({len(sample_names)}) "
                        f"doesn't match number of files ({len(files)})")
    
    # Read files in parallel
    if verbose:
        print("Reading files...")
    
    results = Parallel(n_jobs=n_jobs)(
        delayed(_read_single_cov)(f, name, min_coverage)
        for f, name in (tqdm(list(zip(files, sample_names)), desc="Loading") 
                       if verbose else zip(files, sample_names))
    )
    
    # Combine results
    if verbose:
        print("Combining data from all samples...")
    
    mdata = _combine_cov_results(results, genome, verbose)
    
    if verbose:
        print(f"Loaded {mdata.n_obs} cells × {mdata.n_vars} sites")
    
    return mdata


def _read_single_cov(filepath: str, sample_name: str, min_coverage: int):
    """Read a single .cov file"""
    try:
        df = pd.read_csv(
            filepath,
            sep='\t',
            header=None,
            names=['chrom', 'start', 'end', 'meth_pct', 'met_count', 'unmet_count'],
            dtype={
                'chrom': str,
                'start': int,
                'end': int,
                'meth_pct': float,
                'met_count': int,
                'unmet_count': int
            }
        )
        
        # Calculate total coverage
        df['total_count'] = df['met_count'] + df['unmet_count']
        
        # Filter by minimum coverage
        df = df[df['total_count'] >= min_coverage].copy()
        
        # Create position key
        df['position'] = df['chrom'] + ':' + df['start'].astype(str)
        
        return {
            'sample_name': sample_name,
            'data': df,
            'n_sites': len(df)
        }
    except Exception as e:
        warnings.warn(f"Error reading {filepath}: {e}")
        return None


def _combine_cov_results(results, genome, verbose):
    """Combine results from multiple samples into MethylationData"""
    # Remove failed reads
    results = [r for r in results if r is not None]
    
    if len(results) == 0:
        raise ValueError("No data could be read from files")
    
    # Get all unique positions across all samples
    all_positions = set()
    for r in results:
        all_positions.update(r['data']['position'].values)
    
    all_positions = sorted(list(all_positions))
    
    if verbose:
        print(f"Found {len(all_positions)} unique genomic positions")
    
    # Create position to index mapping
    pos_to_idx = {pos: idx for idx, pos in enumerate(all_positions)}
    
    # Initialize matrices
    n_cells = len(results)
    n_sites = len(all_positions)
    
    met_matrix = np.zeros((n_cells, n_sites), dtype=np.int32)
    total_matrix = np.zeros((n_cells, n_sites), dtype=np.int32)
    
    # Fill matrices
    for cell_idx, result in enumerate(results):
        df = result['data']
        for _, row in df.iterrows():
            site_idx = pos_to_idx[row['position']]
            met_matrix[cell_idx, site_idx] = row['met_count']
            total_matrix[cell_idx, site_idx] = row['total_count']
    
    # Calculate methylation rate
    with np.errstate(divide='ignore', invalid='ignore'):
        rate_matrix = met_matrix.astype(float) / total_matrix
        rate_matrix[~np.isfinite(rate_matrix)] = np.nan
    
    # Create var (site) annotations
    var_df = pd.DataFrame({'position': all_positions})
    var_df[['chrom', 'start']] = var_df['position'].str.split(':', expand=True)
    var_df['start'] = var_df['start'].astype(int)
    var_df['end'] = var_df['start'] + 1
    var_df = var_df.drop(columns=['position'])
    var_df.index = all_positions
    
    # Create obs (cell) annotations
    obs_df = pd.DataFrame({
        'sample_name': [r['sample_name'] for r in results],
        'n_sites_raw': [r['n_sites'] for r in results]
    })
    obs_df.index = obs_df['sample_name'].values
    
    # Create MethylationData object
    mdata = MethylationData(
        X=rate_matrix,
        obs=obs_df,
        var=var_df,
        layers={
            'met': met_matrix,
            'total': total_matrix,
            'rate': rate_matrix
        },
        uns={'genome': genome}
    )
    
    return mdata


def read_allc_files(
    file_paths: Union[str, List[str]],
    sample_names: Optional[List[str]] = None,
    context: str = 'CG',
    genome: str = 'hg38',
    min_coverage: int = 1,
    n_jobs: int = -1,
    verbose: bool = True
) -> MethylationData:
    """
    Read ALLC format files (methylpy output)
    
    Format: chr pos strand context mc_count total_count methylation_level
    
    Parameters:
        file_paths: Path pattern or list of file paths
        sample_names: Sample names for each file
        context: Methylation context ('CG', 'CHG', 'CHH', or 'all')
        genome: Genome assembly name
        min_coverage: Minimum read coverage
        n_jobs: Number of parallel jobs
        verbose: Show progress
        
    Returns:
        MethylationData object
    """
    # Parse file paths
    if isinstance(file_paths, str):
        files = sorted(glob(file_paths))
    else:
        files = file_paths
    
    if len(files) == 0:
        raise ValueError(f"No files found: {file_paths}")
    
    if verbose:
        print(f"Found {len(files)} ALLC files")
    
    # Generate sample names
    if sample_names is None:
        sample_names = [Path(f).stem for f in files]
    
    # Read files
    results = Parallel(n_jobs=n_jobs)(
        delayed(_read_single_allc)(f, name, context, min_coverage)
        for f, name in (tqdm(list(zip(files, sample_names)), desc="Loading") 
                       if verbose else zip(files, sample_names))
    )
    
    # Combine results (similar to cov files)
    mdata = _combine_cov_results(results, genome, verbose)
    
    return mdata


def _read_single_allc(filepath: str, sample_name: str, context: str, min_coverage: int):
    """Read a single ALLC file"""
    try:
        df = pd.read_csv(
            filepath,
            sep='\t',
            header=None,
            names=['chrom', 'start', 'strand', 'context', 'met_count', 'total_count', 'meth_rate']
        )
        
        # Filter by context if specified
        if context != 'all':
            df = df[df['context'].str.startswith(context)].copy()
        
        # Filter by coverage
        df = df[df['total_count'] >= min_coverage].copy()
        
        # Calculate unmethylated count
        df['unmet_count'] = df['total_count'] - df['met_count']
        
        # Add end position
        df['end'] = df['start'] + 1
        
        # Create position key
        df['position'] = df['chrom'] + ':' + df['start'].astype(str)
        
        return {
            'sample_name': sample_name,
            'data': df,
            'n_sites': len(df)
        }
    except Exception as e:
        warnings.warn(f"Error reading {filepath}: {e}")
        return None


def read_bedgraph(
    file_paths: Union[str, List[str]],
    sample_names: Optional[List[str]] = None,
    genome: str = 'hg38',
    n_jobs: int = -1,
    verbose: bool = True
) -> MethylationData:
    """
    Read bedGraph format files
    
    Format: chr start end methylation_value
    
    Parameters:
        file_paths: Path pattern or list of file paths
        sample_names: Sample names
        genome: Genome assembly
        n_jobs: Parallel jobs
        verbose: Show progress
        
    Returns:
        MethylationData object
    """
    if isinstance(file_paths, str):
        files = sorted(glob(file_paths))
    else:
        files = file_paths
    
    if sample_names is None:
        sample_names = [Path(f).stem for f in files]
    
    if verbose:
        print(f"Reading {len(files)} bedGraph files...")
    
    # Read files
    dfs = Parallel(n_jobs=n_jobs)(
        delayed(pd.read_csv)(f, sep='\t', header=None, 
                            names=['chrom', 'start', 'end', 'value'])
        for f in (tqdm(files, desc="Loading") if verbose else files)
    )
    
    # Combine into matrix
    # ... implementation similar to above
    
    raise NotImplementedError("bedGraph reader coming soon")


def read_h5ad(filename: str, backed: Optional[str] = None) -> MethylationData:
    """
    Read MethylationData from HDF5 file
    
    Parameters:
        filename: Path to .h5ad file
        backed: 'r' for read-only backed mode, None for loading into memory
        
    Returns:
        MethylationData object
    """
    from anndata import read_h5ad as anndata_read
    
    adata = anndata_read(filename, backed=backed)
    
    # Convert to MethylationData
    mdata = MethylationData(
        X=adata.X,
        obs=adata.obs,
        var=adata.var,
        uns=adata.uns,
        obsm=adata.obsm if hasattr(adata, 'obsm') else None,
        varm=adata.varm if hasattr(adata, 'varm') else None,
        layers=adata.layers if hasattr(adata, 'layers') else None,
        obsp=adata.obsp if hasattr(adata, 'obsp') else None,
        varp=adata.varp if hasattr(adata, 'varp') else None,
    )
    
    return mdata
