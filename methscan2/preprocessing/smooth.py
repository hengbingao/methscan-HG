"""
Methylation Smoothing Functions
"""

import numpy as np
import pandas as pd
from scipy.ndimage import uniform_filter1d
from joblib import Parallel, delayed
from typing import Optional


def smooth_methylation(
    mdata,
    bandwidth: int = 2000,
    chromosomes: Optional[list] = None,
    n_jobs: int = -1,
    inplace: bool = True
):
    """
    Calculate smoothed methylation curve for each chromosome
    
    Uses kernel smoothing to create a pseudo-bulk methylation profile.
    This is used as a reference for VMR detection.
    
    Parameters:
        mdata: MethylationData object
        bandwidth: Smoothing bandwidth in base pairs
        chromosomes: List of chromosomes to process. None = all
        n_jobs: Number of parallel jobs
        inplace: Store result in mdata.uns['smooth_methylation']
        
    Returns:
        Dictionary of smoothed curves if inplace=False
    """
    if chromosomes is None:
        chromosomes = mdata.chromosomes
    
    print(f"Smoothing methylation on {len(chromosomes)} chromosomes...")
    
    # Process each chromosome in parallel
    results = Parallel(n_jobs=n_jobs)(
        delayed(_smooth_chromosome)(mdata, chrom, bandwidth)
        for chrom in chromosomes
    )
    
    # Combine results
    smooth_dict = dict(zip(chromosomes, results))
    
    if inplace:
        mdata.uns['smooth_methylation'] = smooth_dict
        print("Smoothed curves stored in mdata.uns['smooth_methylation']")
        return None
    else:
        return smooth_dict


def _smooth_chromosome(mdata, chrom: str, bandwidth: int):
    """Smooth methylation on a single chromosome"""
    # Extract chromosome data
    chrom_mask = mdata.var['chrom'] == chrom
    
    if chrom_mask.sum() == 0:
        return None
    
    # Get positions
    positions = mdata.var.loc[chrom_mask, 'start'].values
    
    # Get methylation data
    meth_rate = mdata.methylation_rate[:, chrom_mask]
    
    # Get coverage
    if 'total' in mdata.layers:
        coverage = mdata.layers['total'][:, chrom_mask]
    else:
        coverage = np.ones_like(meth_rate)
    
    # Calculate weighted pseudobulk
    # Weight by coverage and number of cells
    valid_mask = ~np.isnan(meth_rate) & (coverage > 0)
    
    # Sum across cells
    total_met = np.nansum(
        np.where(valid_mask, meth_rate * coverage, 0), 
        axis=0
    )
    total_cov = np.nansum(
        np.where(valid_mask, coverage, 0),
        axis=0
    )
    
    # Calculate pseudobulk methylation
    with np.errstate(divide='ignore', invalid='ignore'):
        pseudobulk = total_met / total_cov
        pseudobulk[~np.isfinite(pseudobulk)] = np.nan
    
    # Apply smoothing
    # Use uniform filter for simplicity
    window_size = max(1, bandwidth // np.median(np.diff(positions)))
    
    smoothed = uniform_filter1d(
        np.nan_to_num(pseudobulk, nan=0),
        size=int(window_size),
        mode='nearest'
    )
    
    return {
        'position': positions,
        'methylation': smoothed,
        'pseudobulk': pseudobulk
    }
