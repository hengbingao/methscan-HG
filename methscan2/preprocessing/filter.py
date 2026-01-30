"""
Filtering Functions
"""

import numpy as np
from typing import Optional


def filter_cells(
    mdata,
    min_sites: Optional[int] = None,
    max_sites: Optional[int] = None,
    min_coverage: Optional[float] = None,
    max_coverage: Optional[float] = None,
    min_methylation: Optional[float] = None,
    max_methylation: Optional[float] = None,
    inplace: bool = True
):
    """
    Filter cells based on QC metrics
    
    Parameters:
        mdata: MethylationData object
        min_sites: Minimum number of covered sites
        max_sites: Maximum number of covered sites
        min_coverage: Minimum mean coverage
        max_coverage: Maximum mean coverage
        min_methylation: Minimum mean methylation level
        max_methylation: Maximum mean methylation level
        inplace: Modify mdata in place or return filtered copy
        
    Returns:
        None if inplace=True, filtered MethylationData if inplace=False
    """
    # Ensure QC metrics are calculated
    required_metrics = []
    if min_sites is not None or max_sites is not None:
        required_metrics.append('n_sites')
    if min_coverage is not None or max_coverage is not None:
        required_metrics.append('mean_coverage')
    if min_methylation is not None or max_methylation is not None:
        required_metrics.append('mean_methylation')
    
    missing = [m for m in required_metrics if m not in mdata.obs.columns]
    if missing:
        raise ValueError(f"QC metrics not found: {missing}. Run pp.calculate_qc_metrics() first.")
    
    # Build filter mask
    mask = np.ones(mdata.n_obs, dtype=bool)
    
    if min_sites is not None:
        mask &= mdata.obs['n_sites'] >= min_sites
    if max_sites is not None:
        mask &= mdata.obs['n_sites'] <= max_sites
        
    if min_coverage is not None:
        mask &= mdata.obs['mean_coverage'] >= min_coverage
    if max_coverage is not None:
        mask &= mdata.obs['mean_coverage'] <= max_coverage
        
    if min_methylation is not None:
        mask &= mdata.obs['mean_methylation'] >= min_methylation
    if max_methylation is not None:
        mask &= mdata.obs['mean_methylation'] <= max_methylation
    
    n_filtered = (~mask).sum()
    n_kept = mask.sum()
    
    print(f"Filtered out {n_filtered} cells, kept {n_kept} cells")
    
    if inplace:
        mdata._inplace_subset_obs(mask)
        return None
    else:
        return mdata[mask, :].copy()


def filter_sites(
    mdata,
    min_cells: Optional[int] = None,
    min_coverage: Optional[int] = None,
    chromosomes: Optional[list] = None,
    inplace: bool = True
):
    """
    Filter genomic sites
    
    Parameters:
        mdata: MethylationData object
        min_cells: Minimum number of cells with coverage
        min_coverage: Minimum coverage in at least min_cells cells
        chromosomes: List of chromosomes to keep (e.g., ['chr1', 'chr2'])
        inplace: Modify in place or return copy
        
    Returns:
        None if inplace=True, filtered MethylationData if inplace=False
    """
    mask = np.ones(mdata.n_vars, dtype=bool)
    
    # Filter by chromosome
    if chromosomes is not None:
        if 'chrom' not in mdata.var.columns:
            raise ValueError("Chromosome information not available in mdata.var")
        mask &= mdata.var['chrom'].isin(chromosomes)
    
    # Filter by coverage
    if 'total' in mdata.layers:
        total = mdata.layers['total']
        
        if min_coverage is not None and min_cells is not None:
            # Sites with at least min_coverage in at least min_cells cells
            sufficient_coverage = total >= min_coverage
            n_cells_covered = sufficient_coverage.sum(axis=0)
            mask &= n_cells_covered >= min_cells
            
        elif min_cells is not None:
            # Sites covered in at least min_cells cells
            covered = total > 0
            n_cells_covered = covered.sum(axis=0)
            mask &= n_cells_covered >= min_cells
    
    n_filtered = (~mask).sum()
    n_kept = mask.sum()
    
    print(f"Filtered out {n_filtered} sites, kept {n_kept} sites")
    
    if inplace:
        mdata._inplace_subset_var(mask)
        return None
    else:
        return mdata[:, mask].copy()
