"""
Quality Control Functions
"""

import numpy as np
import pandas as pd
from typing import Optional, List


def calculate_qc_metrics(
    mdata,
    inplace: bool = True
) -> Optional[pd.DataFrame]:
    """
    Calculate quality control metrics for each cell
    
    Metrics calculated:
    - n_sites: Number of covered CpG sites
    - total_coverage: Total number of reads
    - mean_coverage: Average coverage per site
    - median_coverage: Median coverage per site
    - mean_methylation: Global methylation level
    - methylation_variance: Variance in methylation
    
    Parameters:
        mdata: MethylationData object
        inplace: Add metrics to mdata.obs if True
        
    Returns:
        DataFrame with QC metrics if inplace=False
    """
    if 'total' not in mdata.layers:
        raise ValueError("'total' layer not found. Cannot calculate coverage metrics.")
    
    total = mdata.layers['total']
    rate = mdata.methylation_rate
    
    # Calculate metrics
    metrics = pd.DataFrame(index=mdata.obs_names)
    
    # Number of covered sites
    metrics['n_sites'] = (total > 0).sum(axis=1)
    
    # Coverage metrics
    metrics['total_coverage'] = total.sum(axis=1)
    
    with np.errstate(divide='ignore', invalid='ignore'):
        # Mean coverage (only non-zero sites)
        total_nonzero = np.where(total > 0, total, np.nan)
        metrics['mean_coverage'] = np.nanmean(total_nonzero, axis=1)
        
        # Median coverage
        metrics['median_coverage'] = np.nanmedian(total_nonzero, axis=1)
    
    # Methylation metrics
    metrics['mean_methylation'] = np.nanmean(rate, axis=1)
    metrics['median_methylation'] = np.nanmedian(rate, axis=1)
    metrics['methylation_variance'] = np.nanvar(rate, axis=1)
    
    if inplace:
        # Add to obs
        for col in metrics.columns:
            mdata.obs[col] = metrics[col].values
        return None
    else:
        return metrics


def plot_qc_metrics(
    mdata,
    metrics: Optional[List[str]] = None,
    save: Optional[str] = None
):
    """
    Plot QC metrics distributions
    
    Parameters:
        mdata: MethylationData object
        metrics: List of metrics to plot. If None, plots all available
        save: Path to save figure
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    if metrics is None:
        metrics = ['n_sites', 'mean_coverage', 'mean_methylation']
    
    # Check metrics exist
    available_metrics = [m for m in metrics if m in mdata.obs.columns]
    
    if len(available_metrics) == 0:
        raise ValueError("No requested metrics found in mdata.obs")
    
    # Create figure
    n_metrics = len(available_metrics)
    fig, axes = plt.subplots(1, n_metrics, figsize=(5*n_metrics, 4))
    
    if n_metrics == 1:
        axes = [axes]
    
    # Plot each metric
    for ax, metric in zip(axes, available_metrics):
        values = mdata.obs[metric].values
        
        # Histogram
        ax.hist(values, bins=50, alpha=0.7, edgecolor='black')
        ax.set_xlabel(metric)
        ax.set_ylabel('Number of cells')
        ax.set_title(f'{metric} distribution')
        
        # Add median line
        median_val = np.median(values)
        ax.axvline(median_val, color='red', linestyle='--', 
                   label=f'Median: {median_val:.2f}')
        ax.legend()
    
    plt.tight_layout()
    
    if save:
        plt.savefig(save, dpi=300, bbox_inches='tight')
        print(f"Figure saved to {save}")
    
    return fig
