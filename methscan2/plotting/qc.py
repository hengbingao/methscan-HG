"""
QC Plotting Functions
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Optional, List


def qc_violin(
    mdata,
    keys: List[str] = ['n_sites', 'mean_coverage', 'mean_methylation'],
    groupby: Optional[str] = None,
    figsize: Optional[tuple] = None,
    save: Optional[str] = None,
    **kwargs
):
    """
    QC metrics violin plot
    
    Parameters:
        mdata: MethylationData object
        keys: Metrics to plot from mdata.obs
        groupby: Column to group by
        figsize: Figure size
        save: Save path
        **kwargs: Additional arguments for seaborn.violinplot
    """
    available = [k for k in keys if k in mdata.obs.columns]
    
    if len(available) == 0:
        raise ValueError(f"None of {keys} found in mdata.obs")
    
    if figsize is None:
        figsize = (len(available) * 4, 4)
    
    fig, axes = plt.subplots(1, len(available), figsize=figsize)
    if len(available) == 1:
        axes = [axes]
    
    for ax, key in zip(axes, available):
        if groupby and groupby in mdata.obs:
            sns.violinplot(
                data=mdata.obs,
                x=groupby,
                y=key,
                ax=ax,
                **kwargs
            )
            ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
        else:
            sns.violinplot(data=mdata.obs, y=key, ax=ax, **kwargs)
        
        ax.set_ylabel(key)
        ax.set_xlabel('')
    
    plt.tight_layout()
    
    if save:
        plt.savefig(save, dpi=300, bbox_inches='tight')
        print(f"Saved to {save}")
    
    return fig


def qc_scatter(
    mdata,
    x: str = 'n_sites',
    y: str = 'mean_methylation',
    color: Optional[str] = None,
    figsize: tuple = (6, 5),
    save: Optional[str] = None,
    **kwargs
):
    """
    QC scatter plot
    
    Parameters:
        mdata: MethylationData object
        x: Column for x-axis
        y: Column for y-axis
        color: Column for coloring points
        figsize: Figure size
        save: Save path
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    if color and color in mdata.obs:
        scatter = ax.scatter(
            mdata.obs[x],
            mdata.obs[y],
            c=mdata.obs[color],
            alpha=0.6,
            **kwargs
        )
        plt.colorbar(scatter, ax=ax, label=color)
    else:
        ax.scatter(
            mdata.obs[x],
            mdata.obs[y],
            alpha=0.6,
            **kwargs
        )
    
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    
    plt.tight_layout()
    
    if save:
        plt.savefig(save, dpi=300, bbox_inches='tight')
        print(f"Saved to {save}")
    
    return fig
