"""
Heatmap plotting
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Optional, List


def heatmap(
    mdata,
    var_names: Optional[List[str]] = None,
    groupby: Optional[str] = None,
    n_vars: int = 50,
    standard_scale: Optional[str] = 'var',
    cmap: str = 'RdBu_r',
    figsize: Optional[tuple] = None,
    save: Optional[str] = None,
    **kwargs
):
    """
    Methylation heatmap
    
    Parameters:
        mdata: MethylationData object
        var_names: Region names to plot
        groupby: Group cells by this column
        n_vars: Number of top variable regions if var_names=None
        standard_scale: 'var' (by column), 'obs' (by row), or None
        cmap: Colormap
        figsize: Figure size
        save: Save path
    """
    # Select regions
    if var_names:
        var_idx = [mdata.var.index.get_loc(v) for v in var_names if v in mdata.var.index]
    else:
        # Top variable regions
        variances = np.nanvar(mdata.X, axis=0)
        var_idx = np.argsort(variances)[-n_vars:]
    
    # Extract data
    plot_data = mdata.X[:, var_idx]
    plot_names = mdata.var.index[var_idx]
    
    # Convert to DataFrame
    df = pd.DataFrame(
        plot_data,
        index=mdata.obs_names,
        columns=plot_names
    )
    
    # Standardize
    if standard_scale == 'var':
        df = (df - df.mean()) / df.std()
    elif standard_scale == 'obs':
        df = df.sub(df.mean(axis=1), axis=0).div(df.std(axis=1), axis=0)
    
    # Sort by group
    if groupby and groupby in mdata.obs:
        df = df.loc[mdata.obs.sort_values(groupby).index]
        row_colors = mdata.obs.loc[df.index, groupby]
        
        if pd.api.types.is_categorical_dtype(row_colors):
            unique_groups = row_colors.cat.categories
            palette = sns.color_palette('tab20', len(unique_groups))
            lut = dict(zip(unique_groups, palette))
            row_colors = row_colors.map(lut)
    else:
        row_colors = None
    
    # Set size
    if figsize is None:
        figsize = (10, 8)
    
    # Plot
    g = sns.clustermap(
        df.T,
        cmap=cmap,
        center=0 if standard_scale else None,
        col_colors=row_colors,
        figsize=figsize,
        **kwargs
    )
    
    g.ax_heatmap.set_xlabel('Cells')
    g.ax_heatmap.set_ylabel('Regions')
    
    if save:
        plt.savefig(save, dpi=300, bbox_inches='tight')
        print(f"Saved to {save}")
    
    return g
