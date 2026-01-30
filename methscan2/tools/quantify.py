"""
Methylation Quantification
"""

import numpy as np
import pandas as pd
from typing import Optional, Union


def quantify_methylation(
    mdata,
    regions: Union[str, pd.DataFrame] = 'vmr',
    method: str = 'shrunken_residuals',
    shrinkage: float = 0.5,
    layer: Optional[str] = None,
    key_added: str = 'X_quantified',
    inplace: bool = True
):
    """
    Quantify methylation levels across genomic regions
    
    Parameters:
        mdata: MethylationData object
        regions: 'vmr' to use detected VMRs, or DataFrame with chrom/start/end
        method: Quantification method
            - 'shrunken_residuals': Shrunken mean of residuals (recommended)
            - 'mean': Simple mean methylation
            - 'median': Median methylation
        shrinkage: Shrinkage factor for residuals method (0-1)
        layer: Which layer to use ('rate', 'met', 'total'). None = X
        key_added: Key to store result in mdata.X or mdata.layers
        inplace: Store result in mdata
        
    Returns:
        Quantified matrix if inplace=False
    """
    # Get regions
    if regions == 'vmr':
        if 'vmr' not in mdata.uns or mdata.uns['vmr'] is None:
            raise ValueError("VMRs not detected. Run tl.detect_vmr() first.")
        region_df = mdata.uns['vmr']
    elif isinstance(regions, pd.DataFrame):
        region_df = regions
    else:
        raise ValueError("regions must be 'vmr' or a DataFrame")
    
    print(f"Quantifying methylation across {len(region_df)} regions...")
    
    # Get data layer
    if layer is None:
        data = mdata.methylation_rate
    elif layer in mdata.layers:
        data = mdata.layers[layer]
    else:
        raise ValueError(f"Layer '{layer}' not found")
    
    # Initialize result matrix
    n_cells = mdata.n_obs
    n_regions = len(region_df)
    result = np.zeros((n_cells, n_regions))
    
    # Process each region
    for region_idx, (_, region) in enumerate(region_df.iterrows()):
        # Find sites in this region
        site_mask = (
            (mdata.var['chrom'] == region['chrom']) &
            (mdata.var['start'] >= region['start']) &
            (mdata.var['end'] <= region['end'])
        )
        
        if site_mask.sum() == 0:
            continue
        
        # Extract region data
        region_data = data[:, site_mask]
        
        # Apply quantification method
        if method == 'mean':
            result[:, region_idx] = np.nanmean(region_data, axis=1)
            
        elif method == 'median':
            result[:, region_idx] = np.nanmedian(region_data, axis=1)
            
        elif method == 'shrunken_residuals':
            # Calculate residuals from smooth curve
            if 'smooth_methylation' not in mdata.uns:
                raise ValueError("Smooth curves not found. Run pp.smooth_methylation() first.")
            
            chrom = region['chrom']
            if chrom not in mdata.uns['smooth_methylation']:
                continue
            
            smooth_data = mdata.uns['smooth_methylation'][chrom]
            
            # Get smooth values for these sites
            positions = mdata.var.loc[site_mask, 'start'].values
            smooth_values = np.interp(
                positions,
                smooth_data['position'],
                smooth_data['methylation']
            )
            
            # Calculate residuals
            residuals = region_data - smooth_values[np.newaxis, :]
            
            # Apply shrinkage
            shrunken = residuals * shrinkage
            
            # Mean of shrunken residuals
            result[:, region_idx] = np.nanmean(shrunken, axis=1)
        
        else:
            raise ValueError(f"Unknown method: {method}")
    
    # Create region names
    region_names = [
        f"{r['chrom']}:{r['start']}-{r['end']}"
        for _, r in region_df.iterrows()
    ]
    
    if inplace:
        # Store in mdata
        if key_added == 'X':
            mdata.X = result
        else:
            if mdata.layers is None:
                mdata.layers = {}
            mdata.layers[key_added] = result
        
        # Update var to reflect new regions
        new_var = region_df.copy()
        new_var.index = region_names
        mdata.var = new_var
        
        print(f"Quantified data stored in mdata.{'X' if key_added == 'X' else 'layers[\'' + key_added + '\']'}")
        print(f"Data shape: {mdata.n_obs} cells × {mdata.n_vars} regions")
        
        return None
    else:
        return result
