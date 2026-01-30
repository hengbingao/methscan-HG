"""
UMAP Dimensionality Reduction
"""

import numpy as np
from typing import Optional


def run_umap(
    mdata,
    min_dist: float = 0.5,
    n_neighbors: int = 15,
    n_components: int = 2,
    metric: str = 'euclidean',
    use_rep: str = 'X_pca',
    random_state: int = 0,
    inplace: bool = True
):
    """
    Run UMAP (Uniform Manifold Approximation and Projection)
    
    Parameters:
        mdata: MethylationData object
        min_dist: Minimum distance between points in embedding
        n_neighbors: Number of neighbors to consider
        n_components: Number of dimensions in embedding
        metric: Distance metric
        use_rep: Representation to use from mdata.obsm (usually 'X_pca')
        random_state: Random seed
        inplace: Store result in mdata.obsm['X_umap']
        
    Returns:
        None if inplace=True, UMAP embedding if inplace=False
    """
    try:
        from umap import UMAP
    except ImportError:
        raise ImportError("umap-learn is required. Install with: pip install umap-learn")
    
    # Get input data
    if use_rep in mdata.obsm:
        X = mdata.obsm[use_rep]
    elif use_rep == 'X':
        X = mdata.X
    else:
        raise ValueError(f"'{use_rep}' not found in mdata.obsm")
    
    print(f"Computing UMAP from {use_rep}...")
    
    # Run UMAP
    reducer = UMAP(
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        n_components=n_components,
        metric=metric,
        random_state=random_state
    )
    
    X_umap = reducer.fit_transform(X)
    
    if inplace:
        mdata.obsm['X_umap'] = X_umap
        mdata.uns['umap'] = {
            'params': {
                'n_neighbors': n_neighbors,
                'min_dist': min_dist,
                'metric': metric
            }
        }
        print(f"UMAP completed. Results stored in mdata.obsm['X_umap']")
        return None
    else:
        return X_umap
