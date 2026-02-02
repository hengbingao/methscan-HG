"""
Clustering Functions
"""

import numpy as np
from typing import Optional


def run_leiden(
    mdata,
    resolution: float = 1.0,
    n_neighbors: int = 15,
    use_rep: str = 'X_pca',
    key_added: str = 'leiden',
    random_state: int = 0,
    inplace: bool = True
):
    """
    Leiden clustering
    
    Parameters:
        mdata: MethylationData object
        resolution: Resolution parameter (higher = more clusters)
        n_neighbors: Number of neighbors for graph construction
        use_rep: Representation to use
        key_added: Key for storing cluster labels in mdata.obs
        random_state: Random seed
        inplace: Store results in mdata.obs
        
    Returns:
        None if inplace=True, cluster labels if inplace=False
    """
    try:
        import igraph as ig
        import leidenalg
    except ImportError:
        raise ImportError("leidenalg and python-igraph required. "
                         "Install with: pip install leidenalg python-igraph")
    
    # Get representation
    if use_rep in mdata.obsm:
        X = mdata.obsm[use_rep]
    else:
        raise ValueError(f"'{use_rep}' not found")
    
    print(f"Running Leiden clustering (resolution={resolution})...")
    
    # Build KNN graph
    from sklearn.neighbors import NearestNeighbors
    
    nn = NearestNeighbors(n_neighbors=n_neighbors, metric='euclidean')
    nn.fit(X)
    distances, indices = nn.kneighbors(X)
    
    # Build igraph
    edges = []
    weights = []
    
    for i in range(len(indices)):
        for j, neighbor in enumerate(indices[i]):
            if i != neighbor:
                edges.append((i, neighbor))
                weights.append(1.0 / (1.0 + distances[i, j]))
    
    g = ig.Graph(n=len(X), edges=edges, directed=False)
    g.es['weight'] = weights
    
    # Run Leiden
    partition = leidenalg.find_partition(
        g,
        leidenalg.RBConfigurationVertexPartition,
        weights='weight',
        resolution_parameter=resolution,
        seed=random_state
    )
    
    # Get cluster labels
    clusters = np.array(partition.membership)
    
    print(f"Found {len(np.unique(clusters))} clusters")
    
    if inplace:
        mdata.obs[key_added] = clusters
        mdata.obs[key_added] = mdata.obs[key_added].astype('category')
        print(f"Cluster labels stored in mdata.obs['{key_added}']")
        return None
    else:
        return clusters


def run_louvain(
    mdata,
    resolution: float = 1.0,
    n_neighbors: int = 15,
    use_rep: str = 'X_pca',
    key_added: str = 'louvain',
    random_state: int = 0,
    inplace: bool = True
):
    """
    Louvain clustering
    
    Parameters similar to run_leiden
    """
    try:
        import igraph as ig
    except ImportError:
        raise ImportError("python-igraph required. Install with: pip install python-igraph")
    
    # Get representation
    if use_rep in mdata.obsm:
        X = mdata.obsm[use_rep]
    else:
        raise ValueError(f"'{use_rep}' not found")
    
    print(f"Running Louvain clustering (resolution={resolution})...")
    
    # Build KNN graph (same as Leiden)
    from sklearn.neighbors import NearestNeighbors
    
    nn = NearestNeighbors(n_neighbors=n_neighbors, metric='euclidean')
    nn.fit(X)
    distances, indices = nn.kneighbors(X)
    
    edges = []
    weights = []
    
    for i in range(len(indices)):
        for j, neighbor in enumerate(indices[i]):
            if i != neighbor:
                edges.append((i, neighbor))
                weights.append(1.0 / (1.0 + distances[i, j]))
    
    g = ig.Graph(n=len(X), edges=edges, directed=False)
    g.es['weight'] = weights
    
    # Run Louvain
    partition = g.community_multilevel(weights='weight', return_levels=False)
    
    clusters = np.array(partition.membership)
    
    print(f"Found {len(np.unique(clusters))} clusters")
    
    if inplace:
        mdata.obs[key_added] = clusters
        mdata.obs[key_added] = mdata.obs[key_added].astype('category')
        print(f"Cluster labels stored in mdata.obs['{key_added}']")
        return None
    else:
        return clusters
