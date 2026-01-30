"""
PCA Dimensionality Reduction
"""

import numpy as np
from sklearn.decomposition import PCA
from typing import Optional


def run_pca(
    mdata,
    n_comps: int = 50,
    zero_center: bool = True,
    svd_solver: str = 'auto',
    random_state: int = 0,
    use_highly_variable: bool = False,
    inplace: bool = True
):
    """
    Run Principal Component Analysis (PCA)
    
    Parameters:
        mdata: MethylationData object
        n_comps: Number of principal components
        zero_center: Zero-center the data before PCA
        svd_solver: SVD solver to use ('auto', 'full', 'arpack', 'randomized')
        random_state: Random seed
        use_highly_variable: Use only highly variable regions (if computed)
        inplace: Store results in mdata.obsm['X_pca']
        
    Returns:
        None if inplace=True, PCA result array if inplace=False
    """
    # Get data
    if use_highly_variable and 'highly_variable' in mdata.var.columns:
        X = mdata.X[:, mdata.var['highly_variable'].values]
        print(f"Using {X.shape[1]} highly variable regions")
    else:
        X = mdata.X
    
    # Handle missing values
    X = np.nan_to_num(X, nan=0.0)
    
    print(f"Computing PCA with {n_comps} components...")
    
    # Run PCA
    pca = PCA(
        n_components=n_comps,
        svd_solver=svd_solver,
        random_state=random_state
    )
    
    if zero_center:
        X_centered = X - X.mean(axis=0)
    else:
        X_centered = X
    
    X_pca = pca.fit_transform(X_centered)
    
    if inplace:
        # Store in obsm
        mdata.obsm['X_pca'] = X_pca
        
        # Store variance explained
        mdata.uns['pca'] = {
            'variance': pca.explained_variance_,
            'variance_ratio': pca.explained_variance_ratio_,
            'components': pca.components_
        }
        
        print(f"PCA completed. Results stored in mdata.obsm['X_pca']")
        print(f"Variance explained by first 10 PCs: {pca.explained_variance_ratio_[:10].sum():.2%}")
        
        return None
    else:
        return X_pca
