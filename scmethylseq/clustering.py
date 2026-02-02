"""
Clustering module
PCA, UMAP, and Leiden clustering for methylation data
"""

import logging
from typing import Optional, List, Tuple, Dict
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.decomposition import PCA as skPCA
from sklearn.preprocessing import StandardScaler
import warnings

logger = logging.getLogger(__name__)


class IterativePCA:
    """
    Iterative PCA for handling missing values in methylation data
    """
    
    def __init__(self, 
                 n_components: int = 20,
                 max_iter: int = 50,
                 min_gain: float = 0.001,
                 random_state: int = 42):
        """
        Initialize IterativePCA
        
        Parameters:
        -----------
        n_components : int
            Number of principal components
        max_iter : int
            Maximum number of iterations
        min_gain : float
            Minimum improvement threshold
        random_state : int
            Random state for reproducibility
        """
        self.n_components = n_components
        self.max_iter = max_iter
        self.min_gain = min_gain
        self.random_state = random_state
        
        self.components_ = None
        self.explained_variance_ = None
        self.explained_variance_ratio_ = None
        self.mean_ = None
        self.mse_iter_ = []
    
    def fit_transform(self, X: np.ndarray) -> np.ndarray:
        """
        Fit the model and transform the data
        
        Parameters:
        -----------
        X : np.ndarray
            Data matrix (cells x features) with NaN for missing values
            
        Returns:
        --------
        np.ndarray
            Transformed data (cells x components)
        """
        logger.info(f"Starting iterative PCA with {self.n_components} components")
        
        # Make a copy to avoid modifying original data
        X_imputed = X.copy()
        
        # Identify missing value locations
        na_mask = np.isnan(X_imputed)
        n_missing = na_mask.sum()
        
        logger.info(f"Found {n_missing} missing values ({100*n_missing/X.size:.2f}%)")
        
        # Initial imputation with 0
        X_imputed[na_mask] = 0
        
        # Center the data
        self.mean_ = np.nanmean(X, axis=0)
        X_centered = X_imputed - self.mean_
        
        # Iterative imputation
        mse_values = []
        for iteration in range(self.max_iter):
            # Store previous imputed values
            prev_imputed = X_centered[na_mask].copy()
            
            # Perform PCA
            pca = skPCA(n_components=self.n_components, random_state=self.random_state)
            X_transformed = pca.fit_transform(X_centered)
            
            # Reconstruct data
            X_reconstructed = X_transformed @ pca.components_
            
            # Update missing values
            X_centered[na_mask] = X_reconstructed[na_mask]
            
            # Calculate MSE
            mse = np.mean((prev_imputed - X_centered[na_mask]) ** 2)
            mse_values.append(mse)
            
            # Check convergence
            if iteration > 0:
                gain = mse / max(mse_values)
                if gain < self.min_gain:
                    logger.info(f"Converged after {iteration + 1} iterations")
                    break
            
            if (iteration + 1) % 10 == 0:
                logger.info(f"Iteration {iteration + 1}/{self.max_iter}, MSE: {mse:.6f}")
        
        # Store final PCA results
        self.components_ = pca.components_
        self.explained_variance_ = pca.explained_variance_
        self.explained_variance_ratio_ = pca.explained_variance_ratio_
        self.mse_iter_ = mse_values
        
        logger.info(f"Final MSE: {mse_values[-1]:.6f}")
        logger.info(f"Explained variance: {self.explained_variance_ratio_[:5]}")
        
        return X_transformed
    
    def transform(self, X: np.ndarray) -> np.ndarray:
        """
        Transform new data using fitted PCA
        
        Parameters:
        -----------
        X : np.ndarray
            Data matrix to transform
            
        Returns:
        --------
        np.ndarray
            Transformed data
        """
        if self.components_ is None:
            raise ValueError("PCA has not been fitted yet")
        
        # Handle missing values by imputing with mean
        X_imputed = X.copy()
        na_mask = np.isnan(X_imputed)
        X_imputed[na_mask] = 0
        
        # Center and transform
        X_centered = X_imputed - self.mean_
        return X_centered @ self.components_.T


class PCAAnalyzer:
    """
    PCA analysis for methylation data
    """
    
    def __init__(self, 
                 n_components: int = 20,
                 center: bool = True,
                 scale: bool = False):
        """
        Initialize PCAAnalyzer
        
        Parameters:
        -----------
        n_components : int
            Number of principal components
        center : bool
            Whether to center the data
        scale : bool
            Whether to scale the data
        """
        self.n_components = n_components
        self.center = center
        self.scale = scale
        self.pca = None
        self.scaler = None
    
    def fit_transform(self, 
                     matrix: pd.DataFrame,
                     use_iterative: bool = True) -> Tuple[np.ndarray, pd.DataFrame]:
        """
        Perform PCA on methylation matrix
        
        Parameters:
        -----------
        matrix : pd.DataFrame
            Methylation matrix (cells x features)
        use_iterative : bool
            Whether to use iterative PCA for missing values
            
        Returns:
        --------
        Tuple[np.ndarray, pd.DataFrame]
            PCA transformed data and metadata DataFrame
        """
        logger.info(f"Running PCA on matrix of shape {matrix.shape}")
        
        # Get cell names
        if isinstance(matrix.index, pd.RangeIndex):
            cell_names = [f"cell_{i}" for i in range(len(matrix))]
        else:
            cell_names = matrix.index.tolist()
        
        # Convert to numpy array
        X = matrix.values
        
        # Handle missing values
        if np.isnan(X).any():
            if use_iterative:
                logger.info("Using iterative PCA for missing value imputation")
                self.pca = IterativePCA(
                    n_components=self.n_components,
                    random_state=42
                )
                X_pca = self.pca.fit_transform(X)
            else:
                logger.warning("Replacing NaN with 0 for PCA")
                X[np.isnan(X)] = 0
                
                # Scale if requested
                if self.scale or self.center:
                    self.scaler = StandardScaler(
                        with_mean=self.center,
                        with_std=self.scale
                    )
                    X = self.scaler.fit_transform(X)
                
                self.pca = skPCA(n_components=self.n_components, random_state=42)
                X_pca = self.pca.fit_transform(X)
        else:
            # Scale if requested
            if self.scale or self.center:
                self.scaler = StandardScaler(
                    with_mean=self.center,
                    with_std=self.scale
                )
                X = self.scaler.fit_transform(X)
            
            self.pca = skPCA(n_components=self.n_components, random_state=42)
            X_pca = self.pca.fit_transform(X)
        
        # Create result DataFrame
        pc_cols = [f"PC{i+1}" for i in range(self.n_components)]
        pca_df = pd.DataFrame(X_pca, columns=pc_cols, index=cell_names)
        pca_df['cell'] = cell_names
        
        logger.info(f"PCA complete. Explained variance: {self.pca.explained_variance_ratio_[:5]}")
        
        return X_pca, pca_df
    
    def get_explained_variance(self) -> pd.DataFrame:
        """
        Get explained variance for each component
        
        Returns:
        --------
        pd.DataFrame
            DataFrame with variance explained by each PC
        """
        if self.pca is None:
            raise ValueError("PCA has not been fitted yet")
        
        var_exp = pd.DataFrame({
            'PC': range(1, len(self.pca.explained_variance_ratio_) + 1),
            'variance_explained': self.pca.explained_variance_ratio_,
            'cumulative_variance': np.cumsum(self.pca.explained_variance_ratio_)
        })
        
        return var_exp


class UMAPAnalyzer:
    """
    UMAP dimensionality reduction for visualization
    """
    
    def __init__(self,
                 n_neighbors: int = 200,
                 min_dist: float = 0.1,
                 n_components: int = 2,
                 metric: str = 'euclidean',
                 random_state: int = 42):
        """
        Initialize UMAPAnalyzer
        
        Parameters:
        -----------
        n_neighbors : int
            Number of neighbors for UMAP
        min_dist : float
            Minimum distance parameter
        n_components : int
            Number of UMAP dimensions
        metric : str
            Distance metric
        random_state : int
            Random state for reproducibility
        """
        self.n_neighbors = n_neighbors
        self.min_dist = min_dist
        self.n_components = n_components
        self.metric = metric
        self.random_state = random_state
        self.umap_obj = None
    
    def fit_transform(self, X: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        Fit UMAP and transform data
        
        Parameters:
        -----------
        X : np.ndarray
            Input data (usually PCA coordinates)
            
        Returns:
        --------
        Tuple[np.ndarray, np.ndarray]
            UMAP coordinates and neighbor indices
        """
        try:
            import umap
        except ImportError:
            raise ImportError("umap-learn is required. Install with: pip install umap-learn")
        
        logger.info(f"Running UMAP with n_neighbors={self.n_neighbors}")
        
        self.umap_obj = umap.UMAP(
            n_neighbors=self.n_neighbors,
            min_dist=self.min_dist,
            n_components=self.n_components,
            metric=self.metric,
            random_state=self.random_state,
            verbose=True
        )
        
        # Fit and transform
        embedding = self.umap_obj.fit_transform(X)
        
        logger.info(f"UMAP complete. Shape: {embedding.shape}")
        
        return embedding, self.umap_obj.graph_
    
    def scan_parameters(self,
                       X: np.ndarray,
                       n_neighbors_list: List[int] = [10, 20, 30, 50, 100, 200, 300, 400]
                       ) -> Dict[int, np.ndarray]:
        """
        Scan different n_neighbors parameters
        
        Parameters:
        -----------
        X : np.ndarray
            Input data
        n_neighbors_list : List[int]
            List of n_neighbors values to try
            
        Returns:
        --------
        Dict[int, np.ndarray]
            Dictionary mapping n_neighbors to UMAP embeddings
        """
        results = {}
        
        for n_neighbors in n_neighbors_list:
            logger.info(f"Testing n_neighbors={n_neighbors}")
            
            temp_analyzer = UMAPAnalyzer(
                n_neighbors=n_neighbors,
                min_dist=self.min_dist,
                random_state=self.random_state
            )
            
            embedding, _ = temp_analyzer.fit_transform(X)
            results[n_neighbors] = embedding
        
        return results


class LeidenClustering:
    """
    Leiden clustering on neighborhood graph
    """
    
    def __init__(self, resolution: float = 0.2):
        """
        Initialize LeidenClustering
        
        Parameters:
        -----------
        resolution : float
            Resolution parameter for Leiden algorithm
        """
        self.resolution = resolution
        self.clusters = None
    
    def cluster(self, 
                neighbor_graph,
                cell_names: Optional[List[str]] = None) -> pd.DataFrame:
        """
        Perform Leiden clustering
        
        Parameters:
        -----------
        neighbor_graph : scipy.sparse matrix or igraph Graph
            Neighbor graph
        cell_names : List[str]
            Cell names
            
        Returns:
        --------
        pd.DataFrame
            DataFrame with cluster assignments
        """
        try:
            import igraph as ig
            import leidenalg
        except ImportError:
            raise ImportError("igraph and leidenalg required. Install with: pip install igraph leidenalg")
        
        logger.info(f"Running Leiden clustering with resolution={self.resolution}")
        
        # Convert neighbor graph to igraph if needed
        if not isinstance(neighbor_graph, ig.Graph):
            # Assume it's a scipy sparse matrix
            sources, targets = neighbor_graph.nonzero()
            edges = list(zip(sources.tolist(), targets.tolist()))
            
            g = ig.Graph(n=neighbor_graph.shape[0])
            g.add_edges(edges)
            
            # Add weights if available
            if hasattr(neighbor_graph, 'data'):
                g.es['weight'] = neighbor_graph.data.tolist()
        else:
            g = neighbor_graph
        
        # Run Leiden clustering
        partition = leidenalg.find_partition(
            g,
            leidenalg.RBConfigurationVertexPartition,
            resolution_parameter=self.resolution,
            seed=42
        )
        
        self.clusters = np.array(partition.membership)
        
        logger.info(f"Found {len(set(self.clusters))} clusters")
        
        # Create result DataFrame
        if cell_names is None:
            cell_names = [f"cell_{i}" for i in range(len(self.clusters))]
        
        cluster_df = pd.DataFrame({
            'cell': cell_names,
            'leiden_cluster': self.clusters.astype(str)
        })
        
        return cluster_df
    
    def scan_resolutions(self,
                        neighbor_graph,
                        resolutions: List[float] = [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5],
                        cell_names: Optional[List[str]] = None) -> pd.DataFrame:
        """
        Scan multiple resolution parameters
        
        Parameters:
        -----------
        neighbor_graph : igraph Graph or scipy sparse matrix
            Neighbor graph
        resolutions : List[float]
            List of resolution values to try
        cell_names : List[str]
            Cell names
            
        Returns:
        --------
        pd.DataFrame
            DataFrame with cluster assignments for each resolution
        """
        if cell_names is None:
            cell_names = [f"cell_{i}" for i in range(neighbor_graph.shape[0])]
        
        results = pd.DataFrame({'cell': cell_names})
        
        for res in resolutions:
            logger.info(f"Testing resolution={res}")
            
            temp_cluster = LeidenClustering(resolution=res)
            cluster_df = temp_cluster.cluster(neighbor_graph, cell_names)
            
            results[f'leiden_res_{res}'] = cluster_df['leiden_cluster']
        
        return results


def build_neighbor_graph_from_umap(umap_obj, 
                                    n_neighbors: int,
                                    cell_names: List[str]) -> pd.DataFrame:
    """
    Build neighbor graph edge list from UMAP object
    
    Parameters:
    -----------
    umap_obj : umap.UMAP
        Fitted UMAP object
    n_neighbors : int
        Number of neighbors
    cell_names : List[str]
        Cell names
        
    Returns:
    --------
    pd.DataFrame
        Edge list with columns: from, to, weight
    """
    # Get neighbor indices and distances
    knn_indices = umap_obj.graph_.tocoo()
    
    edges = []
    for i, j, dist in zip(knn_indices.row, knn_indices.col, knn_indices.data):
        if i != j:  # Exclude self-loops
            edges.append({
                'from': cell_names[i],
                'to': cell_names[j],
                'weight': dist
            })
    
    edge_df = pd.DataFrame(edges)
    
    logger.info(f"Built neighbor graph with {len(edge_df)} edges")
    
    return edge_df
