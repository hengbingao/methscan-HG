"""
AnnData builder module
Construct AnnData objects from methylation matrices and metadata
"""

import logging
from pathlib import Path
from typing import Optional, Dict, List
import numpy as np
import pandas as pd
from anndata import AnnData

logger = logging.getLogger(__name__)


class AnnDataBuilder:
    """
    Build AnnData object from methylation data
    """
    
    def __init__(self,
                 cg_matrix_file: str,
                 ch_matrix_file: str,
                 metadata_file: str,
                 output_file: str):
        """
        Initialize AnnDataBuilder
        
        Parameters:
        -----------
        cg_matrix_file : str
            Path to CG methylation matrix (genes x cells)
        ch_matrix_file : str
            Path to CH methylation matrix (genes x cells)
        metadata_file : str
            Path to cell metadata CSV (with UMAP, clusters, etc.)
        output_file : str
            Output h5ad file path
        """
        self.cg_matrix_file = cg_matrix_file
        self.ch_matrix_file = ch_matrix_file
        self.metadata_file = metadata_file
        self.output_file = output_file
    
    def load_matrices(self) -> tuple:
        """
        Load CG and CH methylation matrices
        
        Returns:
        --------
        tuple
            (cg_df, ch_df, gene_names, cell_names)
        """
        logger.info("Loading methylation matrices...")
        
        # Load CG matrix
        cg_df = pd.read_csv(self.cg_matrix_file, sep='\t', header=None)
        
        # Load CH matrix
        ch_df = pd.read_csv(self.ch_matrix_file, sep='\t', header=None)
        
        logger.info(f"CG matrix shape: {cg_df.shape}")
        logger.info(f"CH matrix shape: {ch_df.shape}")
        
        # Transpose so rows are cells and columns are genes
        # First row contains cell names, first column contains gene names
        gene_names = cg_df.iloc[:, 0].values
        
        # Get cell names from first row, excluding first element (which is empty or gene label)
        cell_names_cg = cg_df.iloc[0, 1:].values
        cell_names_ch = ch_df.iloc[0, 1:].values
        
        # Get data (excluding first row and first column)
        cg_data = cg_df.iloc[1:, 1:].T  # Transpose so cells are rows
        ch_data = ch_df.iloc[1:, 1:].T
        
        # Remove first gene name (header)
        gene_names = gene_names[1:]
        
        logger.info(f"Genes: {len(gene_names)}, Cells (CG): {len(cell_names_cg)}, Cells (CH): {len(cell_names_ch)}")
        
        return cg_data, ch_data, gene_names, cell_names_cg, cell_names_ch
    
    def load_metadata(self, cell_names: np.ndarray) -> pd.DataFrame:
        """
        Load and process cell metadata
        
        Parameters:
        -----------
        cell_names : np.ndarray
            Cell names from matrix
            
        Returns:
        --------
        pd.DataFrame
            Cell metadata
        """
        logger.info("Loading cell metadata...")
        
        metadata = pd.read_csv(self.metadata_file)
        
        # Ensure cell column exists
        if 'cell' not in metadata.columns:
            raise ValueError("Metadata must contain 'cell' column")
        
        # Filter to common cells
        metadata = metadata[metadata['cell'].astype(str).isin(cell_names.astype(str))]
        
        # Sort to match cell_names order
        metadata = metadata.set_index('cell')
        metadata = metadata.loc[cell_names]
        metadata = metadata.reset_index()
        
        logger.info(f"Loaded metadata for {len(metadata)} cells")
        
        return metadata
    
    def build_adata(self,
                   store_ch_in_layers: bool = False) -> AnnData:
        """
        Build AnnData object
        
        Parameters:
        -----------
        store_ch_in_layers : bool
            If True, store CH in layers. Otherwise, store in obsm.
            
        Returns:
        --------
        AnnData
            Complete AnnData object
        """
        logger.info("Building AnnData object...")
        
        # Load data
        cg_data, ch_data, gene_names, cell_names_cg, cell_names_ch = self.load_matrices()
        
        # Verify cell names match
        if not np.array_equal(cell_names_cg, cell_names_ch):
            logger.warning("CG and CH cell names don't match exactly, finding intersection...")
            common_cells = np.intersect1d(cell_names_cg, cell_names_ch)
            logger.info(f"Using {len(common_cells)} common cells")
            
            # Filter data
            cg_mask = np.isin(cell_names_cg, common_cells)
            ch_mask = np.isin(cell_names_ch, common_cells)
            
            cg_data = cg_data.iloc[cg_mask, :]
            ch_data = ch_data.iloc[ch_mask, :]
            cell_names = common_cells
        else:
            cell_names = cell_names_cg
        
        # Load metadata
        metadata = self.load_metadata(cell_names)
        
        # Further filter to cells in metadata
        common_cells = np.intersect1d(cell_names, metadata['cell'].values)
        logger.info(f"Final cell count after metadata merge: {len(common_cells)}")
        
        cg_mask = np.isin(cell_names, common_cells)
        ch_mask = np.isin(cell_names, common_cells)
        meta_mask = metadata['cell'].isin(common_cells)
        
        cg_data = cg_data.iloc[cg_mask, :]
        ch_data = ch_data.iloc[ch_mask, :]
        metadata = metadata[meta_mask]
        cell_names = common_cells
        
        # Convert to float and handle NA
        logger.info("Converting to float and handling missing values...")
        cg_matrix = cg_data.values.astype(float)
        ch_matrix = ch_data.values.astype(float)
        
        # Replace NA with 2 (as in original code)
        cg_matrix[np.isnan(cg_matrix)] = 2
        ch_matrix[np.isnan(ch_matrix)] = 2
        
        # Create AnnData object
        logger.info("Creating AnnData object...")
        adata = AnnData(
            X=cg_matrix,
            var=pd.DataFrame(index=gene_names),
            obs=pd.DataFrame(index=cell_names)
        )
        
        # Add metadata to obs
        logger.info("Adding metadata to obs...")
        for col in metadata.columns:
            if col != 'cell':
                adata.obs[col] = metadata[col].values
        
        # Store UMAP coordinates in obsm if available
        if 'UMAP1' in metadata.columns and 'UMAP2' in metadata.columns:
            logger.info("Adding UMAP coordinates to obsm...")
            umap_coords = metadata[['UMAP1', 'UMAP2']].values
            adata.obsm['X_umap'] = umap_coords
        
        # Store CG and CH matrices in obsm
        logger.info("Storing methylation matrices in obsm...")
        adata.obsm['X_cg'] = cg_matrix
        adata.obsm['X_ch'] = ch_matrix
        
        # Optionally store CH in layers
        if store_ch_in_layers:
            adata.layers['mCH'] = ch_matrix
        
        logger.info(f"AnnData object created: {adata}")
        logger.info(f"  Cells: {adata.n_obs}")
        logger.info(f"  Genes: {adata.n_vars}")
        logger.info(f"  obs columns: {list(adata.obs.columns)}")
        logger.info(f"  obsm keys: {list(adata.obsm.keys())}")
        
        return adata
    
    def save_adata(self, adata: AnnData) -> None:
        """
        Save AnnData object to file
        
        Parameters:
        -----------
        adata : AnnData
            AnnData object to save
        """
        logger.info(f"Saving AnnData to {self.output_file}...")
        adata.write_h5ad(self.output_file)
        logger.info("AnnData saved successfully")
    
    def build_and_save(self) -> AnnData:
        """
        Build and save AnnData object
        
        Returns:
        --------
        AnnData
            Built AnnData object
        """
        adata = self.build_adata()
        self.save_adata(adata)
        return adata


def add_pca_to_adata(adata: AnnData,
                     pca_coords: np.ndarray,
                     variance_explained: np.ndarray) -> AnnData:
    """
    Add PCA coordinates to existing AnnData object
    
    Parameters:
    -----------
    adata : AnnData
        AnnData object
    pca_coords : np.ndarray
        PCA coordinates (cells x components)
    variance_explained : np.ndarray
        Variance explained by each component
        
    Returns:
    --------
    AnnData
        Updated AnnData object
    """
    logger.info("Adding PCA to AnnData...")
    
    adata.obsm['X_pca'] = pca_coords
    adata.uns['pca'] = {
        'variance': variance_explained,
        'variance_ratio': variance_explained / variance_explained.sum()
    }
    
    return adata


def add_clustering_to_adata(adata: AnnData,
                            clusters: pd.Series,
                            key: str = 'leiden') -> AnnData:
    """
    Add clustering results to AnnData object
    
    Parameters:
    -----------
    adata : AnnData
        AnnData object
    clusters : pd.Series
        Cluster assignments
    key : str
        Key for storing clusters in obs
        
    Returns:
    --------
    AnnData
        Updated AnnData object
    """
    logger.info(f"Adding {key} clustering to AnnData...")
    
    adata.obs[key] = clusters.astype('category')
    
    return adata


def integrate_scanpy_workflow(adata: AnnData,
                              n_pcs: int = 20,
                              n_neighbors: int = 200,
                              resolution: float = 0.2) -> AnnData:
    """
    Run standard scanpy workflow on AnnData
    
    Parameters:
    -----------
    adata : AnnData
        Input AnnData object
    n_pcs : int
        Number of PCs for analysis
    n_neighbors : int
        Number of neighbors
    resolution : float
        Clustering resolution
        
    Returns:
    --------
    AnnData
        Processed AnnData object
    """
    try:
        import scanpy as sc
    except ImportError:
        raise ImportError("scanpy required. Install with: pip install scanpy")
    
    logger.info("Running scanpy workflow...")
    
    # Handle missing values for scanpy
    if np.isnan(adata.X).any():
        logger.info("Imputing missing values for scanpy...")
        adata.X[np.isnan(adata.X)] = 0
    
    # PCA
    logger.info("Computing PCA...")
    sc.tl.pca(adata, n_comps=n_pcs, svd_solver='arpack')
    
    # Neighborhood graph
    logger.info("Computing neighborhood graph...")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    
    # UMAP (if not already present)
    if 'X_umap' not in adata.obsm:
        logger.info("Computing UMAP...")
        sc.tl.umap(adata)
    
    # Leiden clustering
    logger.info("Running Leiden clustering...")
    sc.tl.leiden(adata, resolution=resolution, key_added='leiden')
    
    logger.info("Scanpy workflow complete")
    
    return adata


def validate_adata(adata: AnnData) -> Dict[str, bool]:
    """
    Validate AnnData object structure
    
    Parameters:
    -----------
    adata : AnnData
        AnnData object to validate
        
    Returns:
    --------
    Dict[str, bool]
        Validation results
    """
    results = {
        'has_X': adata.X is not None,
        'has_obs': len(adata.obs) > 0,
        'has_var': len(adata.var) > 0,
        'has_umap': 'X_umap' in adata.obsm,
        'has_pca': 'X_pca' in adata.obsm,
        'has_cg': 'X_cg' in adata.obsm,
        'has_ch': 'X_ch' in adata.obsm,
        'has_clusters': 'leiden' in adata.obs.columns or 'leiden_cluster' in adata.obs.columns,
    }
    
    logger.info("AnnData validation results:")
    for key, value in results.items():
        logger.info(f"  {key}: {value}")
    
    return results
