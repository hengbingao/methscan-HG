#!/usr/bin/env python
"""
Example script showing how to use scMethylSeq programmatically

This script demonstrates the complete workflow from methylation matrices
to AnnData object with clustering.
"""

import numpy as np
import pandas as pd
from pathlib import Path

from scmethylseq import (
    PCAAnalyzer,
    UMAPAnalyzer,
    LeidenClustering,
    AnnDataBuilder,
)


def example_complete_workflow(
    cg_matrix_file: str,
    ch_matrix_file: str,
    output_dir: str,
):
    """
    Complete workflow example
    
    Parameters:
    -----------
    cg_matrix_file : str
        Path to CG methylation matrix
    ch_matrix_file : str
        Path to CH methylation matrix
    output_dir : str
        Output directory
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 70)
    print("scMethylSeq Complete Workflow Example")
    print("=" * 70)
    
    # ========================================================================
    # Step 1: Load methylation matrix
    # ========================================================================
    print("\n[Step 1] Loading CG methylation matrix...")
    
    # Load matrix (genes x cells format from your pipeline)
    cg_matrix = pd.read_csv(cg_matrix_file, sep='\t', header=None)
    
    # Transpose so rows are cells and columns are genes
    gene_names = cg_matrix.iloc[:, 0].values[1:]  # Skip header
    cell_names = cg_matrix.iloc[0, 1:].values
    cg_data = cg_matrix.iloc[1:, 1:].T
    cg_data.columns = gene_names
    cg_data.index = cell_names
    
    print(f"Matrix shape: {cg_data.shape}")
    print(f"Cells: {len(cell_names)}, Genes: {len(gene_names)}")
    
    # ========================================================================
    # Step 2: PCA Analysis
    # ========================================================================
    print("\n[Step 2] Running PCA...")
    
    pca_analyzer = PCAAnalyzer(n_components=20, center=True, scale=False)
    pca_coords, pca_df = pca_analyzer.fit_transform(cg_data, use_iterative=True)
    
    # Get variance explained
    var_df = pca_analyzer.get_explained_variance()
    var_df.to_csv(output_dir / "pca_variance_explained.csv", index=False)
    
    print(f"PCA complete. Explained variance (PC1-5): {var_df['variance_explained'].head().values}")
    
    # ========================================================================
    # Step 3: UMAP Embedding
    # ========================================================================
    print("\n[Step 3] Running UMAP...")
    
    # Use first 10 PCs for UMAP
    umap_analyzer = UMAPAnalyzer(n_neighbors=200, min_dist=0.1, random_state=42)
    umap_coords, umap_graph = umap_analyzer.fit_transform(pca_coords[:, :10])
    
    print(f"UMAP complete. Shape: {umap_coords.shape}")
    
    # ========================================================================
    # Step 4: Leiden Clustering
    # ========================================================================
    print("\n[Step 4] Running Leiden clustering...")
    
    # Scan multiple resolutions
    leiden_cluster = LeidenClustering(resolution=0.2)
    
    # Scan resolutions
    cluster_results = leiden_cluster.scan_resolutions(
        neighbor_graph=umap_graph,
        resolutions=[0.01, 0.05, 0.1, 0.2, 0.3, 0.5],
        cell_names=cell_names.tolist()
    )
    
    cluster_results.to_csv(output_dir / "clustering_multiple_resolutions.csv", index=False)
    
    # Use resolution 0.2 for final clustering
    final_clusters = cluster_results['leiden_res_0.2']
    
    print(f"Found {len(final_clusters.unique())} clusters at resolution 0.2")
    
    # ========================================================================
    # Step 5: Create metadata DataFrame
    # ========================================================================
    print("\n[Step 5] Creating metadata...")
    
    metadata = pd.DataFrame({
        'cell': cell_names,
        'leiden_cluster': final_clusters,
        'UMAP1': umap_coords[:, 0],
        'UMAP2': umap_coords[:, 1],
    })
    
    # Add PCA coordinates
    for i in range(20):
        metadata[f'PC{i+1}'] = pca_coords[:, i]
    
    metadata.to_csv(output_dir / "cell_metadata.csv", index=False)
    
    print(f"Metadata created for {len(metadata)} cells")
    
    # ========================================================================
    # Step 6: Build AnnData object
    # ========================================================================
    print("\n[Step 6] Building AnnData object...")
    
    adata_builder = AnnDataBuilder(
        cg_matrix_file=cg_matrix_file,
        ch_matrix_file=ch_matrix_file,
        metadata_file=str(output_dir / "cell_metadata.csv"),
        output_file=str(output_dir / "final_adata.h5ad")
    )
    
    adata = adata_builder.build_and_save()
    
    print(f"\nAnnData object created:")
    print(f"  Shape: {adata.shape}")
    print(f"  obs columns: {list(adata.obs.columns)}")
    print(f"  obsm keys: {list(adata.obsm.keys())}")
    
    # ========================================================================
    # Step 7: Visualization (optional)
    # ========================================================================
    print("\n[Step 7] Creating visualizations...")
    
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
        
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        
        # UMAP colored by cluster
        scatter = axes[0].scatter(
            metadata['UMAP1'],
            metadata['UMAP2'],
            c=metadata['leiden_cluster'].astype(int),
            cmap='tab20',
            s=1,
            alpha=0.7
        )
        axes[0].set_title('UMAP colored by Leiden cluster')
        axes[0].set_xlabel('UMAP1')
        axes[0].set_ylabel('UMAP2')
        axes[0].axis('equal')
        
        # PCA plot
        axes[1].scatter(
            metadata['PC1'],
            metadata['PC2'],
            c=metadata['leiden_cluster'].astype(int),
            cmap='tab20',
            s=1,
            alpha=0.7
        )
        axes[1].set_title('PCA (PC1 vs PC2)')
        axes[1].set_xlabel('PC1')
        axes[1].set_ylabel('PC2')
        axes[1].axis('equal')
        
        plt.tight_layout()
        plt.savefig(output_dir / "clustering_visualization.pdf", dpi=300)
        plt.close()
        
        print(f"Visualization saved to {output_dir / 'clustering_visualization.pdf'}")
        
    except ImportError:
        print("Matplotlib not available, skipping visualization")
    
    print("\n" + "=" * 70)
    print("Workflow complete!")
    print(f"Results saved to: {output_dir}")
    print("=" * 70)
    
    return adata


def example_load_and_analyze(h5ad_file: str):
    """
    Example of loading and analyzing existing AnnData
    
    Parameters:
    -----------
    h5ad_file : str
        Path to h5ad file
    """
    import scanpy as sc
    
    print("\n" + "=" * 70)
    print("Loading and analyzing existing AnnData")
    print("=" * 70)
    
    # Load AnnData
    adata = sc.read_h5ad(h5ad_file)
    print(f"\nLoaded AnnData: {adata}")
    
    # Basic statistics
    print("\n[Cell Statistics]")
    print(f"Total cells: {adata.n_obs}")
    print(f"Total genes: {adata.n_vars}")
    
    if 'leiden_cluster' in adata.obs.columns:
        print(f"\nCells per cluster:")
        print(adata.obs['leiden_cluster'].value_counts().sort_index())
    
    # Methylation statistics
    if 'X_cg' in adata.obsm:
        cg_data = adata.obsm['X_cg']
        print(f"\n[CG Methylation]")
        print(f"Mean: {np.nanmean(cg_data[cg_data != 2]):.3f}")
        print(f"Median: {np.nanmedian(cg_data[cg_data != 2]):.3f}")
    
    if 'X_ch' in adata.obsm:
        ch_data = adata.obsm['X_ch']
        print(f"\n[CH Methylation]")
        print(f"Mean: {np.nanmean(ch_data[ch_data != 2]):.3f}")
        print(f"Median: {np.nanmedian(ch_data[ch_data != 2]):.3f}")
    
    return adata


if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 4:
        print("Usage:")
        print("  python example_usage.py <cg_matrix> <ch_matrix> <output_dir>")
        print("\nOr to load existing h5ad:")
        print("  python example_usage.py --load <h5ad_file>")
        sys.exit(1)
    
    if sys.argv[1] == "--load":
        example_load_and_analyze(sys.argv[2])
    else:
        cg_matrix = sys.argv[1]
        ch_matrix = sys.argv[2]
        output_dir = sys.argv[3]
        
        adata = example_complete_workflow(cg_matrix, ch_matrix, output_dir)
