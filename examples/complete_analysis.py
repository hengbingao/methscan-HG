#!/usr/bin/env python3
"""
MethSCAn2 Complete Analysis Example

This script demonstrates a complete single-cell methylation analysis workflow.
"""

import methscan2 as ms2
import os

def main():
    """Run complete analysis pipeline"""
    
    # Create output directory
    os.makedirs('results', exist_ok=True)
    
    print("=" * 60)
    print("MethSCAn2 - Complete Analysis Pipeline")
    print("=" * 60)
    
    # ====================================================================
    # STEP 1: Load Data
    # ====================================================================
    print("\n[1/8] Loading methylation data...")
    
    # Read .cov files from a directory
    # Format: chr start end meth% met_count unmet_count
    mdata = ms2.read_cov_files(
        'data/*.cov',  # Path pattern
        genome='mm10',  # Genome assembly
        min_coverage=1,  # Minimum coverage
        n_jobs=-1  # Use all CPUs
    )
    
    print(f"Loaded: {mdata.n_obs} cells × {mdata.n_vars} sites")
    
    # ====================================================================
    # STEP 2: Quality Control
    # ====================================================================
    print("\n[2/8] Running quality control...")
    
    # Calculate QC metrics
    ms2.pp.calculate_qc_metrics(mdata)
    
    # Plot QC distributions
    ms2.pl.qc_violin(
        mdata,
        keys=['n_sites', 'mean_coverage', 'mean_methylation'],
        save='results/qc_metrics.pdf'
    )
    
    # Filter low-quality cells
    print("Filtering cells...")
    ms2.pp.filter_cells(
        mdata,
        min_sites=50000,      # At least 50k CpG sites
        min_coverage=5,       # Mean coverage ≥ 5
        min_methylation=0.2,  # Global methylation ≥ 20%
        max_methylation=0.8   # Global methylation ≤ 80%
    )
    
    # Filter sites
    print("Filtering sites...")
    ms2.pp.filter_sites(
        mdata,
        min_cells=10,         # Present in at least 10 cells
        min_coverage=5        # With coverage ≥ 5
    )
    
    print(f"After filtering: {mdata.n_obs} cells × {mdata.n_vars} sites")
    
    # ====================================================================
    # STEP 3: VMR Detection
    # ====================================================================
    print("\n[3/8] Detecting variably methylated regions (VMRs)...")
    
    # First, calculate smoothed methylation curves
    ms2.pp.smooth_methylation(
        mdata,
        bandwidth=2000,  # 2kb smoothing window
        n_jobs=-1
    )
    
    # Detect VMRs
    ms2.tl.detect_vmr(
        mdata,
        bandwidth=2000,
        stepsize=100,
        var_threshold=0.02,  # Variance threshold
        min_coverage=5,
        min_cells=10,
        n_jobs=-1
    )
    
    # Export VMRs to BED file
    if mdata.uns['vmr'] is not None:
        ms2.tl.export_vmr_bed(
            mdata.uns['vmr'],
            'results/vmrs.bed'
        )
    
    # ====================================================================
    # STEP 4: Methylation Quantification
    # ====================================================================
    print("\n[4/8] Quantifying methylation across VMRs...")
    
    # Quantify using shrunken residuals method
    ms2.tl.quantify_methylation(
        mdata,
        regions='vmr',
        method='shrunken_residuals',
        shrinkage=0.5
    )
    
    print(f"Quantified: {mdata.n_obs} cells × {mdata.n_vars} VMRs")
    
    # ====================================================================
    # STEP 5: Dimensionality Reduction
    # ====================================================================
    print("\n[5/8] Performing dimensionality reduction...")
    
    # PCA
    ms2.tl.run_pca(
        mdata,
        n_comps=50,
        zero_center=True
    )
    
    # UMAP
    ms2.tl.run_umap(
        mdata,
        min_dist=0.5,
        n_neighbors=15,
        use_rep='X_pca'
    )
    
    # Visualize PCA
    ms2.pl.pca(
        mdata,
        color='mean_methylation',
        save='results/pca.pdf'
    )
    
    # ====================================================================
    # STEP 6: Clustering
    # ====================================================================
    print("\n[6/8] Clustering cells...")
    
    # Leiden clustering
    ms2.tl.run_leiden(
        mdata,
        resolution=0.8,
        use_rep='X_pca'
    )
    
    # Also try Louvain for comparison
    ms2.tl.run_louvain(
        mdata,
        resolution=0.8,
        use_rep='X_pca'
    )
    
    # ====================================================================
    # STEP 7: Visualization
    # ====================================================================
    print("\n[7/8] Creating visualizations...")
    
    # UMAP colored by clusters
    ms2.pl.umap(
        mdata,
        color='leiden',
        save='results/umap_leiden.pdf'
    )
    
    # UMAP colored by QC metrics
    ms2.pl.umap(
        mdata,
        color=['n_sites', 'mean_coverage', 'mean_methylation'],
        ncols=3,
        save='results/umap_qc.pdf'
    )
    
    # Heatmap of top variable VMRs
    ms2.pl.heatmap(
        mdata,
        groupby='leiden',
        n_vars=50,
        standard_scale='var',
        save='results/heatmap_vmr.pdf'
    )
    
    # QC scatter plot
    ms2.pl.qc_scatter(
        mdata,
        x='n_sites',
        y='mean_methylation',
        color='leiden',
        save='results/qc_scatter.pdf'
    )
    
    # ====================================================================
    # STEP 8: Save Results
    # ====================================================================
    print("\n[8/8] Saving results...")
    
    # Save complete analysis
    mdata.write('results/analyzed_data.h5ad')
    
    # Export cluster assignments
    mdata.obs[['leiden', 'louvain']].to_csv('results/clusters.csv')
    
    # Print summary
    print("\n" + "=" * 60)
    print("Analysis Complete!")
    print("=" * 60)
    print(f"Final dataset: {mdata.n_obs} cells × {mdata.n_vars} VMRs")
    print(f"Number of clusters (Leiden): {mdata.obs['leiden'].nunique()}")
    print(f"\nResults saved to: results/")
    print("  - analyzed_data.h5ad    : Complete analysis object")
    print("  - vmrs.bed              : Detected VMRs")
    print("  - clusters.csv          : Cluster assignments")
    print("  - *.pdf                 : Visualizations")
    print("=" * 60)


if __name__ == '__main__':
    main()
