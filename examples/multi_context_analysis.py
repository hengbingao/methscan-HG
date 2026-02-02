#!/usr/bin/env python3
"""
Complete Multi-Context Methylation Analysis Pipeline

This script demonstrates:
1. Loading ALLC files and separating contexts (CG, CH, all)
2. Calculating methylation across genomic features (TSS, gene body, etc.)
3. Storing all contexts in separate layers of one AnnData object
4. Running analysis on specific contexts
5. Comparing contexts
6. Visualization

Usage:
    python multi_context_analysis.py --allc raw/allc --gtf genes.gtf --output results
"""

import methscan2 as ms2
import numpy as np
import pandas as pd
from pathlib import Path
import argparse


def main(args):
    """Run complete multi-context analysis"""
    
    print("=" * 70)
    print("Multi-Context Methylation Analysis Pipeline")
    print("=" * 70)
    
    # ====================================================================
    # STEP 1: Load Genomic Annotations
    # ====================================================================
    print("\n[1/8] Loading genomic annotations...")
    
    from methscan2.preprocessing.genomic_annotation import GenomicRegionAnnotator
    
    annotator = GenomicRegionAnnotator(
        gtf_file=args.gtf,
        genome=args.genome
    )
    
    # Define genomic regions
    tss_regions = annotator.define_tss_regions(
        upstream=2000,
        downstream=500
    )
    
    gene_body_regions = annotator.define_gene_body_regions()
    
    up_down_regions = annotator.define_upstream_downstream_regions(
        upstream=5000,
        downstream=5000
    )
    
    regions_dict = {
        'tss': tss_regions,
        'gene_body': gene_body_regions,
        'upstream': up_down_regions['upstream'],
        'downstream': up_down_regions['downstream']
    }
    
    print(f"Defined {len(regions_dict)} region types:")
    for name, regions in regions_dict.items():
        print(f"  {name}: {len(regions)} regions")
    
    # ====================================================================
    # STEP 2: Create Multi-Context Data from ALLC Files
    # ====================================================================
    print("\n[2/8] Creating multi-context methylation data...")
    
    from methscan2.core.multi_context_data import create_multi_context_data_from_allc
    
    # Process for each region type
    mdata_dict = {}
    
    for region_name in args.region_types:
        print(f"\n  Processing {region_name} regions...")
        
        # Convert PyRanges to DataFrame
        regions_pr = regions_dict[region_name]
        regions_df = regions_pr.df[['Chromosome', 'Start', 'End']].copy()
        regions_df.columns = ['chr', 'start', 'end']
        
        # Create multi-context data
        mdata = create_multi_context_data_from_allc(
            allc_files=f'{args.allc}/*.allc.tsv.gz',
            regions=regions_df,
            contexts=args.contexts,
            region_type=region_name,
            min_coverage=args.min_coverage,
            n_jobs=args.jobs,
            verbose=True
        )
        
        mdata_dict[region_name] = mdata
    
    # Use TSS regions as primary for analysis
    mdata = mdata_dict[args.primary_region]
    
    print(f"\nPrimary analysis dataset: {args.primary_region}")
    print(mdata)
    
    # ====================================================================
    # STEP 3: Calculate QC Metrics for Each Context
    # ====================================================================
    print("\n[3/8] Calculating QC metrics for each context...")
    
    for ctx in mdata.available_contexts:
        print(f"\n  Context: {ctx}")
        
        # Get context data
        rate = mdata.get_context(ctx, 'rate')
        total = mdata.get_context(ctx, 'total')
        
        # Calculate metrics
        n_sites = (total > 0).sum(axis=1)
        mean_cov = np.nanmean(np.where(total > 0, total, np.nan), axis=1)
        mean_meth = np.nanmean(rate, axis=1)
        
        # Add to obs
        mdata.obs[f'{ctx}_n_sites'] = n_sites
        mdata.obs[f'{ctx}_mean_coverage'] = mean_cov
        mdata.obs[f'{ctx}_mean_methylation'] = mean_meth
        
        print(f"    Mean sites: {n_sites.mean():.0f}")
        print(f"    Mean methylation: {mean_meth.mean():.3f}")
    
    # ====================================================================
    # STEP 4: Filter Cells Based on QC
    # ====================================================================
    print("\n[4/8] Filtering cells...")
    
    print(f"Before filtering: {mdata.n_obs} cells")
    
    # Filter based on primary context (usually CG)
    primary_ctx = args.contexts[0]
    
    keep_mask = (
        (mdata.obs[f'{primary_ctx}_n_sites'] >= args.min_sites) &
        (mdata.obs[f'{primary_ctx}_mean_methylation'] >= 0.3) &
        (mdata.obs[f'{primary_ctx}_mean_methylation'] <= 0.8)
    )
    
    mdata = mdata[keep_mask, :].copy()
    
    print(f"After filtering: {mdata.n_obs} cells")
    
    # ====================================================================
    # STEP 5: Run Analysis on Each Context
    # ====================================================================
    print("\n[5/8] Running dimensionality reduction and clustering...")
    
    for ctx in args.analysis_contexts:
        print(f"\n  Analyzing context: {ctx}")
        
        # Set as active context
        mdata.set_active_context(ctx)
        
        # PCA
        ms2.tl.run_pca(
            mdata,
            n_comps=50,
            zero_center=True
        )
        
        # Store PCA results with context name
        mdata.obsm[f'X_pca_{ctx}'] = mdata.obsm['X_pca'].copy()
        
        # UMAP
        ms2.tl.run_umap(
            mdata,
            use_rep='X_pca',
            min_dist=0.5,
            n_neighbors=15
        )
        
        # Store UMAP results
        mdata.obsm[f'X_umap_{ctx}'] = mdata.obsm['X_umap'].copy()
        
        # Clustering
        ms2.tl.run_leiden(
            mdata,
            resolution=args.resolution,
            use_rep='X_pca',
            key_added=f'leiden_{ctx}'
        )
        
        print(f"    Found {mdata.obs[f'leiden_{ctx}'].nunique()} clusters")
    
    # ====================================================================
    # STEP 6: Compare Contexts
    # ====================================================================
    print("\n[6/8] Comparing methylation contexts...")
    
    # Calculate correlations between contexts
    corr_df = mdata.calculate_context_correlations()
    corr_df.to_csv(f'{args.output}/context_correlations.csv', index=False)
    
    print("\nContext correlations:")
    print(corr_df)
    
    # ====================================================================
    # STEP 7: Visualization
    # ====================================================================
    print("\n[7/8] Creating visualizations...")
    
    Path(args.output).mkdir(parents=True, exist_ok=True)
    
    # QC plots
    print("  Generating QC plots...")
    ms2.pl.qc_violin(
        mdata,
        keys=[f'{ctx}_mean_methylation' for ctx in mdata.available_contexts],
        save=f'{args.output}/qc_methylation_by_context.pdf'
    )
    
    # Context comparison
    print("  Comparing contexts...")
    mdata.plot_context_comparison(
        metric='mean',
        save=f'{args.output}/context_comparison.pdf'
    )
    
    # UMAP for each context
    for ctx in args.analysis_contexts:
        print(f"  Plotting UMAP for {ctx}...")
        
        # Set active context for visualization
        mdata.set_active_context(ctx)
        
        # UMAP colored by clusters
        ms2.pl.umap(
            mdata,
            color=f'leiden_{ctx}',
            use_rep=f'X_umap_{ctx}',
            save=f'{args.output}/umap_{ctx}_clusters.pdf'
        )
        
        # UMAP colored by methylation level
        ms2.pl.umap(
            mdata,
            color=f'{ctx}_mean_methylation',
            use_rep=f'X_umap_{ctx}',
            save=f'{args.output}/umap_{ctx}_methylation.pdf'
        )
    
    # Multi-panel UMAP comparison
    print("  Creating multi-context UMAP comparison...")
    import matplotlib.pyplot as plt
    
    n_contexts = len(args.analysis_contexts)
    fig, axes = plt.subplots(1, n_contexts, figsize=(6*n_contexts, 5))
    
    if n_contexts == 1:
        axes = [axes]
    
    for ax, ctx in zip(axes, args.analysis_contexts):
        coords = mdata.obsm[f'X_umap_{ctx}']
        clusters = mdata.obs[f'leiden_{ctx}']
        
        scatter = ax.scatter(
            coords[:, 0], coords[:, 1],
            c=clusters.cat.codes,
            cmap='tab20',
            s=10, alpha=0.8
        )
        
        ax.set_title(f'{ctx} Context')
        ax.set_xlabel('UMAP1')
        ax.set_ylabel('UMAP2')
        ax.set_aspect('equal')
    
    plt.tight_layout()
    plt.savefig(f'{args.output}/umap_all_contexts.pdf', dpi=300)
    plt.close()
    
    # Heatmap for primary context
    print("  Creating heatmap...")
    mdata.set_active_context(args.contexts[0])
    
    ms2.pl.heatmap(
        mdata,
        groupby=f'leiden_{args.contexts[0]}',
        n_vars=50,
        save=f'{args.output}/heatmap_{args.contexts[0]}.pdf'
    )
    
    # ====================================================================
    # STEP 8: Save Results
    # ====================================================================
    print("\n[8/8] Saving results...")
    
    # Save main object
    mdata.write(f'{args.output}/multi_context_data.h5ad')
    
    # Save all region-specific objects
    for region_name, region_mdata in mdata_dict.items():
        region_mdata.write(f'{args.output}/mdata_{region_name}.h5ad')
    
    # Export cluster assignments
    cluster_df = mdata.obs[[f'leiden_{ctx}' for ctx in args.analysis_contexts]]
    cluster_df.to_csv(f'{args.output}/cluster_assignments.csv')
    
    # Export methylation statistics
    stats_df = mdata.obs[[
        f'{ctx}_mean_methylation' for ctx in mdata.available_contexts
    ]]
    stats_df.to_csv(f'{args.output}/methylation_statistics.csv')
    
    # ====================================================================
    # Summary
    # ====================================================================
    print("\n" + "=" * 70)
    print("Analysis Complete!")
    print("=" * 70)
    
    print(f"\nDataset summary:")
    print(f"  Cells: {mdata.n_obs}")
    print(f"  Regions: {mdata.n_vars} ({args.primary_region})")
    print(f"  Contexts: {', '.join(mdata.available_contexts)}")
    
    print(f"\nClusters found:")
    for ctx in args.analysis_contexts:
        n_clusters = mdata.obs[f'leiden_{ctx}'].nunique()
        print(f"  {ctx}: {n_clusters} clusters")
    
    print(f"\nOutput files:")
    print(f"  {args.output}/")
    print(f"    ├── multi_context_data.h5ad         (Main analysis object)")
    for region in args.region_types:
        print(f"    ├── mdata_{region}.h5ad             ({region} regions)")
    print(f"    ├── cluster_assignments.csv          (Cluster labels)")
    print(f"    ├── methylation_statistics.csv       (Methylation stats)")
    print(f"    ├── context_correlations.csv         (Context correlations)")
    print(f"    └── *.pdf                            (Visualizations)")
    
    print("\n" + "=" * 70)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Multi-context methylation analysis'
    )
    
    # Input
    parser.add_argument(
        '--allc', '-i',
        required=True,
        help='Directory with ALLC files'
    )
    
    parser.add_argument(
        '--gtf', '-g',
        required=True,
        help='GTF annotation file'
    )
    
    parser.add_argument(
        '--genome',
        default='hg38',
        help='Genome assembly (default: hg38)'
    )
    
    # Output
    parser.add_argument(
        '--output', '-o',
        default='results',
        help='Output directory (default: results)'
    )
    
    # Contexts
    parser.add_argument(
        '--contexts',
        nargs='+',
        default=['CG', 'CH', 'all'],
        help='Methylation contexts (default: CG CH all)'
    )
    
    parser.add_argument(
        '--analysis-contexts',
        nargs='+',
        default=['CG', 'CH'],
        help='Contexts to analyze (default: CG CH)'
    )
    
    # Regions
    parser.add_argument(
        '--region-types',
        nargs='+',
        default=['tss', 'gene_body', 'upstream', 'downstream'],
        help='Region types to process'
    )
    
    parser.add_argument(
        '--primary-region',
        default='tss',
        help='Primary region for analysis (default: tss)'
    )
    
    # Filtering
    parser.add_argument(
        '--min-coverage',
        type=int,
        default=1,
        help='Minimum coverage (default: 1)'
    )
    
    parser.add_argument(
        '--min-sites',
        type=int,
        default=1000,
        help='Minimum sites per cell (default: 1000)'
    )
    
    # Analysis
    parser.add_argument(
        '--resolution',
        type=float,
        default=0.8,
        help='Clustering resolution (default: 0.8)'
    )
    
    parser.add_argument(
        '--jobs', '-j',
        type=int,
        default=-1,
        help='Parallel jobs (default: all CPUs)'
    )
    
    args = parser.parse_args()
    
    main(args)
