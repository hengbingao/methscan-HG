#!/usr/bin/env python3
"""
Complete ALLC Processing Pipeline Example

This script demonstrates how to:
1. Convert ALLC files to COV format
2. Extract CG, CH contexts separately
3. Calculate global methylation levels
4. Generate quality control reports

Replicates the bash workflow in pure Python with parallelization.
"""

import methscan2 as ms2
import os
from pathlib import Path


def main():
    """Run complete ALLC processing pipeline"""
    
    print("=" * 70)
    print("MethSCAn2 - ALLC Processing Pipeline")
    print("=" * 70)
    
    # ====================================================================
    # Configuration
    # ====================================================================
    
    # Input: Directory containing .allc.tsv.gz files
    allc_dir = 'raw/allc'
    
    # Output: Base directory for all outputs
    output_dir = 'processed'
    
    # Contexts to extract
    contexts = ['CG', 'CH', 'all']  # CG, CH, and all methylation
    
    # Number of parallel jobs (-1 = use all CPUs)
    n_jobs = 100  # Match your xargs -P 100
    
    # ====================================================================
    # OPTION 1: Complete Automated Pipeline
    # ====================================================================
    
    print("\n" + "=" * 70)
    print("Running Complete Automated Pipeline")
    print("=" * 70)
    
    # This single function call does everything:
    # 1. Converts ALLC → COV
    # 2. Separates contexts (CG, CH, all)
    # 3. Organizes into directories
    # 4. Calculates global methylation levels
    
    global_stats = ms2.pp.batch_process_allc(
        allc_dir=allc_dir,
        output_base_dir=output_dir,
        contexts=contexts,
        calculate_global=True,
        n_jobs=n_jobs,
        verbose=True
    )
    
    # Results are automatically organized as:
    # processed/
    # ├── CG_cov/
    # │   └── *.cov files
    # ├── CH_cov/
    # │   └── *.cov files
    # ├── all_cov/
    # │   └── *.cov files
    # └── global_level/
    #     ├── all_cell_global_average_mCG.tsv
    #     ├── all_cell_global_average_mCH.tsv
    #     └── all_cell_global_average_mall.tsv
    
    print("\n✓ Automated pipeline complete!")
    
    # ====================================================================
    # OPTION 2: Step-by-Step Manual Control
    # ====================================================================
    
    print("\n" + "=" * 70)
    print("Alternative: Step-by-Step Processing")
    print("=" * 70)
    
    # Step 1: Convert ALLC to COV with context separation
    print("\n[Step 1] Converting ALLC to COV...")
    
    output_files = ms2.pp.allc_to_cov(
        allc_files='raw/allc/*.allc.tsv.gz',
        output_dir='manual_output',
        contexts=['CG', 'CH', 'all'],
        separate_contexts=True,
        compress=False,  # Set to True if you want .cov.gz files
        n_jobs=n_jobs,
        verbose=True
    )
    
    # Step 2: Calculate global methylation for each context
    print("\n[Step 2] Calculating global methylation levels...")
    
    for context in ['CG', 'CH', 'all']:
        print(f"\n  Processing {context} context...")
        
        stats = ms2.pp.calculate_global_methylation(
            cov_files=f'manual_output/*_{context}.cov',
            output_file=f'manual_output/global_{context}.tsv',
            context_name=context,
            n_jobs=n_jobs,
            verbose=True
        )
        
        # Display summary
        print(f"\n  {context} Statistics:")
        print(f"    Cells: {len(stats)}")
        print(f"    Mean methylation: {stats[f'global_{context}_level'].mean():.3f}")
        print(f"    Mean sites: {stats[f'{context}_sites'].mean():.0f}")
    
    # Step 3: Merge all context statistics
    print("\n[Step 3] Merging statistics...")
    
    merged_stats = ms2.pp.merge_global_stats(
        stats_files={
            'CG': 'manual_output/global_CG.tsv',
            'CH': 'manual_output/global_CH.tsv',
            'all': 'manual_output/global_all.tsv'
        },
        output_file='manual_output/merged_global_stats.tsv'
    )
    
    print(f"\n✓ Merged statistics saved!")
    print(f"  Columns: {list(merged_stats.columns)}")
    
    # ====================================================================
    # Quality Control and Filtering
    # ====================================================================
    
    print("\n" + "=" * 70)
    print("Quality Control and Filtering")
    print("=" * 70)
    
    # Load CG methylation data for QC
    print("\n[QC] Loading CG COV files...")
    
    mdata_cg = ms2.read_cov_files(
        f'{output_dir}/CG_cov/*.cov',
        genome='hg38',
        n_jobs=n_jobs,
        verbose=True
    )
    
    # Calculate QC metrics
    print("\n[QC] Calculating quality metrics...")
    ms2.pp.calculate_qc_metrics(mdata_cg)
    
    # Plot QC distributions
    print("\n[QC] Creating QC plots...")
    ms2.pl.qc_violin(
        mdata_cg,
        keys=['n_sites', 'mean_coverage', 'mean_methylation'],
        save=f'{output_dir}/qc_metrics_CG.pdf'
    )
    
    ms2.pl.qc_scatter(
        mdata_cg,
        x='n_sites',
        y='mean_methylation',
        save=f'{output_dir}/qc_scatter_CG.pdf'
    )
    
    # Filter low-quality cells
    print("\n[QC] Filtering cells...")
    print(f"Before filtering: {mdata_cg.n_obs} cells")
    
    ms2.pp.filter_cells(
        mdata_cg,
        min_sites=50000,      # At least 50k CG sites
        min_coverage=5,        # Mean coverage ≥ 5
        min_methylation=0.3,   # Typical CG methylation 30-80%
        max_methylation=0.8
    )
    
    print(f"After filtering: {mdata_cg.n_obs} cells")
    
    # Save filtered cell list
    mdata_cg.obs[['n_sites', 'mean_coverage', 'mean_methylation']].to_csv(
        f'{output_dir}/filtered_cells_CG.csv'
    )
    
    # ====================================================================
    # Advanced: Using allcools (if installed)
    # ====================================================================
    
    print("\n" + "=" * 70)
    print("Advanced: Using allcools for Context Extraction")
    print("=" * 70)
    
    try:
        # This requires allcools to be installed
        # pip install allcools
        
        print("\n[allcools] Extracting contexts with allcools...")
        
        for allc_file in Path(allc_dir).glob('*.allc.tsv.gz'):
            output_files = ms2.pp.extract_context_from_allc(
                allc_file=str(allc_file),
                contexts=['CG', 'CH'],
                output_prefix=str(allc_file.stem),
                chrom_sizes='genome/hg38.chrom.sizes',  # Optional
                cpu=15
            )
            
            print(f"  Extracted: {allc_file.name}")
            for f in output_files:
                print(f"    → {f}")
    
    except ImportError:
        print("\n  allcools not installed. Skipping this method.")
        print("  Install with: pip install allcools")
    
    # ====================================================================
    # Summary
    # ====================================================================
    
    print("\n" + "=" * 70)
    print("Pipeline Summary")
    print("=" * 70)
    
    print("\n✓ All processing complete!")
    print(f"\nOutput directory structure:")
    print(f"  {output_dir}/")
    print(f"    ├── CG_cov/              (CG methylation)")
    print(f"    ├── CH_cov/              (CH methylation)")
    print(f"    ├── all_cov/             (All methylation)")
    print(f"    ├── global_level/        (Global methylation stats)")
    print(f"    ├── qc_metrics_CG.pdf    (QC plots)")
    print(f"    ├── qc_scatter_CG.pdf")
    print(f"    └── filtered_cells_CG.csv")
    
    print("\nNext steps:")
    print("  1. Review global methylation levels in global_level/")
    print("  2. Check QC plots to identify problematic cells")
    print("  3. Use filtered COV files for downstream analysis")
    print("  4. Run: python examples/complete_analysis.py")
    
    print("\n" + "=" * 70)


if __name__ == '__main__':
    main()
