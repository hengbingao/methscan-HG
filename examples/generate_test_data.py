#!/usr/bin/env python3
"""
Generate synthetic test data for MethSCAn2

This creates small synthetic .cov files for testing the pipeline.
"""

import numpy as np
import pandas as pd
from pathlib import Path


def generate_test_data(
    n_cells=20,
    n_sites_per_chr=1000,
    chromosomes=['chr1', 'chr2', 'chr3'],
    output_dir='data',
    coverage_mean=10,
    coverage_std=3
):
    """
    Generate synthetic methylation data
    
    Parameters:
        n_cells: Number of cells to generate
        n_sites_per_chr: Number of CpG sites per chromosome
        chromosomes: List of chromosome names
        output_dir: Output directory
        coverage_mean: Mean coverage depth
        coverage_std: Standard deviation of coverage
    """
    
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True, parents=True)
    
    print(f"Generating test data for {n_cells} cells...")
    print(f"Chromosomes: {chromosomes}")
    print(f"Sites per chromosome: {n_sites_per_chr}")
    
    # Generate cell types for realistic clustering
    n_types = 3
    cell_types = np.random.choice(n_types, n_cells)
    
    for cell_idx in range(n_cells):
        cell_type = cell_types[cell_idx]
        
        # Each cell type has different methylation patterns
        base_meth_rate = 0.3 + (cell_type * 0.15)  # 0.3, 0.45, 0.6
        
        rows = []
        
        for chrom in chromosomes:
            # Generate CpG positions
            positions = np.sort(np.random.choice(
                100000000,  # Chromosome length
                n_sites_per_chr,
                replace=False
            ))
            
            for pos in positions:
                # Coverage
                coverage = max(1, int(np.random.normal(coverage_mean, coverage_std)))
                
                # Methylation rate varies by cell type and adds noise
                meth_rate = np.clip(
                    base_meth_rate + np.random.normal(0, 0.1),
                    0, 1
                )
                
                # Methylated and unmethylated counts
                met_count = int(np.random.binomial(coverage, meth_rate))
                unmet_count = coverage - met_count
                
                # Calculate percentage
                meth_pct = (met_count / coverage * 100) if coverage > 0 else 0
                
                rows.append({
                    'chrom': chrom,
                    'start': pos,
                    'end': pos + 1,
                    'meth_pct': meth_pct,
                    'met_count': met_count,
                    'unmet_count': unmet_count
                })
        
        # Create DataFrame and save
        df = pd.DataFrame(rows)
        df = df.sort_values(['chrom', 'start'])
        
        filename = output_path / f'cell_{cell_idx:03d}.cov'
        df.to_csv(
            filename,
            sep='\t',
            header=False,
            index=False
        )
        
        if (cell_idx + 1) % 5 == 0:
            print(f"  Generated {cell_idx + 1}/{n_cells} cells...")
    
    print(f"\n✓ Generated {n_cells} .cov files in {output_dir}/")
    print(f"  Total sites per cell: ~{len(chromosomes) * n_sites_per_chr}")
    print(f"  Expected cell types: {n_types}")
    print(f"\nYou can now run:")
    print(f"  python examples/complete_analysis.py")


if __name__ == '__main__':
    generate_test_data(
        n_cells=20,
        n_sites_per_chr=500,  # Small for quick testing
        output_dir='data'
    )
