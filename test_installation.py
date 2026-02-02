#!/usr/bin/env python3
"""
Quick test to verify MethSCAn2 installation

This script creates synthetic data and runs a minimal analysis.
"""

import numpy as np
import pandas as pd

def test_installation():
    """Test if MethSCAn2 is properly installed"""
    
    print("=" * 60)
    print("MethSCAn2 Installation Test")
    print("=" * 60)
    
    # Test 1: Import
    print("\n[1/5] Testing import...")
    try:
        import methscan2 as ms2
        print(f"✓ MethSCAn2 version {ms2.__version__} imported successfully")
    except ImportError as e:
        print(f"✗ Failed to import MethSCAn2: {e}")
        return False
    
    # Test 2: Create data
    print("\n[2/5] Creating test data...")
    try:
        n_cells = 10
        n_sites = 100
        
        # Simulate methylation data
        X = np.random.rand(n_cells, n_sites) * 0.5 + 0.25  # 25-75% methylation
        
        var = pd.DataFrame({
            'chrom': ['chr1'] * n_sites,
            'start': np.arange(n_sites) * 1000,
            'end': (np.arange(n_sites) + 1) * 1000
        })
        
        met = np.random.binomial(10, X)
        total = np.ones((n_cells, n_sites), dtype=int) * 10
        
        mdata = ms2.MethylationData(
            X=X,
            var=var,
            layers={'met': met, 'total': total}
        )
        
        print(f"✓ Created MethylationData: {mdata.n_obs} cells × {mdata.n_vars} sites")
    except Exception as e:
        print(f"✗ Failed to create data: {e}")
        return False
    
    # Test 3: QC
    print("\n[3/5] Testing QC functions...")
    try:
        ms2.pp.calculate_qc_metrics(mdata)
        print(f"✓ QC metrics calculated")
        print(f"  Mean sites per cell: {mdata.obs['n_sites'].mean():.1f}")
        print(f"  Mean methylation: {mdata.obs['mean_methylation'].mean():.2%}")
    except Exception as e:
        print(f"✗ QC failed: {e}")
        return False
    
    # Test 4: Dimensionality reduction
    print("\n[4/5] Testing dimensionality reduction...")
    try:
        ms2.tl.run_pca(mdata, n_comps=5)
        ms2.tl.run_umap(mdata, n_components=2)
        print(f"✓ PCA and UMAP completed")
        print(f"  PCA shape: {mdata.obsm['X_pca'].shape}")
        print(f"  UMAP shape: {mdata.obsm['X_umap'].shape}")
    except Exception as e:
        print(f"✗ Dimensionality reduction failed: {e}")
        return False
    
    # Test 5: Clustering
    print("\n[5/5] Testing clustering...")
    try:
        ms2.tl.run_leiden(mdata, resolution=0.5)
        n_clusters = mdata.obs['leiden'].nunique()
        print(f"✓ Clustering completed")
        print(f"  Number of clusters: {n_clusters}")
    except Exception as e:
        print(f"✗ Clustering failed: {e}")
        return False
    
    # Success
    print("\n" + "=" * 60)
    print("✓ All tests passed! MethSCAn2 is ready to use.")
    print("=" * 60)
    print("\nNext steps:")
    print("  1. Generate test data: python examples/generate_test_data.py")
    print("  2. Run full analysis: python examples/complete_analysis.py")
    print("  3. Check documentation: README.md")
    print("=" * 60)
    
    return True


if __name__ == '__main__':
    success = test_installation()
    exit(0 if success else 1)
