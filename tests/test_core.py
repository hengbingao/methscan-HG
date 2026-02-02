"""
Basic tests for MethSCAn2
"""

import pytest
import numpy as np
import pandas as pd
from methscan2.core.methylation_data import MethylationData


def test_methylation_data_creation():
    """Test MethylationData object creation"""
    # Create test data
    n_cells = 10
    n_sites = 100
    
    X = np.random.rand(n_cells, n_sites)
    obs = pd.DataFrame({'cell_id': [f'cell_{i}' for i in range(n_cells)]})
    var = pd.DataFrame({
        'chrom': ['chr1'] * n_sites,
        'start': np.arange(n_sites) * 100,
        'end': (np.arange(n_sites) + 1) * 100
    })
    
    mdata = MethylationData(X=X, obs=obs, var=var)
    
    assert mdata.n_obs == n_cells
    assert mdata.n_vars == n_sites
    assert 'chrom' in mdata.var.columns


def test_methylation_rate():
    """Test methylation rate calculation"""
    n_cells = 5
    n_sites = 10
    
    met = np.array([[1, 0, 2], [0, 1, 1], [2, 2, 0], [1, 1, 1], [0, 0, 2]])
    total = np.array([[2, 1, 4], [1, 2, 2], [4, 4, 1], [2, 2, 2], [1, 1, 4]])
    
    X = np.zeros((n_cells, 3))
    var = pd.DataFrame({
        'chrom': ['chr1'] * 3,
        'start': [0, 100, 200],
        'end': [1, 101, 201]
    })
    
    mdata = MethylationData(
        X=X,
        var=var,
        layers={'met': met, 'total': total}
    )
    
    rate = mdata.methylation_rate
    
    # Check dimensions
    assert rate.shape == (n_cells, 3)
    
    # Check values
    expected = met.astype(float) / total
    np.testing.assert_array_almost_equal(rate, expected)


def test_subset_by_region():
    """Test region subsetting"""
    n_cells = 5
    n_sites = 100
    
    X = np.random.rand(n_cells, n_sites)
    var = pd.DataFrame({
        'chrom': ['chr1'] * 50 + ['chr2'] * 50,
        'start': list(range(0, 5000, 100)) + list(range(0, 5000, 100)),
        'end': list(range(100, 5100, 100)) + list(range(100, 5100, 100))
    })
    
    mdata = MethylationData(X=X, var=var)
    
    # Subset chr1
    subset = mdata.subset_by_region(chrom='chr1')
    assert subset.n_vars == 50
    assert (subset.var['chrom'] == 'chr1').all()
    
    # Subset region
    subset2 = mdata.subset_by_region(chrom='chr1', start=1000, end=2000)
    assert subset2.n_vars < 50


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
