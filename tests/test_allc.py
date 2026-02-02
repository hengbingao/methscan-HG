"""
Tests for ALLC processing functions
"""

import pytest
import numpy as np
import pandas as pd
import tempfile
import gzip
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from methscan2.preprocessing.allc_convert import (
    allc_to_cov,
    calculate_global_methylation,
    _process_single_allc,
    _calculate_single_global
)


def create_test_allc(filepath, n_sites=100):
    """Create a test ALLC file"""
    data = []
    for i in range(n_sites):
        # Simulate CG and CH sites
        if i % 2 == 0:
            context = 'CGN'
        else:
            context = 'CHN'
        
        mc = np.random.randint(0, 10)
        cov = np.random.randint(mc, 20)
        meth = mc / cov if cov > 0 else 0
        
        data.append({
            'chr': 'chr1',
            'pos': i * 100,
            'strand': '+',
            'context': context,
            'mc': mc,
            'cov': cov,
            'meth': meth
        })
    
    df = pd.DataFrame(data)
    
    # Write as gzipped file
    with gzip.open(filepath, 'wt') as f:
        df.to_csv(f, sep='\t', header=False, index=False)


def create_test_cov(filepath, n_sites=100):
    """Create a test COV file"""
    data = []
    for i in range(n_sites):
        met = np.random.randint(0, 10)
        unmet = np.random.randint(0, 10)
        meth_pct = 100 if met > 0 else 0
        
        data.append({
            'chr': 'chr1',
            'start': i * 100,
            'end': i * 100 + 1,
            'meth_pct': meth_pct,
            'met': met,
            'unmet': unmet
        })
    
    df = pd.DataFrame(data)
    df.to_csv(filepath, sep='\t', header=False, index=False)


def test_process_single_allc():
    """Test single ALLC file processing"""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create test ALLC file
        allc_file = Path(tmpdir) / 'test.allc.tsv.gz'
        create_test_allc(allc_file, n_sites=50)
        
        # Process it
        result = _process_single_allc(
            filepath=str(allc_file),
            output_dir=tmpdir,
            contexts=['CG', 'CH'],
            compress=False,
            verbose=False
        )
        
        # Check outputs
        assert 'CG' in result
        assert 'CH' in result
        
        # Check files exist
        assert Path(result['CG']).exists()
        assert Path(result['CH']).exists()
        
        # Check file contents
        cg_df = pd.read_csv(
            result['CG'],
            sep='\t',
            header=None,
            names=['chr', 'start', 'end', 'meth_pct', 'met', 'unmet']
        )
        
        assert len(cg_df) > 0
        assert all(cg_df['chr'] == 'chr1')
        assert all(cg_df['meth_pct'].isin([0, 100]))


def test_allc_to_cov():
    """Test batch ALLC to COV conversion"""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create multiple test ALLC files
        n_files = 5
        for i in range(n_files):
            allc_file = Path(tmpdir) / f'cell_{i}.allc.tsv.gz'
            create_test_allc(allc_file, n_sites=30)
        
        # Convert all files
        output_files = allc_to_cov(
            allc_files=str(Path(tmpdir) / '*.allc.tsv.gz'),
            output_dir=tmpdir,
            contexts=['CG', 'CH'],
            n_jobs=2,
            verbose=False
        )
        
        # Check outputs
        assert 'CG' in output_files
        assert 'CH' in output_files
        assert len(output_files['CG']) == n_files
        assert len(output_files['CH']) == n_files


def test_calculate_single_global():
    """Test global methylation calculation for single file"""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create test COV file
        cov_file = Path(tmpdir) / 'test.cov'
        create_test_cov(cov_file, n_sites=100)
        
        # Calculate global methylation
        cell_name, global_level, n_sites = _calculate_single_global(str(cov_file))
        
        assert cell_name == 'test'
        assert 0 <= global_level <= 1
        assert n_sites == 100


def test_calculate_global_methylation():
    """Test batch global methylation calculation"""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create multiple test COV files
        n_files = 5
        for i in range(n_files):
            cov_file = Path(tmpdir) / f'cell_{i}.cov'
            create_test_cov(cov_file, n_sites=50)
        
        # Calculate global methylation
        output_file = Path(tmpdir) / 'global_stats.tsv'
        stats = calculate_global_methylation(
            cov_files=str(Path(tmpdir) / '*.cov'),
            output_file=str(output_file),
            context_name='CG',
            n_jobs=2,
            verbose=False
        )
        
        # Check results
        assert len(stats) == n_files
        assert 'cell' in stats.columns
        assert 'global_CG_level' in stats.columns
        assert 'CG_sites' in stats.columns
        
        # Check file was created
        assert output_file.exists()
        
        # Check values are reasonable
        assert all(stats['global_CG_level'] >= 0)
        assert all(stats['global_CG_level'] <= 1)
        assert all(stats['CG_sites'] > 0)


def test_context_separation():
    """Test that CG and CH contexts are properly separated"""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create ALLC with mixed contexts
        allc_file = Path(tmpdir) / 'test.allc.tsv.gz'
        
        data = []
        # Add CG sites
        for i in range(10):
            data.append(['chr1', i*100, '+', 'CGN', 5, 10, 0.5])
        # Add CH sites
        for i in range(10, 20):
            data.append(['chr1', i*100, '+', 'CAT', 2, 10, 0.2])
        
        df = pd.DataFrame(data)
        with gzip.open(allc_file, 'wt') as f:
            df.to_csv(f, sep='\t', header=False, index=False)
        
        # Process
        result = _process_single_allc(
            filepath=str(allc_file),
            output_dir=tmpdir,
            contexts=['CG', 'CH'],
            compress=False,
            verbose=False
        )
        
        # Check CG file has only CG sites
        cg_df = pd.read_csv(result['CG'], sep='\t', header=None)
        assert len(cg_df) == 10
        
        # Check CH file has only CH sites
        ch_df = pd.read_csv(result['CH'], sep='\t', header=None)
        assert len(ch_df) == 10


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
