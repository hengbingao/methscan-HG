# Installation and Testing Guide

## System Requirements

### Hardware
- RAM: Minimum 16GB, recommended 32GB+ for large datasets
- CPU: Multi-core processor recommended (8+ cores optimal)
- Storage: Depends on dataset size, recommend 100GB+ free space

### Software
- Python 3.8 or higher
- bedtools (system dependency)
- Optional: allcools (for allc file processing)

## Installation Steps

### 1. Install System Dependencies

#### Ubuntu/Debian
```bash
sudo apt-get update
sudo apt-get install -y bedtools samtools tabix
```

#### macOS (using Homebrew)
```bash
brew install bedtools samtools htslib
```

#### Conda (cross-platform)
```bash
conda install -c bioconda bedtools samtools tabix
```

### 2. Install allcools (optional)

If you need to process raw allc files:

```bash
# Via conda
conda install -c bioconda allcools

# Or via pip
pip install allcools
```

### 3. Install scMethylSeq

#### From source (recommended for development)
```bash
git clone https://github.com/yourusername/scMethylSeq.git
cd scMethylSeq
pip install -e .
```

#### From PyPI (when available)
```bash
pip install scmethylseq
```

### 4. Verify Installation

```bash
# Check if command-line tool is available
scmethyl --help

# Check Python import
python -c "import scmethylseq; print(scmethylseq.__version__)"

# Check dependencies
python -c "from scmethylseq.utils import print_dependency_status; print_dependency_status()"
```

## Quick Test

### Test with Example Data

Create a test script:

```python
# test_installation.py
import numpy as np
import pandas as pd
from scmethylseq import PCAAnalyzer, UMAPAnalyzer

# Create dummy methylation data
n_cells, n_genes = 100, 1000
data = np.random.rand(n_cells, n_genes)
data[np.random.rand(n_cells, n_genes) < 0.1] = np.nan  # Add missing values

# Test PCA
print("Testing PCA...")
pca = PCAAnalyzer(n_components=10)
pca_coords, pca_df = pca.fit_transform(pd.DataFrame(data))
print(f"PCA output shape: {pca_coords.shape}")
print("✓ PCA test passed")

# Test UMAP
print("\nTesting UMAP...")
umap = UMAPAnalyzer(n_neighbors=15, random_state=42)
umap_coords, graph = umap.fit_transform(pca_coords)
print(f"UMAP output shape: {umap_coords.shape}")
print("✓ UMAP test passed")

print("\n✓ All tests passed!")
```

Run the test:
```bash
python test_installation.py
```

## Troubleshooting

### Common Issues

#### 1. "bedtools not found"
**Solution**: Install bedtools using your system package manager (see Step 1)

#### 2. "ImportError: No module named 'umap'"
**Solution**: Install umap-learn
```bash
pip install umap-learn
```

#### 3. "ImportError: No module named 'leidenalg'"
**Solution**: Install leidenalg and igraph
```bash
pip install leidenalg python-igraph
```

#### 4. Memory Error during processing
**Solutions**:
- Reduce `n_jobs` parameter
- Process smaller batches of cells
- Use a machine with more RAM
- Enable `cleanup_intermediate: true` in config

#### 5. "Permission denied" errors
**Solution**: Check file permissions and ensure output directories are writable
```bash
chmod -R 755 /path/to/output
```

## Testing the Complete Pipeline

### 1. Prepare Test Data

You'll need:
- A few allc.gz files (or start from cov files)
- Chromosome sizes file
- Gene annotation BED file

### 2. Create Test Configuration

```yaml
# test_config.yaml
run_allc_processing: false  # Skip if starting from cov
run_cov_conversion: false
run_global_meth: true
run_qc: true
run_gene_meth: true
run_matrix_building: true
run_adata_building: true

n_jobs: 4  # Use fewer jobs for testing

cg_cov_dir: "/path/to/test/CG_cov"
ch_cov_dir: "/path/to/test/CH_cov"
gene_bed: "/path/to/test/genes.bed"
gene_names_file: "/path/to/test/gene_names.txt"
output_base: "/path/to/test/output"

# ... rest of config
```

### 3. Run Test Pipeline

```bash
scmethyl run-pipeline --config test_config.yaml
```

### 4. Verify Results

Check that output files are created:
```bash
# Check QC results
ls /path/to/test/output/qc_results/

# Check matrices
ls /path/to/test/output/matrices/

# Check final AnnData
ls /path/to/test/output/*.h5ad
```

Load and inspect the results:
```python
import scanpy as sc

adata = sc.read_h5ad('/path/to/test/output/final.h5ad')
print(adata)
print(adata.obs.columns)
print(adata.obsm.keys())
```

## Performance Benchmarks

Typical processing times on a 32-core, 64GB RAM machine:

| Step | 1,000 cells | 10,000 cells | 50,000 cells |
|------|-------------|--------------|--------------|
| allc processing | 10 min | 1.5 hours | 8 hours |
| cov conversion | 2 min | 20 min | 2 hours |
| Global methylation | 1 min | 10 min | 1 hour |
| Gene body calculation | 15 min | 2.5 hours | 12 hours |
| Matrix building | 1 min | 5 min | 30 min |
| PCA + UMAP + clustering | 30 sec | 5 min | 30 min |
| **Total** | **~30 min** | **~5 hours** | **~24 hours** |

Note: Times vary based on:
- Number of genomic features
- Coverage depth
- System specifications
- Number of parallel jobs

## Optimization Tips

### 1. Maximize Parallel Processing
```yaml
n_jobs: 30  # Use most of your CPU cores
```

### 2. Use SSD Storage
Process data on SSD drives for faster I/O

### 3. Filter Cells Early
Run QC first and filter cells before gene body calculation

### 4. Process Contexts Separately
Run CG and CH pipelines independently to reduce peak memory

### 5. Use Checkpoints
Enable progress tracking to resume from failures:
```python
from scmethylseq.utils import ProgressTracker

tracker = ProgressTracker('pipeline_progress.txt')
# Check completion before running steps
if not tracker.is_complete('step1'):
    # run step1
    tracker.mark_complete('step1')
```

## Getting Help

### Documentation
- README.md: Overview and quick start
- This file: Detailed installation guide
- examples/: Example scripts

### Issues
Report bugs or request features on GitHub:
https://github.com/yourusername/scMethylSeq/issues

### Support
Contact: your.email@example.com

## Next Steps

After successful installation:
1. Read the [User Guide](USER_GUIDE.md) for detailed usage
2. Review [examples/](examples/) for common workflows
3. Customize config_example.yaml for your data
4. Start with a small test dataset before processing full data
