# Quick Start Guide

This guide will help you get started with scMethylSeq in 15 minutes.

## Prerequisites

- Python 3.8+
- bedtools installed
- At least 16GB RAM
- Methylation data (allc or cov files)

## Installation (2 minutes)

```bash
# Clone repository
git clone https://github.com/yourusername/scMethylSeq.git
cd scMethylSeq

# Install
pip install -e .

# Verify
scmethyl --help
```

## Scenario 1: You Have Methylation Matrices (5 minutes)

If you already have CG and CH methylation matrices:

### Step 1: Prepare metadata

Create a CSV file with cell information:
```csv
cell,leiden_cluster,UMAP1,UMAP2,batch
cell_001,0,1.23,-2.45,batch1
cell_002,1,-0.34,3.21,batch1
cell_003,0,1.45,-1.89,batch2
```

### Step 2: Build AnnData

```bash
scmethyl build-adata \
    --cg-matrix path/to/cell_mCG_matrix \
    --ch-matrix path/to/cell_mCH_matrix \
    --metadata path/to/metadata.csv \
    --output output/adata.h5ad
```

### Step 3: Load and analyze

```python
import scanpy as sc

adata = sc.read_h5ad('output/adata.h5ad')
print(adata)

# Plot UMAP
sc.pl.umap(adata, color='leiden_cluster')
```

## Scenario 2: You Have Cov Files (30 minutes)

If you have cov files for each cell:

### Step 1: Calculate global methylation

```bash
# For CG
scmethyl calc-global-meth \
    --cov-dir path/to/CG_cov \
    --output-file output/global_CG.txt \
    --context CG \
    --n-jobs 20

# For CH
scmethyl calc-global-meth \
    --cov-dir path/to/CH_cov \
    --output-file output/global_CH.txt \
    --context CH \
    --n-jobs 20
```

### Step 2: Filter cells

```bash
scmethyl filter-cells \
    --cg-file output/global_CG.txt \
    --ch-file output/global_CH.txt \
    --output-dir output/qc \
    --cg-min 0.5 \
    --ch-max 0.1
```

This creates `cell_filtered.csv` with passing cells.

### Step 3: Calculate gene body methylation

You need a gene BED file (e.g., genebody_extend2kb.bed):
```
chr1	1000	5000	GENE1
chr1	10000	15000	GENE2
```

```bash
# For CG
scmethyl calc-gene-meth \
    --cov-dir path/to/CG_cov \
    --gene-bed path/to/genes.bed \
    --output-dir output/CG_gene_meth \
    --context CG \
    --n-jobs 20

# For CH
scmethyl calc-gene-meth \
    --cov-dir path/to/CH_cov \
    --gene-bed path/to/genes.bed \
    --output-dir output/CH_gene_meth \
    --context CH \
    --n-jobs 20
```

### Step 4: Build matrices

You need a gene names file (one gene per line):
```
GENE1
GENE2
GENE3
```

```bash
# For CG
scmethyl build-matrix \
    --gene-meth-dir output/CG_gene_meth \
    --gene-names-file path/to/gene_names.txt \
    --output-file output/cell_mCG_matrix \
    --context CG \
    --cleanup

# For CH
scmethyl build-matrix \
    --gene-meth-dir output/CH_gene_meth \
    --gene-names-file path/to/gene_names.txt \
    --output-file output/cell_mCH_matrix \
    --context CH \
    --cleanup
```

### Step 5: Continue with Scenario 1

Now you have matrices and can follow Scenario 1.

## Scenario 3: Complete Pipeline with Config (Easiest!)

### Step 1: Create config file

```bash
cp config_example.yaml my_project.yaml
```

Edit the file with your paths:
```yaml
run_allc_processing: false
run_cov_conversion: false
run_global_meth: true
run_qc: true
run_gene_meth: true
run_matrix_building: true
run_adata_building: true

n_jobs: 20

cg_cov_dir: "/path/to/CG_cov"
ch_cov_dir: "/path/to/CH_cov"
gene_bed: "/path/to/genes.bed"
gene_names_file: "/path/to/gene_names.txt"
output_base: "/path/to/output"

# ... rest of paths
```

### Step 2: Run pipeline

```bash
scmethyl run-pipeline --config my_project.yaml
```

That's it! The pipeline will:
1. Calculate global methylation
2. Filter cells by QC
3. Calculate gene body methylation
4. Build matrices
5. Create AnnData object

### Step 3: Analyze results

```python
import scanpy as sc

adata = sc.read_h5ad('output/final_adata.h5ad')

# Basic stats
print(f"Cells: {adata.n_obs}")
print(f"Genes: {adata.n_vars}")

# Run analysis if not done
if 'X_pca' not in adata.obsm:
    sc.pp.neighbors(adata, n_pcs=20)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.2)

# Visualize
sc.pl.umap(adata, color=['leiden', 'global_CG', 'global_CH'])
```

## Python API Usage

For more control, use the Python API:

```python
from scmethylseq import (
    GlobalMethylationCalculator,
    QualityControl,
    GeneBodyMethylation,
    MethylationMatrixBuilder,
    AnnDataBuilder
)

# Calculate global methylation
calc = GlobalMethylationCalculator(
    cov_dir="path/to/CG_cov",
    output_file="global_CG.txt",
    context="CG",
    n_jobs=20
)
cg_df = calc.calculate_all()

# Filter cells
qc = QualityControl(
    cg_file="global_CG.txt",
    ch_file="global_CH.txt",
    output_dir="qc_results"
)
filtered = qc.filter_cells()

# Calculate gene body methylation
gene_meth = GeneBodyMethylation(
    cov_dir="path/to/CG_cov",
    gene_bed="genes.bed",
    output_dir="CG_gene_meth",
    context="CG",
    n_jobs=20
)
gene_meth.process_all_cells()

# Build matrix
builder = MethylationMatrixBuilder(
    gene_meth_dir="CG_gene_meth",
    gene_names_file="gene_names.txt",
    output_file="cell_mCG_matrix",
    context="CG"
)
matrix = builder.build_matrix()

# Build AnnData
adata_builder = AnnDataBuilder(
    cg_matrix_file="cell_mCG_matrix",
    ch_matrix_file="cell_mCH_matrix",
    metadata_file="metadata.csv",
    output_file="adata.h5ad"
)
adata = adata_builder.build_and_save()
```

## Common Workflows

### Add clustering to existing AnnData

```python
import scanpy as sc

adata = sc.read_h5ad('adata.h5ad')

# Run analysis
sc.pp.neighbors(adata, n_pcs=20, n_neighbors=200)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.2)

# Save
adata.write('adata_clustered.h5ad')
```

### Export cluster-specific data

```python
import scanpy as sc

adata = sc.read_h5ad('adata.h5ad')

# Get cells from cluster 0
cluster0_cells = adata[adata.obs['leiden'] == '0']

# Export
cluster0_cells.write('cluster0.h5ad')
```

### Merge multiple datasets

```python
import scanpy as sc

adata1 = sc.read_h5ad('dataset1.h5ad')
adata2 = sc.read_h5ad('dataset2.h5ad')

# Concatenate
adata_merged = adata1.concatenate(adata2, batch_key='dataset')

# Batch correction if needed
import scanpy.external as sce
sce.pp.harmony_integrate(adata_merged, 'dataset')
```

## Next Steps

1. **Read the full documentation**: README.md and INSTALL.md
2. **Explore examples**: Check examples/ directory
3. **Customize parameters**: Adjust QC thresholds, n_jobs, etc.
4. **Integrate with scanpy**: Use full scanpy ecosystem

## Tips for Success

1. **Start small**: Test with 100-1000 cells first
2. **Monitor memory**: Use `top` or `htop` to watch RAM usage
3. **Use progress tracking**: Check log files for errors
4. **Adjust n_jobs**: More is not always better (I/O bottleneck)
5. **Save intermediate results**: Don't recompute everything

## Getting Help

- Check logs for error messages
- Review INSTALL.md for troubleshooting
- Open an issue on GitHub
- Email: your.email@example.com

Happy analyzing! 🧬
