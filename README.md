# scMethylSeq: Single-cell DNA Methylation Analysis Toolkit

A comprehensive Python toolkit for analyzing single-cell DNA methylation sequencing data, from allc files to AnnData objects with clustering and visualization.

## Features

- **Fast parallel processing**: Multi-threaded/multi-process execution for all computationally intensive steps
- **Complete pipeline**: From raw allc files to final AnnData object
- **Quality control**: Automatic cell filtering based on methylation metrics
- **Flexible**: Run entire pipeline or individual steps
- **Methylation contexts**: Support for both CG and CH methylation
- **Integration with scanpy**: Compatible with standard single-cell analysis workflows

## Installation

### Requirements

- Python >= 3.8
- bedtools (system dependency)
- allcools (optional, for allc processing)

### Install from source

```bash
git clone https://github.com/yourusername/scMethylSeq.git
cd scMethylSeq
pip install -e .
```

### Dependencies

The package will automatically install Python dependencies:
- numpy, pandas, scipy
- scikit-learn
- scanpy, anndata
- pybedtools, pysam
- umap-learn, leidenalg, igraph
- joblib, tqdm, click

## Quick Start

### 1. Complete Pipeline with Configuration File

The easiest way to run the entire pipeline:

```bash
# Edit the configuration file
cp config_example.yaml my_config.yaml
# ... edit my_config.yaml with your paths ...

# Run the pipeline
scmethyl run-pipeline --config my_config.yaml
```

### 2. Step-by-Step Execution

You can also run individual steps:

#### Step 1: Process allc files
```bash
scmethyl process-allc \
    --allc-dir /path/to/allc/files \
    --output-dir /path/to/output/processed_allc \
    --chrom-sizes /path/to/hg38.chrom.sizes \
    --n-jobs 10
```

#### Step 2: Convert to cov format
```bash
scmethyl convert-to-cov \
    --allc-dir /path/to/processed_allc \
    --output-dir /path/to/output/CG_cov \
    --context CG \
    --n-jobs 20
```

#### Step 3: Calculate global methylation
```bash
scmethyl calc-global-meth \
    --cov-dir /path/to/CG_cov \
    --output-file /path/to/output/global_CG.txt \
    --context CG \
    --n-jobs 20
```

#### Step 4: Filter cells based on QC
```bash
scmethyl filter-cells \
    --cg-file /path/to/global_CG.txt \
    --ch-file /path/to/global_CH.txt \
    --output-dir /path/to/output/qc_results \
    --cg-min 0.5 \
    --ch-max 0.1
```

#### Step 5: Calculate gene body methylation
```bash
scmethyl calc-gene-meth \
    --cov-dir /path/to/CG_cov \
    --gene-bed /path/to/genebody_extend2kb.bed \
    --output-dir /path/to/output/CG_gene_meth \
    --context CG \
    --n-jobs 20
```

#### Step 6: Build methylation matrix
```bash
scmethyl build-matrix \
    --gene-meth-dir /path/to/CG_gene_meth \
    --gene-names-file /path/to/gene_names.txt \
    --output-file /path/to/output/cell_mCG_matrix \
    --context CG \
    --cleanup
```

#### Step 7: Build AnnData object
```bash
scmethyl build-adata \
    --cg-matrix /path/to/cell_mCG_matrix \
    --ch-matrix /path/to/cell_mCH_matrix \
    --metadata /path/to/cell_metadata.csv \
    --output /path/to/output/adata.h5ad
```

## Python API Usage

You can also use scMethylSeq as a Python library:

```python
from scmethylseq import (
    AllcProcessor,
    CovFileGenerator,
    GlobalMethylationCalculator,
    QualityControl,
    GeneBodyMethylation,
    MethylationMatrixBuilder,
    PCAAnalyzer,
    UMAPAnalyzer,
    LeidenClustering,
    AnnDataBuilder
)

# Example: Process allc files
processor = AllcProcessor(
    allc_dir="/path/to/allc",
    output_dir="/path/to/output",
    chrom_sizes="/path/to/chrom.sizes",
    n_jobs=10
)
processor.process_all()

# Example: Build methylation matrix
builder = MethylationMatrixBuilder(
    gene_meth_dir="/path/to/gene_meth",
    gene_names_file="/path/to/genes.txt",
    output_file="/path/to/matrix.txt",
    context="CG"
)
matrix = builder.build_matrix()

# Example: Build AnnData
adata_builder = AnnDataBuilder(
    cg_matrix_file="/path/to/CG_matrix",
    ch_matrix_file="/path/to/CH_matrix",
    metadata_file="/path/to/metadata.csv",
    output_file="/path/to/adata.h5ad"
)
adata = adata_builder.build_and_save()
```

## Configuration File Format

The configuration file (YAML format) controls all pipeline parameters:

```yaml
# Enable/disable pipeline steps
run_allc_processing: true
run_cov_conversion: true
run_global_meth: true
run_qc: true
run_gene_meth: true
run_matrix_building: true
run_adata_building: true

# System resources
n_jobs: 20

# Input paths
allc_dir: "/path/to/allc"
chrom_sizes: "/path/to/hg38.chrom.sizes"
gene_bed: "/path/to/genebody.bed"

# QC parameters
qc_params:
  cg_min: 0.5
  cg_sites_min: 100000
  ch_max: 0.1
```

See `config_example.yaml` for a complete example.

## Input File Formats

### allc files
Standard allc format (gzipped):
```
chr1	3000827	+	CGN	0	1	1
chr1	3001007	+	CGN	0	1	1
```

### Gene BED file
BED format with gene coordinates (can include 2kb extensions):
```
chr1	1000	5000	GENE1
chr1	10000	15000	GENE2
```

### Metadata CSV
Must include a 'cell' column and optionally UMAP coordinates:
```
cell,leiden_cluster,UMAP1,UMAP2,batch
cell_001,0,1.23,-2.45,batch1
cell_002,1,-0.34,3.21,batch1
```

## Output Files

### Cell-level metrics
- `cell_raw.csv`: All cells with QC metrics
- `cell_filtered.csv`: Filtered cells passing QC

### Methylation matrices
- `cell_mCG_matrix`: Genes × Cells CG methylation matrix
- `cell_mCH_matrix`: Genes × Cells CH methylation matrix

### AnnData object
- `adata.h5ad`: Complete AnnData object with:
  - `X`: CG methylation matrix
  - `obsm['X_cg']`: CG methylation
  - `obsm['X_ch']`: CH methylation
  - `obsm['X_umap']`: UMAP coordinates
  - `obs`: Cell metadata

## Performance Optimization

The toolkit uses multiple strategies for fast processing:

1. **Parallel processing**: All I/O and computation-intensive steps use multiprocessing
2. **Efficient file I/O**: Uses system tools (bedtools, paste) where possible
3. **Memory management**: Streaming processing for large files
4. **Configurable resources**: Adjust `n_jobs` based on your system

### Recommended settings:

- **Small dataset** (<1000 cells): `n_jobs=10`
- **Medium dataset** (1000-10000 cells): `n_jobs=20`
- **Large dataset** (>10000 cells): `n_jobs=30-50`

## Troubleshooting

### Memory issues
If you encounter memory errors:
1. Reduce `n_jobs`
2. Process CG and CH contexts separately
3. Use `cleanup_intermediate: true` to remove temporary files

### bedtools not found
Install bedtools on your system:
```bash
# Ubuntu/Debian
sudo apt-get install bedtools

# macOS
brew install bedtools

# Conda
conda install -c bioconda bedtools
```

## Citation

If you use scMethylSeq in your research, please cite:

```
[Your citation information]
```

## License

MIT License

## Contact

For questions and support, please open an issue on GitHub or contact:
[Your contact information]

## Acknowledgments

This toolkit was developed based on analysis workflows for single-cell DNA methylation data.
