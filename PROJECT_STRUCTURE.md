# Project Structure

```
scMethylSeq/
├── README.md                      # Main documentation
├── QUICKSTART.md                  # Quick start guide (15 min tutorial)
├── INSTALL.md                     # Detailed installation guide
├── LICENSE                        # MIT License
├── .gitignore                     # Git ignore patterns
├── setup.py                       # Package installation script
├── config_example.yaml            # Example configuration file
│
├── scmethylseq/                   # Main package directory
│   ├── __init__.py               # Package initialization
│   ├── preprocessing.py          # Step 1-2: allc/cov processing, QC
│   ├── matrix_builder.py         # Step 3: Gene body methylation, matrix building
│   ├── clustering.py             # PCA, UMAP, Leiden clustering
│   ├── adata_builder.py          # Final AnnData object construction
│   ├── cli.py                    # Command-line interface
│   └── utils.py                  # Utility functions
│
└── examples/                      # Example scripts and notebooks
    └── example_usage.py          # Complete workflow example

```

## Module Descriptions

### Core Modules

#### preprocessing.py
Classes for early-stage data processing:
- `AllcProcessor`: Extract CG/CH contexts from allc files
- `CovFileGenerator`: Convert allc to cov format
- `GlobalMethylationCalculator`: Calculate global CG/CH levels
- `QualityControl`: Filter cells based on QC metrics

**Key Features:**
- Parallel processing with configurable workers
- Automatic file format handling
- QC visualization

#### matrix_builder.py
Classes for building methylation matrices:
- `GeneBodyMethylation`: Calculate methylation for gene bodies
- `MethylationMatrixBuilder`: Build cell × gene matrices
- `VMRMethylationExtractor`: Extract VMR methylation (planned)
- `TSSMethylationExtractor`: Extract TSS methylation (planned)

**Key Features:**
- Bedtools integration for efficient genomic operations
- Memory-efficient streaming processing
- Automatic cleanup of intermediate files

#### clustering.py
Classes for dimensionality reduction and clustering:
- `IterativePCA`: PCA with missing value imputation
- `PCAAnalyzer`: Standard PCA analysis
- `UMAPAnalyzer`: UMAP dimensionality reduction
- `LeidenClustering`: Graph-based clustering

**Key Features:**
- Handles missing values in methylation data
- Multiple resolution parameter scanning
- Integration with scanpy workflows

#### adata_builder.py
Classes for AnnData construction:
- `AnnDataBuilder`: Build complete AnnData from matrices
- Helper functions for adding PCA, UMAP, clustering results
- Validation and integration utilities

**Key Features:**
- Automatic metadata merging
- Multiple methylation context support (CG, CH)
- Scanpy-compatible format

#### cli.py
Command-line interface with commands:
- `process-allc`: Process allc files
- `convert-to-cov`: Convert to cov format
- `calc-global-meth`: Calculate global methylation
- `filter-cells`: QC filtering
- `calc-gene-meth`: Gene body methylation
- `build-matrix`: Build methylation matrix
- `build-adata`: Create AnnData object
- `run-pipeline`: Run complete pipeline from config

#### utils.py
Utility functions:
- File validation and management
- Dependency checking
- Memory estimation
- Progress tracking
- Statistics computation

### Configuration

#### config_example.yaml
Template configuration file with:
- Pipeline step toggles
- Input/output paths
- QC parameters
- Resource allocation
- All configurable options

### Documentation

#### README.md
- Project overview
- Feature list
- Installation instructions
- Usage examples
- API documentation

#### QUICKSTART.md
Three quick-start scenarios:
1. From methylation matrices (5 min)
2. From cov files (30 min)
3. Complete pipeline with config (easiest)

#### INSTALL.md
Comprehensive installation guide:
- System requirements
- Dependency installation
- Testing procedures
- Troubleshooting
- Performance benchmarks

### Examples

#### example_usage.py
Complete workflow demonstration:
- Loading matrices
- PCA analysis
- UMAP embedding
- Leiden clustering
- Metadata creation
- AnnData construction
- Visualization

## Data Flow

```
Raw Data (allc.gz files)
    ↓
[AllcProcessor]
    ↓
CG/CH allc files
    ↓
[CovFileGenerator]
    ↓
Cov files (.cov)
    ↓
[GlobalMethylationCalculator]
    ↓
Global methylation metrics
    ↓
[QualityControl]
    ↓
Filtered cell list
    ↓
[GeneBodyMethylation]
    ↓
Gene methylation per cell
    ↓
[MethylationMatrixBuilder]
    ↓
Cell × Gene matrices (CG, CH)
    ↓
[PCA/UMAP/Leiden]
    ↓
Clustering & embeddings
    ↓
[AnnDataBuilder]
    ↓
Final AnnData (.h5ad)
```

## Key Design Principles

### 1. Modularity
Each processing step is a separate class that can be used independently or as part of the pipeline.

### 2. Parallel Processing
All computationally intensive steps support multi-threading/multi-processing via `n_jobs` parameter.

### 3. Memory Efficiency
- Streaming file processing
- Automatic cleanup of intermediate files
- Configurable memory usage

### 4. Flexibility
- Run complete pipeline or individual steps
- Python API and CLI interface
- Compatible with scanpy ecosystem

### 5. Robustness
- Input validation
- Comprehensive error handling
- Progress tracking for long-running jobs
- Automatic dependency checking

## Extending the Package

### Adding New Features

1. **New methylation context**
   - Add to `preprocessing.py` (extraction)
   - Update `matrix_builder.py` (calculation)
   - Modify `adata_builder.py` (storage)

2. **New genomic feature**
   - Create new class in `matrix_builder.py`
   - Follow pattern of `GeneBodyMethylation`
   - Add CLI command in `cli.py`

3. **New clustering method**
   - Add class to `clustering.py`
   - Implement `fit_transform` interface
   - Update `adata_builder.py` integration

### Testing

Create test files in `tests/` directory:
```python
# tests/test_preprocessing.py
import pytest
from scmethylseq import GlobalMethylationCalculator

def test_global_methylation():
    # Test code here
    pass
```

Run tests:
```bash
pytest tests/
```

## Version History

- **0.1.0** (Current): Initial release
  - Basic pipeline implementation
  - CG/CH methylation support
  - PCA, UMAP, Leiden clustering
  - CLI and Python API

## Future Enhancements

Planned features:
- VMR (Variable Methylated Regions) analysis
- TSS (Transcription Start Site) methylation
- Differential methylation analysis
- Integration with RNA-seq data
- GPU acceleration for large datasets
- Web-based visualization interface
