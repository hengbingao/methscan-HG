# MethSCAn2

**Single-cell DNA Methylation Analysis Toolkit**

A comprehensive Python package for analyzing single-cell DNA methylation data, reimagined as a pure Python implementation with no R dependencies.

## 🚀 Quick Start

### Installation

```bash
# Install from GitHub
pip install git+https://github.com/yourusername/methscan2.git

# Or install from source
git clone https://github.com/yourusername/methscan2.git
cd methscan2
pip install -e .
```

### Basic Usage

```python
import methscan2 as ms2

# 1. Load data
mdata = ms2.read_cov_files('data/*.cov')

# 2. Quality control
ms2.pp.calculate_qc_metrics(mdata)
ms2.pp.filter_cells(mdata, min_sites=50000)

# 3. Detect VMRs
ms2.pp.smooth_methylation(mdata)
ms2.tl.detect_vmr(mdata)

# 4. Quantify methylation
ms2.tl.quantify_methylation(mdata)

# 5. Dimensionality reduction
ms2.tl.run_pca(mdata, n_comps=50)
ms2.tl.run_umap(mdata)

# 6. Clustering
ms2.tl.run_leiden(mdata, resolution=0.8)

# 7. Visualization
ms2.pl.umap(mdata, color='leiden', save='umap.pdf')
ms2.pl.heatmap(mdata, groupby='leiden', save='heatmap.pdf')

# 8. Save results
mdata.write('results.h5ad')
```

## 📦 Features

- **Data I/O**: Read Bismark cov, ALLC, and bedGraph formats
- **Quality Control**: Comprehensive QC metrics and filtering
- **VMR Detection**: Identify variably methylated regions
- **Dimensionality Reduction**: PCA and UMAP
- **Clustering**: Leiden and Louvain algorithms
- **Visualization**: Beautiful matplotlib/seaborn plots
- **AnnData Integration**: Compatible with scanpy ecosystem

## 📊 Example

See `examples/` directory for complete analysis examples.

## 🛠️ Requirements

- Python ≥ 3.8
- numpy, scipy, pandas
- anndata
- scikit-learn
- matplotlib, seaborn
- umap-learn
- leidenalg, python-igraph

## 📖 Documentation

Coming soon at [methscan2.readthedocs.io](https://methscan2.readthedocs.io)

## 🤝 Contributing

Contributions welcome! Please open an issue or pull request.

## 📄 License

MIT License

## 🙏 Acknowledgments

Based on the original [MethSCAn](https://github.com/anders-biostat/MethSCAn) by Kremer et al. (Nature Methods, 2024).

## 📧 Contact

For questions and support, please open a GitHub issue.
