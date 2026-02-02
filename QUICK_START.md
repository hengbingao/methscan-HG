# MethSCAn2 - Complete Package 🎉

## 📦 What's Included

This is a **complete, ready-to-use** Python package for single-cell DNA methylation analysis!

### Package Contents:
```
methscan2/
├── methscan2/              # Main package
│   ├── core/              # Data structures (MethylationData)
│   ├── preprocessing/     # Data loading, QC, filtering
│   ├── tools/             # Analysis tools (VMR, PCA, UMAP, clustering)
│   └── plotting/          # Visualization functions
├── examples/              # Example scripts
├── tests/                 # Unit tests
├── setup.py              # Installation configuration
├── README.md             # Documentation
├── LICENSE               # MIT license
└── .github/workflows/    # CI/CD configuration
```

## 🚀 Quick Installation

### Option 1: Install from extracted folder (Immediate Use)

```bash
# Extract the archive
tar -xzf methscan2-project.tar.gz
cd methscan2-project

# Install
pip install -e .

# Test installation
python test_installation.py
```

### Option 2: Upload to GitHub (Recommended)

Follow instructions in `GITHUB_SETUP.md` to:
1. Create a GitHub repository
2. Push the code
3. Allow others to install with: `pip install git+https://github.com/yourusername/methscan2.git`

## ⚡ Quick Start (3 minutes!)

### Step 1: Install
```bash
cd methscan2-project
pip install -e .
```

### Step 2: Generate test data
```bash
python examples/generate_test_data.py
```

### Step 3: Run analysis
```bash
python examples/complete_analysis.py
```

### Step 4: Check results
```bash
ls results/
# Should see: analyzed_data.h5ad, *.pdf, clusters.csv, vmrs.bed
```

## 📖 Basic Usage

```python
import methscan2 as ms2

# Load your .cov files
mdata = ms2.read_cov_files('path/to/*.cov')

# Quality control
ms2.pp.calculate_qc_metrics(mdata)
ms2.pp.filter_cells(mdata, min_sites=50000)

# Detect VMRs
ms2.pp.smooth_methylation(mdata)
ms2.tl.detect_vmr(mdata)
ms2.tl.quantify_methylation(mdata)

# Analyze
ms2.tl.run_pca(mdata)
ms2.tl.run_umap(mdata)
ms2.tl.run_leiden(mdata)

# Visualize
ms2.pl.umap(mdata, color='leiden', save='umap.pdf')
ms2.pl.heatmap(mdata, groupby='leiden', save='heatmap.pdf')

# Save
mdata.write('results.h5ad')
```

## 🔧 System Requirements

- **Python**: ≥ 3.8
- **OS**: Linux, macOS, or Windows
- **RAM**: 4GB minimum (8GB+ recommended for large datasets)
- **Dependencies**: Auto-installed with pip

## 📊 Input Data Formats

MethSCAn2 accepts:

1. **Bismark .cov files** (Recommended)
   ```
   chr1  10468  10469  100  1  0
   chr1  10470  10471  100  0  1
   ```

2. **ALLC files** (methylpy output)
   ```
   chr1  3000827  +  CCC  0  1  1
   ```

3. **bedGraph files**
   ```
   chr1  10468  10469  1.0
   ```

## 🎯 Key Features

✅ **No R dependency** - Pure Python implementation  
✅ **Fast** - Parallel processing with joblib  
✅ **Scalable** - Handles 1000s of cells  
✅ **Compatible** - Integrates with scanpy/AnnData ecosystem  
✅ **Complete** - Data loading → Analysis → Visualization  
✅ **Well-tested** - Unit tests and CI/CD  
✅ **Documented** - Clear examples and docstrings  

## 📚 Documentation

- **README.md** - Overview and quick start
- **GITHUB_SETUP.md** - How to upload to GitHub
- **examples/complete_analysis.py** - Full analysis example
- **examples/generate_test_data.py** - Create synthetic data
- **test_installation.py** - Verify installation

## 🤝 Contributing

To contribute:
1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests
5. Submit a pull request

## 🐛 Troubleshooting

### Installation fails?
```bash
# Upgrade pip first
pip install --upgrade pip

# Install with verbose output
pip install -e . -v
```

### Missing dependencies?
```bash
# Install manually
pip install numpy pandas scipy matplotlib seaborn
pip install anndata scikit-learn umap-learn leidenalg python-igraph
```

### Import errors?
```bash
# Make sure you're in the right directory
cd methscan2-project
pip install -e .

# Test
python -c "import methscan2; print(methscan2.__version__)"
```

## 📧 Support

- **Issues**: Open a GitHub issue
- **Questions**: Check examples/ directory first
- **Feature requests**: Open a GitHub discussion

## 📄 License

MIT License - See LICENSE file

## 🙏 Acknowledgments

Based on [MethSCAn](https://github.com/anders-biostat/MethSCAn) by Kremer et al. (2024).  
Inspired by [snapATAC2](https://github.com/kaizhang/SnapATAC2) architecture.

## ⭐ Citation

If you use MethSCAn2 in your research, please cite:

```bibtex
@software{methscan2_2024,
  title={MethSCAn2: A Python toolkit for single-cell DNA methylation analysis},
  author={Your Name},
  year={2024},
  url={https://github.com/yourusername/methscan2}
}
```

---

**Ready to analyze single-cell methylation data?** 🚀

Run `python test_installation.py` to get started!
