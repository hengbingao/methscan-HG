# MethSCAn2 v0.2.0 - Multi-Context Analysis Update

## 🎉 Major Feature Release

MethSCAn2 now supports **comprehensive multi-context methylation analysis** with genomic region annotation!

---

## ✨ What's New in v0.2.0

### 1. Multi-Context Data Structure

**Store all methylation contexts in one object:**

```python
from methscan2.core.multi_context_data import MultiContextMethylationData

# Single object with multiple contexts
mdata = MultiContextMethylationData(...)
mdata.available_contexts  # ['CG', 'CH', 'CHG', 'CHH', 'all']

# Access any context
cg_rate = mdata.get_context('CG', 'rate')
ch_rate = mdata.get_context('CH', 'rate')

# Switch active context for analysis
mdata.set_active_context('CG')  # Run analysis on CG
mdata.set_active_context('CH')  # Run analysis on CH
```

### 2. Genomic Region Annotation

**Calculate methylation across genomic features:**

```python
from methscan2.preprocessing.genomic_annotation import GenomicRegionAnnotator

# Load gene annotations
annotator = GenomicRegionAnnotator(gtf_file='genes.gtf')

# Define regions
tss = annotator.define_tss_regions(upstream=2000, downstream=500)
gene_body = annotator.define_gene_body_regions()
promoters = annotator.define_promoter_regions()
upstream_downstream = annotator.define_upstream_downstream_regions()
```

### 3. Direct ALLC to Multi-Context Analysis

**One function creates complete multi-context data:**

```python
from methscan2.core.multi_context_data import create_multi_context_data_from_allc

mdata = create_multi_context_data_from_allc(
    allc_files='raw/allc/*.allc.tsv.gz',
    regions=tss_regions_df,  # TSS, gene body, etc.
    contexts=['CG', 'CH', 'CHG', 'CHH', 'all'],
    region_type='tss',
    n_jobs=100
)
```

### 4. Context-Specific Analysis

**Run separate analyses for each context:**

```python
# Analyze CG
mdata.set_active_context('CG')
ms2.tl.run_pca(mdata)
ms2.tl.run_umap(mdata)
ms2.tl.run_leiden(mdata, key_added='leiden_CG')

# Analyze CH
mdata.set_active_context('CH')
ms2.tl.run_pca(mdata)
ms2.tl.run_umap(mdata)
ms2.tl.run_leiden(mdata, key_added='leiden_CH')

# Compare contexts
corr = mdata.calculate_context_correlations()
```

---

## 🏗️ Architecture

### Data Storage

```
MultiContextMethylationData
├── X                    # Active context (switchable)
├── obs                  # Cell metadata
├── var                  # Region metadata (chr, start, end)
├── layers
│   ├── CG_rate         # CG methylation rates
│   ├── CG_met          # CG methylated counts
│   ├── CG_total        # CG total coverage
│   ├── CH_rate         # CH methylation rates
│   ├── CH_met          # CH methylated counts
│   ├── CH_total        # CH total coverage
│   ├── all_rate        # All mC rates
│   ├── all_met         # All mC counts
│   └── all_total       # All mC coverage
├── obsm
│   ├── X_pca_CG        # CG PCA results
│   ├── X_umap_CG       # CG UMAP results
│   ├── X_pca_CH        # CH PCA results
│   ├── X_umap_CH       # CH UMAP results
│   └── ...
└── uns
    ├── contexts        # ['CG', 'CH', 'all']
    ├── active_context  # Currently active context
    └── region_type     # 'tss', 'gene_body', etc.
```

---

## 🚀 Complete Workflow

### Workflow 1: Quick Analysis

```python
import methscan2 as ms2

# Load with multiple contexts
mdata = ms2.create_multi_context_data_from_allc(
    allc_files='raw/allc/*.allc.tsv.gz',
    regions=regions_df,
    contexts=['CG', 'CH', 'all'],
    n_jobs=100
)

# Analyze each context
for ctx in ['CG', 'CH']:
    mdata.set_active_context(ctx)
    ms2.tl.run_pca(mdata)
    ms2.tl.run_umap(mdata)
    ms2.pl.umap(mdata, color=f'leiden_{ctx}', save=f'umap_{ctx}.pdf')

# Compare contexts
mdata.plot_context_comparison(save='context_comparison.pdf')
```

### Workflow 2: Complete Pipeline (Command Line)

```bash
python examples/multi_context_analysis.py \
    --allc raw/allc \
    --gtf gencode.v38.gtf \
    --output results \
    --contexts CG CH CHG CHH all \
    --analysis-contexts CG CH \
    --region-types tss gene_body upstream downstream \
    --jobs 100
```

**Output:**
```
results/
├── multi_context_data.h5ad           # Main object (all contexts)
├── mdata_tss.h5ad                    # TSS regions
├── mdata_gene_body.h5ad              # Gene body regions
├── mdata_upstream.h5ad               # Upstream regions
├── mdata_downstream.h5ad             # Downstream regions
├── cluster_assignments.csv           # Clusters for all contexts
├── methylation_statistics.csv        # Methylation stats
├── context_correlations.csv          # Context correlations
└── *.pdf                             # Visualizations
```

---

## 📊 Use Cases

### Use Case 1: Compare CG and CH Methylation

```python
# CG typically shows cell type identity
# CH shows finer cell state differences

# Load data with both contexts
mdata = ms2.create_multi_context_data_from_allc(
    allc_files='brain_cells/*.allc.tsv.gz',
    regions=tss_df,
    contexts=['CG', 'CH']
)

# Cluster on CG
mdata.set_active_context('CG')
ms2.tl.run_leiden(mdata, key_added='leiden_CG')
# → Broad cell types (neurons, glia, etc.)

# Cluster on CH
mdata.set_active_context('CH')
ms2.tl.run_leiden(mdata, key_added='leiden_CH')
# → Finer subtypes within cell types

# Compare
from sklearn.metrics import adjusted_rand_score
ari = adjusted_rand_score(
    mdata.obs['leiden_CG'],
    mdata.obs['leiden_CH']
)
print(f"ARI: {ari:.3f}")  # How well do contexts agree?
```

### Use Case 2: Analyze Multiple Genomic Regions

```python
annotator = GenomicRegionAnnotator(gtf_file='genes.gtf')

# Process multiple region types
region_types = {
    'tss': annotator.define_tss_regions(),
    'gene_body': annotator.define_gene_body_regions(),
    'promoters': annotator.define_promoter_regions()
}

# Create dataset for each
mdata_dict = {}
for name, regions in region_types.items():
    mdata_dict[name] = ms2.create_multi_context_data_from_allc(
        allc_files='*.allc.tsv.gz',
        regions=regions.df[['Chromosome', 'Start', 'End']].rename(
            columns={'Chromosome': 'chr', 'Start': 'start', 'End': 'end'}
        ),
        contexts=['CG', 'CH'],
        region_type=name
    )

# Compare clustering across regions
for name, mdata in mdata_dict.items():
    mdata.set_active_context('CG')
    ms2.tl.run_leiden(mdata)
    print(f"{name}: {mdata.obs['leiden'].nunique()} clusters")
```

### Use Case 3: Integration with scRNA-seq

```python
# Load methylation (CG context)
mdata = ms2.create_multi_context_data_from_allc(...)
mdata.set_active_context('CG')

# Load RNA-seq
import scanpy as sc
rna = sc.read_h5ad('scrna.h5ad')

# Match cells
common = set(mdata.obs_names) & set(rna.obs_names)
mdata_sub = mdata[list(common), :].copy()
rna_sub = rna[list(common), :].copy()

# Compare cell types
ari = adjusted_rand_score(
    mdata_sub.obs['leiden_CG'],
    rna_sub.obs['leiden']
)
```

---

## 🎯 Key Features

### 1. Unified Storage
✅ All contexts in one AnnData object  
✅ No data duplication  
✅ Easy context switching  
✅ Consistent metadata across contexts  

### 2. Genomic Annotation
✅ TSS regions  
✅ Gene body  
✅ Upstream/Downstream  
✅ Promoters  
✅ Custom regions (BED files)  

### 3. Flexible Analysis
✅ Run analysis on any context  
✅ Compare contexts directly  
✅ Context-specific clustering  
✅ Context correlation analysis  

### 4. High Performance
✅ Parallel processing (100+ jobs)  
✅ Memory efficient  
✅ Handles large datasets  
✅ Fast region quantification  

---

## 📚 Documentation

### New Documentation Files

1. **docs/MULTI_CONTEXT_GUIDE.md**
   - Complete guide to multi-context analysis
   - Step-by-step workflows
   - Advanced usage examples

2. **examples/multi_context_analysis.py**
   - Complete working example
   - Command-line tool
   - Fully documented

### Updated Files

- `README.md` - Added multi-context features
- `methscan2/__init__.py` - New imports
- `methscan2/core/__init__.py` - Multi-context classes
- `methscan2/preprocessing/__init__.py` - Genomic annotation

---

## 🔧 New Modules

### Core Modules

1. **methscan2/core/multi_context_data.py** (500+ lines)
   - `MultiContextMethylationData` class
   - `create_multi_context_data_from_allc()` function
   - Context switching and comparison

2. **methscan2/preprocessing/genomic_annotation.py** (400+ lines)
   - `GenomicRegionAnnotator` class
   - Region definition functions
   - Methylation quantification across regions

### Example Scripts

3. **examples/multi_context_analysis.py** (400+ lines)
   - Complete pipeline
   - Command-line interface
   - Multiple region types

---

## 📦 Installation

```bash
# Update to v0.2.0
cd methscan2-project
git pull  # If using Git
pip install -e .

# Install new dependencies
pip install pyranges
```

---

## 🎓 Quick Start

### Minimal Example

```python
import methscan2 as ms2

# Create multi-context data
mdata = ms2.create_multi_context_data_from_allc(
    allc_files='*.allc.tsv.gz',
    regions=None,  # Use genome-wide bins
    contexts=['CG', 'CH', 'all'],
    n_jobs=100
)

# Analyze CG
mdata.set_active_context('CG')
ms2.tl.run_umap(mdata)
ms2.pl.umap(mdata, save='cg.pdf')

# Analyze CH
mdata.set_active_context('CH')
ms2.tl.run_umap(mdata)
ms2.pl.umap(mdata, save='ch.pdf')

# Compare
corr = mdata.calculate_context_correlations()
print(corr)
```

---

## 🔄 Migration from v0.1.0

### Old Way (v0.1.0)
```python
# Separate analysis for each context
mdata_cg = ms2.read_cov_files('CG_cov/*.cov')
mdata_ch = ms2.read_cov_files('CH_cov/*.cov')

# Analyze separately
ms2.tl.run_umap(mdata_cg)
ms2.tl.run_umap(mdata_ch)

# Hard to compare
```

### New Way (v0.2.0)
```python
# Unified multi-context object
mdata = ms2.create_multi_context_data_from_allc(
    allc_files='*.allc.tsv.gz',
    contexts=['CG', 'CH']
)

# Easy to switch and compare
mdata.set_active_context('CG')
ms2.tl.run_umap(mdata)

mdata.set_active_context('CH')
ms2.tl.run_umap(mdata)

mdata.calculate_context_correlations()
```

---

## ✅ Testing

```bash
# Run tests
pytest tests/ -v

# Test multi-context features
pytest tests/test_multi_context.py -v
```

---

## 🎊 Summary

**MethSCAn2 v0.2.0** provides:

✅ **Multi-context data structure** - All contexts in one object  
✅ **Genomic region annotation** - TSS, gene body, promoters, custom  
✅ **Context switching** - Analyze any context easily  
✅ **Context comparison** - Direct correlation analysis  
✅ **High performance** - 100+ parallel jobs  
✅ **Complete integration** - Works with existing tools  

**From ALLC to multi-context analysis in one step!** 🚀

---

## 📞 Support

- **Documentation**: `docs/MULTI_CONTEXT_GUIDE.md`
- **Examples**: `examples/multi_context_analysis.py`
- **Issues**: GitHub Issues
- **Questions**: GitHub Discussions

**Happy analyzing!** 🎉
