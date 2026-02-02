# Multi-Context Methylation Analysis Guide

## 🎯 Overview

MethSCAn2 now supports **comprehensive multi-context methylation analysis**, allowing you to:

✅ Store **CG, CH, CHG, CHH, and all** methylation contexts in a single AnnData object  
✅ Calculate methylation across **TSS, gene body, upstream/downstream** regions  
✅ Run **separate analyses** for each context (UMAP, clustering)  
✅ **Compare contexts** within the same cells  
✅ **Switch contexts** dynamically for different analyses  

---

## 🏗️ Architecture

### Multi-Layer Storage

All methylation contexts are stored in separate layers:

```python
mdata = MultiContextMethylationData(...)

# Data structure:
mdata.X                  # Active context (default: CG)
mdata.layers['CG_rate']  # CG methylation rates
mdata.layers['CG_met']   # CG methylated counts
mdata.layers['CG_total'] # CG total coverage
mdata.layers['CH_rate']  # CH methylation rates
mdata.layers['CH_met']   # CH methylated counts
mdata.layers['CH_total'] # CH total coverage
mdata.layers['all_rate'] # All mC rates
# etc.
```

### Genomic Region Support

Calculate methylation across multiple genomic features:

- **TSS** (Transcription Start Sites)
- **Gene body** (from TSS to TES)
- **Upstream** regions (e.g., -5kb from TSS)
- **Downstream** regions (e.g., +5kb from TES)
- **Promoters** (TSS ± defined window)
- **Custom regions** (from BED files)

---

## 🚀 Quick Start

### Option 1: Complete Automated Pipeline

```bash
python examples/multi_context_analysis.py \
    --allc raw/allc \
    --gtf genes.gtf \
    --output results \
    --contexts CG CH all \
    --jobs 100
```

### Option 2: Python Script

```python
import methscan2 as ms2
from methscan2.core.multi_context_data import create_multi_context_data_from_allc
from methscan2.preprocessing.genomic_annotation import GenomicRegionAnnotator

# 1. Load genomic annotations
annotator = GenomicRegionAnnotator(gtf_file='genes.gtf')
tss_regions = annotator.define_tss_regions(upstream=2000, downstream=500)

# Convert to DataFrame
regions_df = tss_regions.df[['Chromosome', 'Start', 'End']]
regions_df.columns = ['chr', 'start', 'end']

# 2. Create multi-context data from ALLC files
mdata = create_multi_context_data_from_allc(
    allc_files='raw/allc/*.allc.tsv.gz',
    regions=regions_df,
    contexts=['CG', 'CH', 'all'],
    region_type='tss',
    n_jobs=100
)

# Result: MultiContextMethylationData with all contexts in layers
print(mdata.available_contexts)  # ['CG', 'CH', 'all']
```

---

## 📚 Detailed Workflow

### Step 1: Define Genomic Regions

```python
from methscan2.preprocessing.genomic_annotation import GenomicRegionAnnotator

# Load gene annotations
annotator = GenomicRegionAnnotator(gtf_file='gencode.v38.gtf', genome='hg38')

# Define TSS regions (±2kb)
tss_regions = annotator.define_tss_regions(
    upstream=2000,
    downstream=500
)

# Define gene body regions
gene_body = annotator.define_gene_body_regions()

# Define upstream/downstream
up_down = annotator.define_upstream_downstream_regions(
    upstream=5000,
    downstream=5000
)

# Define promoters
promoters = annotator.define_promoter_regions(
    upstream=2000,
    downstream=500
)
```

### Step 2: Process ALLC Files for Each Region Type

```python
from methscan2.core.multi_context_data import create_multi_context_data_from_allc

# Process TSS regions
mdata_tss = create_multi_context_data_from_allc(
    allc_files='raw/allc/*.allc.tsv.gz',
    regions=tss_regions.df[['Chromosome', 'Start', 'End']].rename(
        columns={'Chromosome': 'chr', 'Start': 'start', 'End': 'end'}
    ),
    contexts=['CG', 'CH', 'CHG', 'CHH', 'all'],
    region_type='tss',
    min_coverage=1,
    n_jobs=100
)

# Process gene body
mdata_gb = create_multi_context_data_from_allc(
    allc_files='raw/allc/*.allc.tsv.gz',
    regions=gene_body.df[['Chromosome', 'Start', 'End']].rename(
        columns={'Chromosome': 'chr', 'Start': 'start', 'End': 'end'}
    ),
    contexts=['CG', 'CH', 'all'],
    region_type='gene_body',
    n_jobs=100
)
```

### Step 3: Quality Control per Context

```python
# Calculate QC metrics for each context
for ctx in mdata_tss.available_contexts:
    rate = mdata_tss.get_context(ctx, 'rate')
    total = mdata_tss.get_context(ctx, 'total')
    
    # Add metrics to obs
    mdata_tss.obs[f'{ctx}_n_sites'] = (total > 0).sum(axis=1)
    mdata_tss.obs[f'{ctx}_mean_cov'] = np.nanmean(
        np.where(total > 0, total, np.nan), axis=1
    )
    mdata_tss.obs[f'{ctx}_mean_meth'] = np.nanmean(rate, axis=1)

# Filter cells based on CG context
keep = (
    (mdata_tss.obs['CG_n_sites'] >= 1000) &
    (mdata_tss.obs['CG_mean_meth'] >= 0.3) &
    (mdata_tss.obs['CG_mean_meth'] <= 0.8)
)
mdata_tss = mdata_tss[keep, :].copy()
```

### Step 4: Run Analysis on Each Context

```python
# Analyze CG context
mdata_tss.set_active_context('CG')
ms2.tl.run_pca(mdata_tss, n_comps=50)
ms2.tl.run_umap(mdata_tss)
ms2.tl.run_leiden(mdata_tss, key_added='leiden_CG')

# Save CG results
mdata_tss.obsm['X_pca_CG'] = mdata_tss.obsm['X_pca'].copy()
mdata_tss.obsm['X_umap_CG'] = mdata_tss.obsm['X_umap'].copy()

# Analyze CH context
mdata_tss.set_active_context('CH')
ms2.tl.run_pca(mdata_tss, n_comps=50)
ms2.tl.run_umap(mdata_tss)
ms2.tl.run_leiden(mdata_tss, key_added='leiden_CH')

# Save CH results
mdata_tss.obsm['X_pca_CH'] = mdata_tss.obsm['X_pca'].copy()
mdata_tss.obsm['X_umap_CH'] = mdata_tss.obsm['X_umap'].copy()

# Analyze all mC
mdata_tss.set_active_context('all')
ms2.tl.run_pca(mdata_tss, n_comps=50)
ms2.tl.run_umap(mdata_tss)
ms2.tl.run_leiden(mdata_tss, key_added='leiden_all')

# Save all results
mdata_tss.obsm['X_pca_all'] = mdata_tss.obsm['X_pca'].copy()
mdata_tss.obsm['X_umap_all'] = mdata_tss.obsm['X_umap'].copy()
```

### Step 5: Compare Contexts

```python
# Calculate correlations between contexts
corr_df = mdata_tss.calculate_context_correlations()
print(corr_df)
#   context1 context2  mean_correlation  std_correlation
# 0       CG       CG          1.000000         0.000000
# 1       CG       CH          0.456789         0.123456
# 2       CH       CH          1.000000         0.000000

# Compare methylation distributions
mdata_tss.plot_context_comparison(
    contexts=['CG', 'CH', 'all'],
    metric='mean',
    save='context_comparison.pdf'
)

# Compare cluster assignments
from sklearn.metrics import adjusted_rand_score

ari = adjusted_rand_score(
    mdata_tss.obs['leiden_CG'],
    mdata_tss.obs['leiden_CH']
)
print(f"ARI between CG and CH clustering: {ari:.3f}")
```

### Step 6: Visualization

```python
# UMAP for each context
for ctx in ['CG', 'CH', 'all']:
    ms2.pl.umap(
        mdata_tss,
        color=f'leiden_{ctx}',
        use_rep=f'X_umap_{ctx}',
        save=f'umap_{ctx}.pdf'
    )

# Multi-panel comparison
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 3, figsize=(18, 5))

for ax, ctx in zip(axes, ['CG', 'CH', 'all']):
    coords = mdata_tss.obsm[f'X_umap_{ctx}']
    clusters = mdata_tss.obs[f'leiden_{ctx}']
    
    ax.scatter(
        coords[:, 0], coords[:, 1],
        c=clusters.cat.codes, cmap='tab20',
        s=10, alpha=0.8
    )
    ax.set_title(f'{ctx} Context')
    ax.set_xlabel('UMAP1')
    ax.set_ylabel('UMAP2')

plt.tight_layout()
plt.savefig('umap_all_contexts.pdf')
```

### Step 7: Save Results

```python
# Save multi-context object
mdata_tss.write('results/multi_context_tss.h5ad')

# Export cluster assignments
clusters = mdata_tss.obs[['leiden_CG', 'leiden_CH', 'leiden_all']]
clusters.to_csv('results/clusters_all_contexts.csv')

# Export methylation statistics
stats = mdata_tss.obs[[
    'CG_mean_meth', 'CH_mean_meth', 'all_mean_meth',
    'CG_n_sites', 'CH_n_sites', 'all_n_sites'
]]
stats.to_csv('results/methylation_stats.csv')
```

---

## 🎨 Advanced Usage

### Custom Genomic Regions

```python
# Load custom regions from BED file
custom_regions = pd.read_csv(
    'enhancers.bed',
    sep='\t',
    header=None,
    names=['chr', 'start', 'end', 'name']
)

# Create multi-context data
mdata_custom = create_multi_context_data_from_allc(
    allc_files='raw/allc/*.allc.tsv.gz',
    regions=custom_regions[['chr', 'start', 'end']],
    contexts=['CG', 'CH'],
    region_type='enhancers',
    n_jobs=100
)
```

### Context-Specific Analysis

```python
# Find markers for CG context
mdata_tss.set_active_context('CG')
cg_markers = ms2.tl.find_markers(
    mdata_tss,
    groupby='leiden_CG',
    method='wilcoxon'
)

# Find markers for CH context
mdata_tss.set_active_context('CH')
ch_markers = ms2.tl.find_markers(
    mdata_tss,
    groupby='leiden_CH',
    method='wilcoxon'
)

# Compare marker regions
common_markers = set(cg_markers['region']) & set(ch_markers['region'])
print(f"Common markers: {len(common_markers)}")
```

### Integration with scRNA-seq

```python
# Load scRNA-seq data
import scanpy as sc
rna_adata = sc.read_h5ad('scrna_data.h5ad')

# Match cells
common_cells = set(mdata_tss.obs_names) & set(rna_adata.obs_names)
mdata_sub = mdata_tss[list(common_cells), :].copy()
rna_sub = rna_adata[list(common_cells), :].copy()

# Compare clustering
from sklearn.metrics import adjusted_rand_score

mdata_sub.set_active_context('CG')
ari_rna_cg = adjusted_rand_score(
    mdata_sub.obs['leiden_CG'],
    rna_sub.obs['leiden']
)

mdata_sub.set_active_context('CH')
ari_rna_ch = adjusted_rand_score(
    mdata_sub.obs['leiden_CH'],
    rna_sub.obs['leiden']
)

print(f"ARI RNA vs CG: {ari_rna_cg:.3f}")
print(f"ARI RNA vs CH: {ari_rna_ch:.3f}")
```

---

## 📊 Output Structure

```
results/
├── multi_context_tss.h5ad              # Main object (TSS regions)
├── mdata_gene_body.h5ad                # Gene body regions
├── mdata_upstream.h5ad                 # Upstream regions
├── mdata_downstream.h5ad               # Downstream regions
├── clusters_all_contexts.csv           # Cluster assignments
├── methylation_stats.csv               # Methylation statistics
├── context_correlations.csv            # Context correlations
├── umap_CG.pdf                         # CG UMAP
├── umap_CH.pdf                         # CH UMAP
├── umap_all.pdf                        # All mC UMAP
└── umap_all_contexts.pdf               # Multi-panel comparison
```

---

## 💡 Best Practices

### 1. Context Selection

- **CG**: Most informative, stable, use for primary analysis
- **CH**: Non-CG contexts, useful for identifying specific cell types
- **CHG/CHH**: More specific CH subtypes
- **all**: Overall methylation, useful for global patterns

### 2. Region Selection

- **TSS**: Gene regulation, cell type identity
- **Gene body**: Gene expression correlation
- **Promoters**: Similar to TSS, standard definition
- **Enhancers**: Distal regulation (requires custom BED)

### 3. QC Thresholds

```python
# Typical filters for human brain
CG:  n_sites >= 50000, methylation 0.3-0.8
CH:  n_sites >= 10000, methylation 0.01-0.05
all: n_sites >= 60000, methylation 0.25-0.75
```

### 4. Computational Efficiency

```python
# For large datasets (1000+ cells):
# 1. Use n_jobs=100 for parallelization
# 2. Process regions separately
# 3. Filter regions before quantification
# 4. Use genome-wide bins instead of all genes
```

---

## 🔬 Use Cases

### 1. Cell Type Identification

```python
# Compare CG and CH for cell type discrimination
# CG: broad cell types
# CH: fine-grained subtypes
```

### 2. Developmental Studies

```python
# Track methylation changes across pseudotime
# Different contexts may show different dynamics
```

### 3. Disease Studies

```python
# Compare contexts between conditions
# Some diseases may affect specific contexts
```

---

## 📞 Support

For questions:
1. Check `examples/multi_context_analysis.py`
2. See API documentation
3. Open GitHub issue

---

## ✅ Summary

MethSCAn2's multi-context analysis enables:

✅ **Unified storage** of all methylation contexts  
✅ **Flexible switching** between contexts  
✅ **Parallel analysis** of multiple contexts  
✅ **Direct comparison** within same cells  
✅ **Genomic region** quantification  
✅ **Full integration** with standard analysis  

**One object, all contexts, complete analysis!** 🎉
