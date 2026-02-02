# ALLC Processing Guide

## 🎯 Overview

MethSCAn2 now includes **complete ALLC file processing** capabilities that replicate and improve upon your bash workflow, with **full multi-threading support**.

### What's New

✅ **ALLC → COV conversion** with context separation (CG, CH, CHG, CHH, all)  
✅ **Global methylation calculation** for each context  
✅ **Parallel processing** (up to 100+ jobs simultaneously)  
✅ **Pure Python** - no bash scripts needed  
✅ **Automatic file organization** into context-specific directories  

---

## 🚀 Quick Start

### One-Line Complete Pipeline

```python
import methscan2 as ms2

# This does everything: convert, separate contexts, calculate global levels
global_stats = ms2.pp.batch_process_allc(
    allc_dir='raw/allc',           # Your ALLC files
    output_base_dir='processed',    # Output directory
    contexts=['CG', 'CH', 'all'],  # Contexts to extract
    n_jobs=100                      # Parallel jobs (like xargs -P 100)
)
```

**Output structure:**
```
processed/
├── CG_cov/              # CG methylation COV files
├── CH_cov/              # CH methylation COV files
├── all_cov/             # All methylation COV files
└── global_level/
    ├── all_cell_global_average_mCG.tsv
    ├── all_cell_global_average_mCH.tsv
    └── all_cell_global_average_mall.tsv
```

---

## 📋 Usage Examples

### Example 1: Command Line Tool

```bash
# Simple conversion
python examples/convert_allc.py --input raw/allc --output processed --jobs 100

# Custom contexts
python examples/convert_allc.py --input raw/allc --contexts CG CHG CHH --jobs 50

# Skip global methylation calculation
python examples/convert_allc.py --input raw/allc --no-global
```

### Example 2: Python Script (Step-by-Step)

```python
import methscan2 as ms2

# Step 1: Convert ALLC to COV
output_files = ms2.pp.allc_to_cov(
    allc_files='raw/allc/*.allc.tsv.gz',
    output_dir='cov_output',
    contexts=['CG', 'CH'],
    n_jobs=100
)
# Returns: {'CG': [...], 'CH': [...]}

# Step 2: Calculate global methylation for CG
cg_stats = ms2.pp.calculate_global_methylation(
    cov_files='cov_output/*_CG.cov',
    output_file='global_CG.tsv',
    context_name='CG',
    n_jobs=100
)

# Step 3: Calculate for CH
ch_stats = ms2.pp.calculate_global_methylation(
    cov_files='cov_output/*_CH.cov',
    output_file='global_CH.tsv',
    context_name='CH',
    n_jobs=100
)

# Step 4: Merge statistics
merged = ms2.pp.merge_global_stats(
    stats_files={'CG': 'global_CG.tsv', 'CH': 'global_CH.tsv'},
    output_file='merged_stats.tsv'
)
```

### Example 3: Complete Workflow

```python
import methscan2 as ms2

# 1. Convert ALLC to COV and calculate global levels
global_stats = ms2.pp.batch_process_allc(
    allc_dir='raw/allc',
    output_base_dir='processed',
    contexts=['CG', 'CH', 'all'],
    n_jobs=100
)

# 2. Load CG data for analysis
mdata = ms2.read_cov_files('processed/CG_cov/*.cov', n_jobs=100)

# 3. Quality control
ms2.pp.calculate_qc_metrics(mdata)
ms2.pp.filter_cells(mdata, min_sites=50000, min_methylation=0.3)

# 4. Analysis
ms2.pp.smooth_methylation(mdata, n_jobs=100)
ms2.tl.detect_vmr(mdata, n_jobs=100)
ms2.tl.quantify_methylation(mdata)

# 5. Visualization
ms2.tl.run_pca(mdata)
ms2.tl.run_umap(mdata)
ms2.pl.umap(mdata, color='mean_methylation', save='umap.pdf')
```

---

## 🔧 Function Reference

### `batch_process_allc()`

Complete automated pipeline.

```python
global_stats = ms2.pp.batch_process_allc(
    allc_dir='raw/allc',              # Input directory
    output_base_dir='processed',       # Output directory
    contexts=['CG', 'CH', 'all'],     # Contexts to extract
    calculate_global=True,             # Calculate global methylation
    n_jobs=100,                        # Parallel jobs
    verbose=True                       # Show progress
)
```

**Returns:** Dictionary of pandas DataFrames with global statistics

### `allc_to_cov()`

Convert ALLC files to COV format.

```python
output_files = ms2.pp.allc_to_cov(
    allc_files='*.allc.tsv.gz',       # Input files
    output_dir='cov_files',            # Output directory
    contexts=['CG', 'CH'],             # Contexts
    separate_contexts=True,            # Separate files per context
    compress=False,                    # Gzip output files
    n_jobs=-1                          # Parallel jobs
)
```

**Returns:** Dictionary mapping contexts to output file lists

### `calculate_global_methylation()`

Calculate global methylation levels.

```python
stats = ms2.pp.calculate_global_methylation(
    cov_files='*.cov',                 # Input COV files
    output_file='global_stats.tsv',    # Output file
    context_name='CG',                 # Context name
    n_jobs=-1                          # Parallel jobs
)
```

**Returns:** pandas DataFrame with columns: cell, global_level, n_sites

### `merge_global_stats()`

Merge statistics from multiple contexts.

```python
merged = ms2.pp.merge_global_stats(
    stats_files={
        'CG': 'global_CG.tsv',
        'CH': 'global_CH.tsv'
    },
    output_file='merged.tsv'
)
```

**Returns:** Merged pandas DataFrame

---

## 📊 File Formats

### Input: ALLC Format

```
chr1  3000827  +  CCC  0  1  1
chr1  3000839  +  CCT  0  1  1
chr1  3001007  +  CGA  1  0  1
```

Columns: chr, pos, strand, context, mc_count, total_count, methylation_level

### Output: COV Format

```
chr1  3000827  3000828  0    0  1
chr1  3000839  3000840  0    0  1
chr1  3001007  3001008  100  1  0
```

Columns: chr, start, end, meth_pct, met_count, unmet_count

### Global Stats Format

```
cell             global_CG_level  CG_sites
cell_001         0.723            1234567
cell_002         0.698            1198765
```

---

## ⚡ Performance Comparison

### Your Original Bash Script:
```bash
# xargs -P 100 for parallel processing
# Multiple commands piped together
# Separate scripts for conversion and stats
```

### MethSCAn2 Python:
```python
# Single function call
# Built-in parallelization (n_jobs=100)
# Integrated pipeline
# Progress bars and logging
```

**Performance:** Similar speed (both use 100 parallel jobs)  
**Advantages:** 
- ✅ Easier to use
- ✅ Better error handling
- ✅ Progress tracking
- ✅ Integrated with analysis pipeline

---

## 🎯 Comparison with Your Workflow

| Your Bash Script | MethSCAn2 Python |
|-----------------|------------------|
| `allcools extract-allc` | `ms2.pp.extract_context_from_allc()` |
| `gunzip + awk` for COV conversion | `ms2.pp.allc_to_cov()` |
| Shell loops for global stats | `ms2.pp.calculate_global_methylation()` |
| Manual file organization | Automatic directory structure |
| `xargs -P 100` | `n_jobs=100` |
| Multiple scripts | Single pipeline |

---

## 💡 Tips & Best Practices

### 1. Memory Management

For large datasets (1000+ cells):
```python
# Process in batches
import glob
files = glob.glob('raw/allc/*.allc.tsv.gz')
batch_size = 100

for i in range(0, len(files), batch_size):
    batch = files[i:i+batch_size]
    ms2.pp.allc_to_cov(
        allc_files=batch,
        output_dir=f'output/batch_{i}',
        n_jobs=50
    )
```

### 2. Custom Contexts

Extract specific contexts:
```python
output_files = ms2.pp.allc_to_cov(
    allc_files='*.allc.tsv.gz',
    contexts=['CG', 'CHG', 'CHH'],  # Separate CH into CHG and CHH
    n_jobs=100
)
```

### 3. Filtering Based on Global Stats

```python
# Calculate global stats
stats = ms2.pp.calculate_global_methylation(
    cov_files='*.cov',
    context_name='CG'
)

# Filter cells
good_cells = stats[
    (stats['global_CG_level'] > 0.3) &
    (stats['global_CG_level'] < 0.8) &
    (stats['CG_sites'] > 50000)
]['cell'].tolist()

# Load only good cells
mdata = ms2.read_cov_files(
    [f'{cell}.cov' for cell in good_cells]
)
```

---

## 🐛 Troubleshooting

### "No files found"
```python
# Check glob pattern
import glob
files = glob.glob('raw/allc/*.allc.tsv.gz')
print(f"Found {len(files)} files")
```

### Out of memory
```python
# Reduce parallel jobs
ms2.pp.batch_process_allc(..., n_jobs=20)

# Or process in batches (see Tips section)
```

### Slow processing
```python
# Increase parallel jobs
ms2.pp.batch_process_allc(..., n_jobs=100)

# Check CPU usage
import os
print(f"Available CPUs: {os.cpu_count()}")
```

---

## 📚 Complete Example

See `examples/allc_processing_pipeline.py` for a complete, documented example that includes:
- ALLC → COV conversion
- Context separation
- Global methylation calculation
- Quality control
- Downstream analysis

Run with:
```bash
python examples/allc_processing_pipeline.py
```

---

## ✅ Summary

MethSCAn2 now provides a **complete, efficient ALLC processing pipeline** that:

✅ Replicates your bash workflow in pure Python  
✅ Supports 100+ parallel jobs (like `xargs -P 100`)  
✅ Automatically organizes outputs  
✅ Integrates with downstream analysis  
✅ Includes progress tracking and error handling  

**One command does everything:**
```python
ms2.pp.batch_process_allc('raw/allc', 'processed', n_jobs=100)
```

That's it! 🎉
