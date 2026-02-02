# Multi-Context Data H5AD Storage Guide

## ✅ Yes! MultiContextMethylationData Fully Supports H5AD Format

MultiContextMethylationData继承自AnnData,完全支持.h5ad文件格式存储和读取。

---

## 🔄 基本存储和读取

### 保存到H5AD

```python
import methscan2 as ms2

# 创建多上下文数据
mdata = ms2.create_multi_context_data_from_allc(
    allc_files='raw/allc/*.allc.tsv.gz',
    regions=regions_df,
    contexts=['CG', 'CH', 'CHG', 'CHH', 'all'],
    region_type='tss',
    n_jobs=100
)

# 运行分析
mdata.set_active_context('CG')
ms2.tl.run_pca(mdata)
ms2.tl.run_umap(mdata)
ms2.tl.run_leiden(mdata, key_added='leiden_CG')

# 保存为H5AD - 所有数据和分析结果都会保存
mdata.write('my_methylation_data.h5ad')

# 也支持压缩
mdata.write('my_methylation_data.h5ad', compression='gzip')
```

### 从H5AD读取

```python
import methscan2 as ms2

# 方法1: 使用MethSCAn2的读取函数 (推荐)
from methscan2.core.multi_context_data import MultiContextMethylationData

mdata = MultiContextMethylationData.read_h5ad('my_methylation_data.h5ad')

# 方法2: 使用AnnData读取,然后转换
import anndata
adata = anndata.read_h5ad('my_methylation_data.h5ad')
mdata = MultiContextMethylationData(
    X=adata.X,
    obs=adata.obs,
    var=adata.var,
    uns=adata.uns,
    obsm=adata.obsm,
    varm=adata.varm,
    layers=adata.layers,
    obsp=adata.obsp,
    varp=adata.varp
)

# 检查加载的数据
print(mdata)
print(f"Available contexts: {mdata.available_contexts}")
print(f"Active context: {mdata.uns.get('active_context', 'Not set')}")

# 继续使用
mdata.set_active_context('CH')
ms2.pl.umap(mdata, color='leiden_CH')
```

---

## 📦 H5AD文件包含的内容

当你保存MultiContextMethylationData到h5ad时,会保存:

### 1. 主数据矩阵 (X)
```python
mdata.X  # 当前活动上下文的甲基化率
```

### 2. 所有上下文的Layers
```python
mdata.layers = {
    'CG_rate': ...,    # CG甲基化率
    'CG_met': ...,     # CG甲基化计数
    'CG_total': ...,   # CG总覆盖度
    'CH_rate': ...,    # CH甲基化率
    'CH_met': ...,     # CH甲基化计数
    'CH_total': ...,   # CH总覆盖度
    'CHG_rate': ...,   # CHG甲基化率
    'CHG_met': ...,    # CHG甲基化计数
    'CHG_total': ...,  # CHG总覆盖度
    # ... 所有上下文
}
```

### 3. 细胞元数据 (obs)
```python
mdata.obs = {
    'CG_n_sites': ...,       # CG位点数
    'CG_mean_coverage': ..., # CG平均覆盖度
    'CG_mean_methylation': ..., # CG平均甲基化
    'CH_n_sites': ...,       # CH位点数
    'CH_mean_coverage': ..., # CH平均覆盖度
    'CH_mean_methylation': ..., # CH平均甲基化
    'leiden_CG': ...,        # CG聚类结果
    'leiden_CH': ...,        # CH聚类结果
    # ... 更多QC和分析结果
}
```

### 4. 区域元数据 (var)
```python
mdata.var = {
    'chr': ...,      # 染色体
    'start': ...,    # 起始位置
    'end': ...,      # 结束位置
    'gene_id': ...,  # 基因ID (如果有)
    'gene_name': ..., # 基因名 (如果有)
    # ... 其他区域注释
}
```

### 5. 降维结果 (obsm)
```python
mdata.obsm = {
    'X_pca': ...,      # 当前PCA结果
    'X_umap': ...,     # 当前UMAP结果
    'X_pca_CG': ...,   # CG的PCA结果
    'X_umap_CG': ...,  # CG的UMAP结果
    'X_pca_CH': ...,   # CH的PCA结果
    'X_umap_CH': ...,  # CH的UMAP结果
    # ... 所有上下文的降维结果
}
```

### 6. 元信息 (uns)
```python
mdata.uns = {
    'contexts': ['CG', 'CH', 'CHG', 'CHH', 'all'],  # 可用上下文
    'active_context': 'CG',                          # 当前活动上下文
    'region_type': 'tss',                            # 区域类型
    'genome': 'hg38',                                # 基因组版本
    'pca': {...},                                    # PCA参数
    'umap': {...},                                   # UMAP参数
    # ... 其他分析参数和结果
}
```

---

## 💾 文件大小和压缩

### 文件大小估算

```python
# 对于典型的人类脑细胞数据:
# 1000 cells × 20000 TSS regions × 3 contexts × 3 layers
# 未压缩: ~2-3 GB
# 压缩 (gzip): ~500 MB - 1 GB

# 保存时使用压缩
mdata.write('data.h5ad', compression='gzip', compression_opts=9)

# 查看文件大小
import os
size_mb = os.path.getsize('data.h5ad') / (1024**2)
print(f"File size: {size_mb:.2f} MB")
```

### 压缩选项

```python
# 1. 不压缩 (最快,最大)
mdata.write('data.h5ad')

# 2. GZIP压缩 (推荐,平衡)
mdata.write('data.h5ad', compression='gzip')

# 3. GZIP最大压缩 (最小,较慢)
mdata.write('data.h5ad', compression='gzip', compression_opts=9)

# 4. LZF压缩 (快速,中等压缩)
mdata.write('data.h5ad', compression='lzf')
```

---

## 🔍 Backed Mode (内存高效读取)

对于大文件,可以使用backed mode只读取需要的部分:

```python
# Backed mode读取 - 不加载全部数据到内存
mdata = MultiContextMethylationData.read_h5ad(
    'large_data.h5ad',
    backed='r'  # 只读模式
)

# 只读取需要的部分
subset = mdata[:100, :1000]  # 前100个细胞,前1000个区域
subset_memory = subset.to_memory()  # 加载到内存

# 进行分析
ms2.tl.run_pca(subset_memory)
```

---

## 📊 完整示例:保存和读取工作流

### 第一次分析:创建并保存

```python
import methscan2 as ms2
import numpy as np

# 1. 创建多上下文数据
print("Creating multi-context data...")
mdata = ms2.create_multi_context_data_from_allc(
    allc_files='raw/allc/*.allc.tsv.gz',
    regions=tss_regions_df,
    contexts=['CG', 'CH', 'all'],
    region_type='tss',
    n_jobs=100
)

# 2. QC和过滤
print("Quality control...")
for ctx in mdata.available_contexts:
    rate = mdata.get_context(ctx, 'rate')
    total = mdata.get_context(ctx, 'total')
    mdata.obs[f'{ctx}_n_sites'] = (total > 0).sum(axis=1)
    mdata.obs[f'{ctx}_mean_meth'] = np.nanmean(rate, axis=1)

mdata = mdata[mdata.obs['CG_n_sites'] >= 1000, :].copy()

# 3. 对每个上下文分析
print("Analyzing contexts...")
for ctx in ['CG', 'CH', 'all']:
    print(f"  Processing {ctx}...")
    mdata.set_active_context(ctx)
    ms2.tl.run_pca(mdata, n_comps=50)
    ms2.tl.run_umap(mdata)
    ms2.tl.run_leiden(mdata, key_added=f'leiden_{ctx}')
    
    # 保存结果到专用的obsm键
    mdata.obsm[f'X_pca_{ctx}'] = mdata.obsm['X_pca'].copy()
    mdata.obsm[f'X_umap_{ctx}'] = mdata.obsm['X_umap'].copy()

# 4. 保存完整对象
print("Saving to H5AD...")
mdata.write('results/multi_context_tss.h5ad', compression='gzip')

print(f"Saved! File size: {os.path.getsize('results/multi_context_tss.h5ad')/1e6:.2f} MB")
```

### 后续分析:读取并继续

```python
import methscan2 as ms2
from methscan2.core.multi_context_data import MultiContextMethylationData

# 1. 读取之前保存的数据
print("Loading data from H5AD...")
mdata = MultiContextMethylationData.read_h5ad('results/multi_context_tss.h5ad')

print(mdata)
print(f"Available contexts: {mdata.available_contexts}")

# 2. 所有之前的分析结果都在!
print("\nAvailable analyses:")
print(f"  obs columns: {list(mdata.obs.columns)}")
print(f"  obsm keys: {list(mdata.obsm.keys())}")
print(f"  layers: {list(mdata.layers.keys())}")

# 3. 直接使用之前的结果
# 可视化CG上下文
ms2.pl.umap(
    mdata,
    color='leiden_CG',
    use_rep='X_umap_CG',
    save='umap_cg_reloaded.pdf'
)

# 可视化CH上下文
ms2.pl.umap(
    mdata,
    color='leiden_CH',
    use_rep='X_umap_CH',
    save='umap_ch_reloaded.pdf'
)

# 4. 继续新的分析
# 比如找marker区域
mdata.set_active_context('CG')
# markers = ms2.tl.find_markers(mdata, groupby='leiden_CG')

# 5. 添加新分析后再次保存
mdata.write('results/multi_context_tss_updated.h5ad', compression='gzip')
```

---

## 🔄 与其他工具的兼容性

### 1. Scanpy兼容

```python
import scanpy as sc

# MultiContextMethylationData可以直接用于scanpy
mdata.set_active_context('CG')

# 使用scanpy的功能
sc.pl.pca(mdata, color='leiden_CG')
sc.pl.umap(mdata, color='CG_mean_methylation')

# Scanpy的分析
sc.tl.rank_genes_groups(mdata, groupby='leiden_CG')
sc.pl.rank_genes_groups(mdata)
```

### 2. 与scRNA-seq整合

```python
import scanpy as sc

# 读取甲基化数据
mdata = MultiContextMethylationData.read_h5ad('methylation.h5ad')

# 读取RNA数据
rna = sc.read_h5ad('rna.h5ad')

# 匹配细胞
common_cells = list(set(mdata.obs_names) & set(rna.obs_names))
mdata_sub = mdata[common_cells, :].copy()
rna_sub = rna[common_cells, :].copy()

# 比较分析
mdata_sub.set_active_context('CG')
# ... 整合分析
```

### 3. 导出到其他格式

```python
# 导出obs (细胞元数据)
mdata.obs.to_csv('cell_metadata.csv')

# 导出特定上下文的数据
cg_rate = mdata.get_context('CG', 'rate')
pd.DataFrame(
    cg_rate,
    index=mdata.obs_names,
    columns=mdata.var_names
).to_csv('cg_methylation.csv')

# 导出聚类结果
clusters = mdata.obs[['leiden_CG', 'leiden_CH', 'leiden_all']]
clusters.to_csv('clusters.csv')
```

---

## 💡 最佳实践

### 1. 定期保存检查点

```python
# 在关键步骤后保存
mdata = ms2.create_multi_context_data_from_allc(...)
mdata.write('checkpoint_1_raw.h5ad')

# QC后保存
mdata = mdata[good_cells, :].copy()
mdata.write('checkpoint_2_filtered.h5ad')

# 分析后保存
# ... 运行分析
mdata.write('checkpoint_3_analyzed.h5ad')
```

### 2. 使用描述性文件名

```python
# 包含关键信息
mdata.write('brain_tss_CG-CH-all_1000cells_20241202.h5ad')

# 或使用目录结构
mdata.write('results/region_tss/contexts_CG_CH_all/data.h5ad')
```

### 3. 保存元数据

```python
# 在uns中添加分析记录
mdata.uns['analysis_date'] = '2024-12-02'
mdata.uns['analyst'] = 'YourName'
mdata.uns['description'] = 'Human brain PFC, TSS regions, CG/CH contexts'
mdata.uns['parameters'] = {
    'min_sites': 1000,
    'resolution': 0.8,
    'n_pcs': 50
}

mdata.write('data_with_metadata.h5ad')
```

---

## ⚠️ 注意事项

### 1. 版本兼容性

```python
# 保存时记录版本
import methscan2 as ms2
mdata.uns['methscan2_version'] = ms2.__version__

# 读取时检查版本
loaded = MultiContextMethylationData.read_h5ad('data.h5ad')
print(f"Created with MethSCAn2 v{loaded.uns.get('methscan2_version', 'unknown')}")
```

### 2. 大文件处理

```python
# 对于非常大的数据集 (>10GB),考虑:

# 选项1: 分别保存不同区域类型
mdata_tss.write('tss.h5ad')
mdata_gb.write('gene_body.h5ad')

# 选项2: 只保存需要的上下文
mdata_cg_only = mdata.copy()
# 删除不需要的layers
for layer in list(mdata_cg_only.layers.keys()):
    if not layer.startswith('CG_'):
        del mdata_cg_only.layers[layer]
mdata_cg_only.write('cg_only.h5ad')

# 选项3: 使用更强的压缩
mdata.write('data.h5ad', compression='gzip', compression_opts=9)
```

---

## ✅ 总结

MultiContextMethylationData **完全支持** H5AD格式:

✅ **完整保存** - 所有上下文、分析结果都会保存  
✅ **快速读取** - 完整恢复所有数据和分析  
✅ **压缩支持** - 大幅减小文件大小  
✅ **Backed模式** - 内存高效处理大文件  
✅ **完全兼容** - 与AnnData/Scanpy生态系统兼容  

**一次保存,随时使用!** 🎉
