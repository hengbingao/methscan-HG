"""
MethSCAn2 Embedding Plotting Functions
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Optional, Union, List, Tuple
import warnings


def umap(
    mdata,
    color: Union[str, List[str], None] = None,
    use_rep: str = 'X_umap',
    size: Optional[float] = None,
    alpha: float = 0.8,
    palette: Optional[Union[str, dict]] = None,
    legend_loc: str = 'right margin',
    legend_fontsize: Optional[int] = None,
    figsize: Optional[Tuple[float, float]] = None,
    ncols: int = 4,
    wspace: float = 0.3,
    title: Optional[Union[str, List[str]]] = None,
    save: Optional[str] = None,
    dpi: int = 300,
    **kwargs
):
    """
    UMAP可视化
    
    Parameters:
        mdata: MethylationData对象
        color: 颜色变量,可以是obs中的列名或列名列表
        use_rep: obsm中使用的坐标键
        size: 点的大小
        alpha: 透明度
        palette: 颜色方案
        legend_loc: 图例位置,'right margin', 'on data', 或None
        legend_fontsize: 图例字体大小
        figsize: 图片大小
        ncols: 多面板时的列数
        wspace: 子图间距
        title: 标题
        save: 保存路径
        dpi: 分辨率
        **kwargs: 传递给scatter的其他参数
    
    Returns:
        matplotlib Figure对象
    """
    # 检查坐标是否存在
    if use_rep not in mdata.obsm:
        raise ValueError(f"{use_rep} not found in mdata.obsm")
    
    coords = mdata.obsm[use_rep]
    
    # 处理color参数
    if color is None:
        color = ['grey']
    elif isinstance(color, str):
        color = [color]
    
    # 计算布局
    n_plots = len(color)
    if figsize is None:
        if n_plots == 1:
            figsize = (6, 5)
        else:
            nrows = int(np.ceil(n_plots / ncols))
            figsize = (ncols * 5, nrows * 4)
    
    # 创建子图
    if n_plots == 1:
        fig, ax = plt.subplots(figsize=figsize)
        axes = [ax]
    else:
        nrows = int(np.ceil(n_plots / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
        axes = axes.flatten() if n_plots > 1 else [axes]
    
    # 绘制每个面板
    for i, (ax, c) in enumerate(zip(axes[:n_plots], color)):
        # 确定颜色值
        if c == 'grey' or c not in mdata.obs.columns:
            # 单色
            scatter = ax.scatter(
                coords[:, 0], coords[:, 1],
                c='grey', s=size if size else 10,
                alpha=alpha, rasterized=True, **kwargs
            )
        else:
            # 从obs获取值
            values = mdata.obs[c]
            
            if pd.api.types.is_categorical_dtype(values) or values.dtype == 'object':
                # 分类变量
                unique_vals = values.unique()
                
                # 生成颜色映射
                if palette is None:
                    if len(unique_vals) <= 20:
                        palette = sns.color_palette('tab20', len(unique_vals))
                    else:
                        palette = sns.color_palette('husl', len(unique_vals))
                elif isinstance(palette, str):
                    palette = sns.color_palette(palette, len(unique_vals))
                
                # 创建颜色字典
                if isinstance(palette, list):
                    color_map = dict(zip(unique_vals, palette))
                else:
                    color_map = palette
                
                # 为每个类别绘制
                for val in unique_vals:
                    mask = values == val
                    ax.scatter(
                        coords[mask, 0], coords[mask, 1],
                        c=[color_map[val]], label=val,
                        s=size if size else 10,
                        alpha=alpha, rasterized=True, **kwargs
                    )
                
                # 添加图例
                if legend_loc == 'right margin':
                    ax.legend(
                        bbox_to_anchor=(1.05, 1), loc='upper left',
                        fontsize=legend_fontsize, frameon=False
                    )
                elif legend_loc == 'on data':
                    ax.legend(
                        loc='best', fontsize=legend_fontsize, frameon=False
                    )
            
            else:
                # 连续变量
                scatter = ax.scatter(
                    coords[:, 0], coords[:, 1],
                    c=values, cmap=palette if palette else 'viridis',
                    s=size if size else 10,
                    alpha=alpha, rasterized=True, **kwargs
                )
                
                # 添加colorbar
                cbar = plt.colorbar(scatter, ax=ax)
                cbar.set_label(c, fontsize=legend_fontsize)
        
        # 设置标签和标题
        ax.set_xlabel('UMAP1')
        ax.set_ylabel('UMAP2')
        
        if title is not None:
            if isinstance(title, list):
                ax.set_title(title[i])
            else:
                ax.set_title(title if n_plots == 1 else f'{title} - {c}')
        else:
            ax.set_title(c)
        
        # 移除坐标轴刻度
        ax.set_xticks([])
        ax.set_yticks([])
        
        # 设置等比例
        ax.set_aspect('equal', 'box')
    
    # 隐藏多余的子图
    for ax in axes[n_plots:]:
        ax.set_visible(False)
    
    plt.tight_layout()
    fig.subplots_adjust(wspace=wspace)
    
    if save:
        plt.savefig(save, dpi=dpi, bbox_inches='tight')
        print(f"Figure saved to: {save}")
    
    return fig


def heatmap(
    mdata,
    var_names: Optional[List[str]] = None,
    groupby: Optional[str] = None,
    n_vars: int = 50,
    use_raw: bool = False,
    standard_scale: Optional[str] = 'var',
    cmap: str = 'RdBu_r',
    dendrogram: bool = True,
    figsize: Optional[Tuple[float, float]] = None,
    save: Optional[str] = None,
    **kwargs
):
    """
    甲基化热图
    
    Parameters:
        mdata: MethylationData对象
        var_names: 要显示的区域名称列表
        groupby: 按此obs列分组
        n_vars: 如果var_names=None,显示top n个高变异区域
        use_raw: 是否使用原始数据
        standard_scale: 标准化方式,'var'(按列), 'obs'(按行), None
        cmap: 颜色映射
        dendrogram: 是否显示树状图
        figsize: 图片大小
        save: 保存路径
        **kwargs: 传递给seaborn.clustermap的其他参数
    
    Returns:
        seaborn ClusterGrid对象
    """
    # 获取数据矩阵
    if use_raw and hasattr(mdata, 'raw'):
        data = mdata.raw.X
        var = mdata.raw.var
    else:
        data = mdata.methylation_rate
        var = mdata.var
    
    # 选择变量
    if var_names is not None:
        var_idx = [var.index.get_loc(v) for v in var_names if v in var.index]
    else:
        # 选择高变异区域
        variances = np.nanvar(data, axis=0)
        var_idx = np.argsort(variances)[-n_vars:]
    
    # 提取子集
    plot_data = data[:, var_idx]
    plot_var_names = var.index[var_idx]
    
    # 转换为DataFrame
    df = pd.DataFrame(
        plot_data,
        index=mdata.obs_names,
        columns=plot_var_names
    )
    
    # 标准化
    if standard_scale == 'var':
        df = (df - df.mean()) / df.std()
    elif standard_scale == 'obs':
        df = df.sub(df.mean(axis=1), axis=0).div(df.std(axis=1), axis=0)
    
    # 处理分组
    if groupby is not None and groupby in mdata.obs:
        # 按组排序
        df = df.loc[mdata.obs.sort_values(groupby).index]
        
        # 创建行颜色
        row_colors = mdata.obs.loc[df.index, groupby]
        
        # 为分类变量创建颜色映射
        if pd.api.types.is_categorical_dtype(row_colors):
            unique_groups = row_colors.cat.categories
            palette = sns.color_palette('tab20', len(unique_groups))
            lut = dict(zip(unique_groups, palette))
            row_colors = row_colors.map(lut)
    else:
        row_colors = None
    
    # 设置图片大小
    if figsize is None:
        figsize = (10, 8)
    
    # 绘制热图
    g = sns.clustermap(
        df.T,  # 转置:区域为行,细胞为列
        cmap=cmap,
        center=0 if standard_scale else None,
        row_cluster=dendrogram,
        col_cluster=dendrogram,
        col_colors=row_colors,
        figsize=figsize,
        cbar_kws={'label': 'Methylation level (scaled)'},
        **kwargs
    )
    
    # 调整标签
    g.ax_heatmap.set_xlabel('Cells')
    g.ax_heatmap.set_ylabel('Regions')
    
    if save:
        plt.savefig(save, dpi=300, bbox_inches='tight')
        print(f"Figure saved to: {save}")
    
    return g


def qc_violin(
    mdata,
    keys: List[str] = ['n_sites', 'mean_coverage', 'mean_methylation'],
    groupby: Optional[str] = None,
    figsize: Optional[Tuple[float, float]] = None,
    save: Optional[str] = None,
    **kwargs
):
    """
    QC指标小提琴图
    
    Parameters:
        mdata: MethylationData对象
        keys: 要绘制的obs列名列表
        groupby: 分组变量
        figsize: 图片大小
        save: 保存路径
        **kwargs: 传递给seaborn.violinplot的其他参数
    
    Returns:
        matplotlib Figure对象
    """
    # 准备数据
    plot_keys = [k for k in keys if k in mdata.obs.columns]
    
    if len(plot_keys) == 0:
        raise ValueError(f"None of {keys} found in mdata.obs")
    
    # 创建子图
    if figsize is None:
        figsize = (len(plot_keys) * 4, 4)
    
    fig, axes = plt.subplots(1, len(plot_keys), figsize=figsize)
    if len(plot_keys) == 1:
        axes = [axes]
    
    # 绘制每个指标
    for ax, key in zip(axes, plot_keys):
        if groupby is not None and groupby in mdata.obs:
            sns.violinplot(
                data=mdata.obs,
                x=groupby,
                y=key,
                ax=ax,
                **kwargs
            )
            ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
        else:
            sns.violinplot(
                data=mdata.obs,
                y=key,
                ax=ax,
                **kwargs
            )
        
        ax.set_ylabel(key)
        ax.set_xlabel('')
    
    plt.tight_layout()
    
    if save:
        plt.savefig(save, dpi=300, bbox_inches='tight')
        print(f"Figure saved to: {save}")
    
    return fig


def genome_track(
    mdata,
    region: str,
    cells: Optional[List[str]] = None,
    group_by: Optional[str] = None,
    show_vmr: bool = True,
    figsize: Optional[Tuple[float, float]] = None,
    save: Optional[str] = None
):
    """
    基因组区域的甲基化轨迹图
    
    Parameters:
        mdata: MethylationData对象
        region: 区域字符串,格式"chr:start-end"
        cells: 要显示的细胞列表
        group_by: 按此列分组显示
        show_vmr: 是否标记VMR
        figsize: 图片大小
        save: 保存路径
    
    Returns:
        matplotlib Figure对象
    """
    # 解析区域
    chrom, coords = region.split(':')
    start, end = map(int, coords.split('-'))
    
    # 提取区域数据
    region_data = mdata.subset_by_region(chrom=chrom, start=start, end=end)
    
    if region_data.n_vars == 0:
        raise ValueError(f"No data found in region {region}")
    
    # 选择细胞
    if cells is not None:
        region_data = region_data[cells, :]
    
    # 获取位置和甲基化数据
    positions = region_data.var['start'].values
    meth_data = region_data.methylation_rate
    
    # 创建图形
    if figsize is None:
        figsize = (12, 6)
    
    fig, ax = plt.subplots(figsize=figsize)
    
    # 绘制每个细胞的轨迹
    if group_by is not None and group_by in region_data.obs:
        groups = region_data.obs[group_by].unique()
        palette = sns.color_palette('tab10', len(groups))
        
        for i, group in enumerate(groups):
            group_mask = region_data.obs[group_by] == group
            group_data = meth_data[group_mask, :]
            
            # 绘制组的平均值
            mean_meth = np.nanmean(group_data, axis=0)
            ax.plot(positions, mean_meth, label=group, color=palette[i], linewidth=2)
            
            # 添加置信区间
            std_meth = np.nanstd(group_data, axis=0)
            ax.fill_between(
                positions,
                mean_meth - std_meth,
                mean_meth + std_meth,
                alpha=0.2,
                color=palette[i]
            )
    else:
        # 绘制所有细胞的平均值
        mean_meth = np.nanmean(meth_data, axis=0)
        ax.plot(positions, mean_meth, color='blue', linewidth=2)
    
    # 标记VMR
    if show_vmr and 'vmr' in mdata.uns and mdata.uns['vmr'] is not None:
        vmr_df = mdata.uns['vmr']
        region_vmrs = vmr_df[
            (vmr_df['chrom'] == chrom) &
            (vmr_df['start'] < end) &
            (vmr_df['end'] > start)
        ]
        
        for _, vmr in region_vmrs.iterrows():
            ax.axvspan(
                max(vmr['start'], start),
                min(vmr['end'], end),
                alpha=0.2, color='red', label='VMR'
            )
    
    ax.set_xlabel(f'Position on {chrom} (bp)')
    ax.set_ylabel('Methylation level')
    ax.set_title(f'Methylation track: {region}')
    ax.set_xlim(start, end)
    ax.set_ylim(0, 1)
    
    if group_by is not None or show_vmr:
        ax.legend()
    
    plt.tight_layout()
    
    if save:
        plt.savefig(save, dpi=300, bbox_inches='tight')
        print(f"Figure saved to: {save}")
    
    return fig
