"""
MethSCAn2 Core Data Structure
基于AnnData的甲基化数据类
"""

from anndata import AnnData
import numpy as np
import pandas as pd
from typing import Optional, Union, List, Tuple
from pathlib import Path
import warnings


class MethylationData(AnnData):
    """
    单细胞甲基化数据类,扩展AnnData以适配甲基化特定需求
    
    主要扩展功能:
    1. 甲基化特定的数据层管理
    2. 基因组坐标处理
    3. VMR/DMR信息存储
    4. 甲基化率计算
    
    Attributes:
        X: 主数据矩阵 (cells × regions),默认存储甲基化率
        obs: 细胞元数据 DataFrame
        var: 区域元数据 DataFrame,必须包含染色体坐标信息
        layers: 数据层字典
            - 'met': 甲基化读数计数
            - 'total': 总读数计数  
            - 'rate': 甲基化率
            - 'residuals': 收缩残差
        obsm: 降维结果等矩阵 (cells × features)
            - 'X_pca': PCA结果
            - 'X_umap': UMAP结果
            - 'X_spectral': 谱嵌入结果
        varm: 区域特征矩阵 (regions × features)
        uns: 非结构化元数据字典
            - 'vmr': VMR检测结果
            - 'dmr': DMR检测结果
            - 'genome': 基因组信息
        obsp: 细胞间关系矩阵 (如邻接矩阵)
        varp: 区域间关系矩阵
    """
    
    def __init__(
        self,
        X: Optional[np.ndarray] = None,
        obs: Optional[pd.DataFrame] = None,
        var: Optional[pd.DataFrame] = None,
        uns: Optional[dict] = None,
        obsm: Optional[dict] = None,
        varm: Optional[dict] = None,
        layers: Optional[dict] = None,
        obsp: Optional[dict] = None,
        varp: Optional[dict] = None,
        dtype: str = 'float32',
        shape: Optional[Tuple[int, int]] = None,
        filename: Optional[Path] = None,
        filemode: Optional[str] = None,
    ):
        """
        初始化MethylationData对象
        
        Parameters:
            X: 主数据矩阵,通常是甲基化率
            obs: 细胞元数据
            var: 区域元数据,应包含chrom, start, end列
            其他参数同AnnData
        """
        super().__init__(
            X=X, obs=obs, var=var, uns=uns, obsm=obsm, varm=varm,
            layers=layers, obsp=obsp, varp=varp, dtype=dtype,
            shape=shape, filename=filename, filemode=filemode
        )
        
        # 验证var是否包含必要的基因组坐标
        if var is not None and len(var) > 0:
            self._validate_genomic_coords()
        
        # 初始化甲基化特定的uns字段
        if 'genome' not in self.uns:
            self.uns['genome'] = {}
        if 'vmr' not in self.uns:
            self.uns['vmr'] = None
        if 'dmr' not in self.uns:
            self.uns['dmr'] = {}
    
    def _validate_genomic_coords(self):
        """验证var是否包含基因组坐标信息"""
        required_cols = ['chrom', 'start', 'end']
        missing = [col for col in required_cols if col not in self.var.columns]
        
        if missing:
            warnings.warn(
                f"var缺少基因组坐标列: {missing}. "
                "某些功能可能无法使用。"
            )
    
    @property
    def methylation_rate(self) -> np.ndarray:
        """
        计算或返回甲基化率矩阵
        
        Returns:
            甲基化率矩阵 (cells × regions)
        """
        # 如果已有rate层,直接返回
        if 'rate' in self.layers:
            return self.layers['rate']
        
        # 如果有met和total层,计算
        if 'met' in self.layers and 'total' in self.layers:
            with np.errstate(divide='ignore', invalid='ignore'):
                rate = self.layers['met'].astype(float) / self.layers['total']
                rate[~np.isfinite(rate)] = np.nan
            return rate
        
        # 否则假设X就是甲基化率
        return self.X
    
    @property
    def n_regions(self) -> int:
        """返回区域数量"""
        return self.n_vars
    
    @property
    def chromosomes(self) -> List[str]:
        """返回包含的染色体列表"""
        if 'chrom' not in self.var.columns:
            return []
        return sorted(self.var['chrom'].unique())
    
    def add_vmr(
        self,
        vmr_df: pd.DataFrame,
        key: str = 'is_vmr',
        overlap_threshold: float = 0.5
    ):
        """
        标记VMR区域
        
        Parameters:
            vmr_df: VMR DataFrame,应包含chrom, start, end列
            key: 添加到var的列名
            overlap_threshold: 重叠比例阈值
        """
        if 'chrom' not in self.var.columns:
            raise ValueError("var中缺少基因组坐标信息")
        
        # 初始化标记列
        self.var[key] = False
        
        # 对每个VMR,找到重叠的区域
        for _, vmr in vmr_df.iterrows():
            # 找到同染色体的区域
            chrom_mask = self.var['chrom'] == vmr['chrom']
            
            # 计算重叠
            overlap_mask = (
                (self.var['end'] > vmr['start']) &
                (self.var['start'] < vmr['end'])
            )
            
            # 计算重叠比例
            overlap_start = np.maximum(self.var['start'], vmr['start'])
            overlap_end = np.minimum(self.var['end'], vmr['end'])
            overlap_len = np.maximum(0, overlap_end - overlap_start)
            region_len = self.var['end'] - self.var['start']
            
            overlap_ratio = overlap_len / region_len
            
            # 标记满足阈值的区域
            mark_mask = chrom_mask & overlap_mask & (overlap_ratio >= overlap_threshold)
            self.var.loc[mark_mask, key] = True
        
        # 保存VMR信息到uns
        self.uns['vmr'] = vmr_df
        
        print(f"标记了 {self.var[key].sum()} 个VMR区域")
    
    def subset_by_region(
        self,
        chrom: Optional[str] = None,
        start: Optional[int] = None,
        end: Optional[int] = None,
        region_str: Optional[str] = None
    ) -> 'MethylationData':
        """
        按基因组区域切片数据
        
        Parameters:
            chrom: 染色体名称
            start: 起始位置
            end: 结束位置
            region_str: 区域字符串,格式为"chr:start-end"
        
        Returns:
            切片后的MethylationData对象
        """
        if region_str is not None:
            # 解析区域字符串
            chrom, coords = region_str.split(':')
            start, end = map(int, coords.split('-'))
        
        if 'chrom' not in self.var.columns:
            raise ValueError("var中缺少基因组坐标信息")
        
        # 构建过滤mask
        mask = self.var['chrom'] == chrom
        
        if start is not None:
            mask &= self.var['end'] > start
        if end is not None:
            mask &= self.var['start'] < end
        
        return self[:, mask].copy()
    
    def subset_by_cells(
        self,
        cell_ids: Optional[List[str]] = None,
        cell_mask: Optional[np.ndarray] = None,
        obs_key: Optional[str] = None,
        obs_values: Optional[List] = None
    ) -> 'MethylationData':
        """
        按细胞切片数据
        
        Parameters:
            cell_ids: 细胞ID列表
            cell_mask: 布尔mask数组
            obs_key: obs中的列名
            obs_values: 该列的值列表
        
        Returns:
            切片后的MethylationData对象
        """
        if cell_ids is not None:
            mask = self.obs_names.isin(cell_ids)
        elif cell_mask is not None:
            mask = cell_mask
        elif obs_key is not None and obs_values is not None:
            mask = self.obs[obs_key].isin(obs_values)
        else:
            raise ValueError("必须提供cell_ids, cell_mask或obs_key+obs_values")
        
        return self[mask, :].copy()
    
    def calculate_coverage_stats(self) -> pd.DataFrame:
        """
        计算覆盖度统计
        
        Returns:
            包含覆盖度统计的DataFrame
        """
        if 'total' not in self.layers:
            warnings.warn("缺少total层,无法计算覆盖度")
            return pd.DataFrame()
        
        total = self.layers['total']
        
        stats = pd.DataFrame({
            'total_sites': (total > 0).sum(axis=1),
            'mean_coverage': np.nanmean(np.where(total > 0, total, np.nan), axis=1),
            'median_coverage': np.nanmedian(np.where(total > 0, total, np.nan), axis=1),
        }, index=self.obs_names)
        
        return stats
    
    def calculate_methylation_stats(self) -> pd.DataFrame:
        """
        计算甲基化统计
        
        Returns:
            包含甲基化统计的DataFrame
        """
        rate = self.methylation_rate
        
        stats = pd.DataFrame({
            'mean_methylation': np.nanmean(rate, axis=1),
            'median_methylation': np.nanmedian(rate, axis=1),
            'methylation_variance': np.nanvar(rate, axis=1),
        }, index=self.obs_names)
        
        return stats
    
    def get_region_methylation(
        self,
        region_idx: Union[int, str],
        return_coverage: bool = False
    ) -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray]]:
        """
        获取特定区域的甲基化值
        
        Parameters:
            region_idx: 区域索引或名称
            return_coverage: 是否同时返回覆盖度
        
        Returns:
            甲基化率数组,如果return_coverage=True,还返回覆盖度数组
        """
        if isinstance(region_idx, str):
            idx = self.var_names.get_loc(region_idx)
        else:
            idx = region_idx
        
        meth_rate = self.methylation_rate[:, idx]
        
        if return_coverage and 'total' in self.layers:
            coverage = self.layers['total'][:, idx]
            return meth_rate, coverage
        
        return meth_rate
    
    def __repr__(self) -> str:
        """自定义字符串表示"""
        descr = f"MethylationData object with n_obs × n_vars = {self.n_obs} × {self.n_vars}"
        
        if 'chrom' in self.var.columns:
            n_chroms = len(self.chromosomes)
            descr += f"\n    {n_chroms} chromosomes"
        
        if 'vmr' in self.uns and self.uns['vmr'] is not None:
            n_vmrs = len(self.uns['vmr'])
            descr += f"\n    {n_vmrs} VMRs detected"
        
        # 显示可用的层
        if self.layers:
            descr += f"\n    layers: {', '.join(self.layers.keys())}"
        
        # 显示降维结果
        if self.obsm:
            descr += f"\n    obsm: {', '.join(self.obsm.keys())}"
        
        return descr


# 便利函数
def concat(
    mdatas: List[MethylationData],
    axis: int = 0,
    join: str = 'inner',
    label: Optional[str] = None,
    keys: Optional[List[str]] = None,
    index_unique: Optional[str] = None
) -> MethylationData:
    """
    合并多个MethylationData对象
    
    Parameters:
        mdatas: MethylationData对象列表
        axis: 合并轴,0为按细胞合并,1为按区域合并
        join: 合并方式,'inner'或'outer'
        label: 添加到obs/var的批次标签列名
        keys: 批次标签值列表
        index_unique: 如何处理重复索引
    
    Returns:
        合并后的MethylationData对象
    """
    from anndata import concat as anndata_concat
    
    # 使用AnnData的concat函数
    result = anndata_concat(
        mdatas,
        axis=axis,
        join=join,
        label=label,
        keys=keys,
        index_unique=index_unique
    )
    
    # 转换为MethylationData
    mdata = MethylationData(
        X=result.X,
        obs=result.obs,
        var=result.var,
        uns=result.uns,
        obsm=result.obsm if hasattr(result, 'obsm') else None,
        varm=result.varm if hasattr(result, 'varm') else None,
        layers=result.layers if hasattr(result, 'layers') else None,
        obsp=result.obsp if hasattr(result, 'obsp') else None,
        varp=result.varp if hasattr(result, 'varp') else None,
    )
    
    return mdata
