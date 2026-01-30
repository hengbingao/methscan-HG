"""
MethSCAn2 VMR Detection Module
变异甲基化区域(Variably Methylated Regions)检测
"""

import numpy as np
import pandas as pd
from scipy.ndimage import uniform_filter1d
from scipy.stats import median_abs_deviation
from joblib import Parallel, delayed
from typing import Optional, Union, Tuple
import warnings
from tqdm.auto import tqdm


def detect_vmr(
    mdata,
    bandwidth: int = 2000,
    stepsize: int = 100,
    var_threshold: float = 0.02,
    min_coverage: int = 5,
    min_cells: int = 10,
    chromosomes: Optional[list] = None,
    n_jobs: int = -1,
    verbose: bool = True,
    inplace: bool = True
) -> Optional[pd.DataFrame]:
    """
    检测变异甲基化区域(VMR)
    
    算法步骤:
    1. 计算平滑的全局甲基化曲线(伪bulk)
    2. 使用滑动窗口计算甲基化方差
    3. 选择高方差窗口
    4. 合并重叠的高方差窗口形成VMR
    
    Parameters:
        mdata: MethylationData对象
        bandwidth: 平滑窗口大小(bp)
        stepsize: 滑动窗口步长(bp)
        var_threshold: 方差阈值,超过此值的窗口被认为是VMR
        min_coverage: 最小覆盖度要求
        min_cells: 区域至少在多少个细胞中有覆盖
        chromosomes: 要处理的染色体列表,None表示全部
        n_jobs: 并行线程数,-1表示使用所有CPU
        verbose: 是否显示进度
        inplace: 是否将结果添加到mdata.uns['vmr']
    
    Returns:
        如果inplace=False,返回VMR DataFrame
    """
    from .smooth import smooth_methylation
    
    # 检查是否已经计算平滑曲线
    if 'smooth_methylation' not in mdata.uns:
        if verbose:
            print("计算平滑甲基化曲线...")
        smooth_methylation(mdata, bandwidth=bandwidth, n_jobs=n_jobs, inplace=True)
    
    # 确定要处理的染色体
    if chromosomes is None:
        chromosomes = mdata.chromosomes
    
    if verbose:
        print(f"在 {len(chromosomes)} 条染色体上检测VMR...")
    
    # 并行处理每条染色体
    vmr_list = Parallel(n_jobs=n_jobs)(
        delayed(_detect_vmr_chromosome)(
            mdata,
            chrom=chrom,
            bandwidth=bandwidth,
            stepsize=stepsize,
            var_threshold=var_threshold,
            min_coverage=min_coverage,
            min_cells=min_cells,
            verbose=verbose
        )
        for chrom in (tqdm(chromosomes, desc="Processing chromosomes") if verbose else chromosomes)
    )
    
    # 合并所有染色体的结果
    all_vmrs = pd.concat(vmr_list, ignore_index=True)
    
    if verbose:
        print(f"共检测到 {len(all_vmrs)} 个VMR")
        print(f"VMR总长度: {(all_vmrs['end'] - all_vmrs['start']).sum() / 1e6:.2f} Mb")
    
    if inplace:
        mdata.uns['vmr'] = all_vmrs
        # 标记VMR区域
        mdata.add_vmr(all_vmrs, key='is_vmr')
        return None
    else:
        return all_vmrs


def _detect_vmr_chromosome(
    mdata,
    chrom: str,
    bandwidth: int,
    stepsize: int,
    var_threshold: float,
    min_coverage: int,
    min_cells: int,
    verbose: bool = False
) -> pd.DataFrame:
    """
    在单条染色体上检测VMR
    
    Parameters:
        mdata: MethylationData对象
        chrom: 染色体名称
        其他参数同detect_vmr
    
    Returns:
        该染色体的VMR DataFrame
    """
    # 提取该染色体的数据
    chrom_mask = mdata.var['chrom'] == chrom
    
    if chrom_mask.sum() == 0:
        return pd.DataFrame(columns=['chrom', 'start', 'end', 'variance', 'n_sites'])
    
    # 获取位点坐标
    positions = mdata.var.loc[chrom_mask, 'start'].values
    
    # 获取平滑的甲基化曲线
    smooth_curve = mdata.uns['smooth_methylation'][chrom]
    
    # 获取甲基化率矩阵
    meth_rate = mdata.methylation_rate[:, chrom_mask]
    
    # 获取覆盖度信息
    if 'total' in mdata.layers:
        coverage = mdata.layers['total'][:, chrom_mask]
    else:
        coverage = np.ones_like(meth_rate)
    
    # 过滤低覆盖和低细胞数的位点
    valid_sites = (
        (coverage >= min_coverage).sum(axis=0) >= min_cells
    )
    
    if valid_sites.sum() == 0:
        return pd.DataFrame(columns=['chrom', 'start', 'end', 'variance', 'n_sites'])
    
    # 只保留有效位点
    positions = positions[valid_sites]
    meth_rate = meth_rate[:, valid_sites]
    smooth_curve_interp = np.interp(positions, smooth_curve['position'], smooth_curve['methylation'])
    
    # 计算残差
    residuals = _calculate_residuals(meth_rate, smooth_curve_interp)
    
    # 滑窗计算方差
    windows, variances = _sliding_window_variance(
        positions=positions,
        residuals=residuals,
        bandwidth=bandwidth,
        stepsize=stepsize
    )
    
    if len(windows) == 0:
        return pd.DataFrame(columns=['chrom', 'start', 'end', 'variance', 'n_sites'])
    
    # 过滤高方差窗口
    high_var_mask = variances > var_threshold
    high_var_windows = windows[high_var_mask]
    high_var_values = variances[high_var_mask]
    
    if len(high_var_windows) == 0:
        return pd.DataFrame(columns=['chrom', 'start', 'end', 'variance', 'n_sites'])
    
    # 合并重叠窗口
    vmrs = _merge_overlapping_windows(
        high_var_windows,
        high_var_values,
        chrom=chrom
    )
    
    return vmrs


def _calculate_residuals(
    meth_rate: np.ndarray,
    smooth_curve: np.ndarray,
    shrinkage: float = 0.5
) -> np.ndarray:
    """
    计算收缩残差
    
    Parameters:
        meth_rate: 甲基化率矩阵 (cells × sites)
        smooth_curve: 平滑曲线 (sites,)
        shrinkage: 收缩系数
    
    Returns:
        残差矩阵 (cells × sites)
    """
    # 计算原始残差
    raw_residuals = meth_rate - smooth_curve[np.newaxis, :]
    
    # 应用收缩
    # 收缩朝向0,减少噪声的影响
    residuals = raw_residuals * shrinkage
    
    # 处理NaN值
    residuals = np.nan_to_num(residuals, nan=0.0)
    
    return residuals


def _sliding_window_variance(
    positions: np.ndarray,
    residuals: np.ndarray,
    bandwidth: int,
    stepsize: int
) -> Tuple[np.ndarray, np.ndarray]:
    """
    使用滑动窗口计算方差
    
    Parameters:
        positions: 位点位置数组
        residuals: 残差矩阵 (cells × sites)
        bandwidth: 窗口大小(bp)
        stepsize: 步长(bp)
    
    Returns:
        windows: 窗口坐标数组 (n_windows, 2)
        variances: 方差数组 (n_windows,)
    """
    if len(positions) == 0:
        return np.array([]), np.array([])
    
    # 生成窗口
    chrom_start = positions.min()
    chrom_end = positions.max()
    
    window_starts = np.arange(chrom_start, chrom_end - bandwidth + 1, stepsize)
    window_ends = window_starts + bandwidth
    
    windows = []
    variances = []
    
    for start, end in zip(window_starts, window_ends):
        # 找到窗口内的位点
        in_window = (positions >= start) & (positions < end)
        
        if in_window.sum() < 5:  # 至少需要5个位点
            continue
        
        # 计算窗口内的方差
        window_residuals = residuals[:, in_window]
        
        # 使用所有细胞和位点计算方差
        var = np.nanvar(window_residuals)
        
        windows.append([start, end])
        variances.append(var)
    
    if len(windows) == 0:
        return np.array([]), np.array([])
    
    return np.array(windows), np.array(variances)


def _merge_overlapping_windows(
    windows: np.ndarray,
    variances: np.ndarray,
    chrom: str,
    min_gap: int = 0
) -> pd.DataFrame:
    """
    合并重叠的窗口
    
    Parameters:
        windows: 窗口坐标数组 (n_windows, 2)
        variances: 方差数组 (n_windows,)
        chrom: 染色体名称
        min_gap: 允许的最小间隔,小于此间隔的窗口会被合并
    
    Returns:
        合并后的VMR DataFrame
    """
    if len(windows) == 0:
        return pd.DataFrame(columns=['chrom', 'start', 'end', 'variance', 'n_windows'])
    
    # 按起始位置排序
    sorted_idx = np.argsort(windows[:, 0])
    windows = windows[sorted_idx]
    variances = variances[sorted_idx]
    
    merged = []
    current_start = windows[0, 0]
    current_end = windows[0, 1]
    current_vars = [variances[0]]
    current_count = 1
    
    for i in range(1, len(windows)):
        # 如果当前窗口与前一个重叠或间隔很小
        if windows[i, 0] <= current_end + min_gap:
            # 扩展当前VMR
            current_end = max(current_end, windows[i, 1])
            current_vars.append(variances[i])
            current_count += 1
        else:
            # 保存当前VMR
            merged.append({
                'chrom': chrom,
                'start': int(current_start),
                'end': int(current_end),
                'variance': np.mean(current_vars),
                'n_windows': current_count
            })
            # 开始新的VMR
            current_start = windows[i, 0]
            current_end = windows[i, 1]
            current_vars = [variances[i]]
            current_count = 1
    
    # 添加最后一个VMR
    merged.append({
        'chrom': chrom,
        'start': int(current_start),
        'end': int(current_end),
        'variance': np.mean(current_vars),
        'n_windows': current_count
    })
    
    return pd.DataFrame(merged)


def filter_vmr(
    vmr_df: pd.DataFrame,
    min_length: Optional[int] = None,
    max_length: Optional[int] = None,
    min_variance: Optional[float] = None,
    chromosomes: Optional[list] = None
) -> pd.DataFrame:
    """
    过滤VMR
    
    Parameters:
        vmr_df: VMR DataFrame
        min_length: 最小长度(bp)
        max_length: 最大长度(bp)
        min_variance: 最小方差
        chromosomes: 保留的染色体列表
    
    Returns:
        过滤后的VMR DataFrame
    """
    filtered = vmr_df.copy()
    
    # 计算长度
    filtered['length'] = filtered['end'] - filtered['start']
    
    # 应用过滤条件
    if min_length is not None:
        filtered = filtered[filtered['length'] >= min_length]
    
    if max_length is not None:
        filtered = filtered[filtered['length'] <= max_length]
    
    if min_variance is not None:
        filtered = filtered[filtered['variance'] >= min_variance]
    
    if chromosomes is not None:
        filtered = filtered[filtered['chrom'].isin(chromosomes)]
    
    # 删除临时列
    filtered = filtered.drop(columns=['length'])
    
    return filtered.reset_index(drop=True)


def export_vmr_bed(
    vmr_df: pd.DataFrame,
    filename: str,
    include_score: bool = True
):
    """
    导出VMR为BED格式
    
    Parameters:
        vmr_df: VMR DataFrame
        filename: 输出文件名
        include_score: 是否包含方差作为score列
    """
    bed_df = vmr_df[['chrom', 'start', 'end']].copy()
    
    if include_score:
        # 使用方差作为score,归一化到0-1000
        max_var = vmr_df['variance'].max()
        bed_df['score'] = (vmr_df['variance'] / max_var * 1000).astype(int)
    
    bed_df.to_csv(filename, sep='\t', header=False, index=False)
    print(f"VMR已导出到: {filename}")
