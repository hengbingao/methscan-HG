"""
Multi-Context Methylation Data Structure

Supports storing CG, CH, and other methylation contexts in separate layers
"""

import numpy as np
import pandas as pd
from typing import Optional, Union, List, Dict
from pathlib import Path
import warnings

from ..core.methylation_data import MethylationData


class MultiContextMethylationData(MethylationData):
    """
    Extended MethylationData with multi-context support
    
    Stores different methylation contexts (CG, CH, CHG, CHH, all) in separate layers:
    - layers['CG_rate']: CG methylation rates
    - layers['CG_met']: CG methylated counts
    - layers['CG_total']: CG total coverage
    - layers['CH_rate']: CH methylation rates
    - layers['CH_met']: CH methylated counts
    - layers['CH_total']: CH total coverage
    - etc.
    
    The main X matrix can store any preferred context for visualization.
    """
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        
        # Track which contexts are available
        if 'contexts' not in self.uns:
            self.uns['contexts'] = []
    
    @property
    def available_contexts(self) -> List[str]:
        """Return list of available methylation contexts"""
        return self.uns.get('contexts', [])
    
    def add_context(
        self,
        context: str,
        met: np.ndarray,
        total: np.ndarray,
        overwrite: bool = False
    ):
        """
        Add a methylation context layer
        
        Parameters:
            context: Context name (e.g., 'CG', 'CH', 'CHG')
            met: Methylated counts matrix (cells × regions)
            total: Total coverage matrix (cells × regions)
            overwrite: Overwrite existing context
        """
        if context in self.available_contexts and not overwrite:
            raise ValueError(
                f"Context '{context}' already exists. "
                f"Use overwrite=True to replace."
            )
        
        # Calculate methylation rate
        with np.errstate(divide='ignore', invalid='ignore'):
            rate = met.astype(float) / total
            rate[~np.isfinite(rate)] = np.nan
        
        # Add to layers
        self.layers[f'{context}_met'] = met
        self.layers[f'{context}_total'] = total
        self.layers[f'{context}_rate'] = rate
        
        # Update context list
        if context not in self.available_contexts:
            self.uns['contexts'].append(context)
        
        print(f"Added context '{context}' with shape {met.shape}")
    
    def get_context(
        self,
        context: str,
        data_type: str = 'rate'
    ) -> np.ndarray:
        """
        Get data for a specific context
        
        Parameters:
            context: Context name
            data_type: 'rate', 'met', or 'total'
            
        Returns:
            Data matrix for the specified context
        """
        if context not in self.available_contexts:
            raise ValueError(
                f"Context '{context}' not found. "
                f"Available: {self.available_contexts}"
            )
        
        layer_name = f'{context}_{data_type}'
        
        if layer_name not in self.layers:
            raise ValueError(
                f"Layer '{layer_name}' not found. "
                f"Valid data_type: 'rate', 'met', 'total'"
            )
        
        return self.layers[layer_name]
    
    def set_active_context(self, context: str):
        """
        Set the active context for X matrix
        
        This determines which context is used for default operations
        like PCA, UMAP, clustering, etc.
        
        Parameters:
            context: Context to set as active
        """
        if context not in self.available_contexts:
            raise ValueError(
                f"Context '{context}' not available. "
                f"Available: {self.available_contexts}"
            )
        
        # Set X to this context's rate
        self.X = self.get_context(context, 'rate')
        self.uns['active_context'] = context
        
        print(f"Set '{context}' as active context for analysis")
    
    def calculate_context_correlations(self) -> pd.DataFrame:
        """
        Calculate correlations between different contexts
        
        Returns:
            DataFrame with pairwise correlations
        """
        correlations = []
        
        contexts = self.available_contexts
        
        for i, ctx1 in enumerate(contexts):
            for ctx2 in contexts[i:]:
                # Get methylation rates
                rate1 = self.get_context(ctx1, 'rate')
                rate2 = self.get_context(ctx2, 'rate')
                
                # Calculate correlation (cell-wise)
                corrs = []
                for cell_idx in range(self.n_obs):
                    r1 = rate1[cell_idx, :]
                    r2 = rate2[cell_idx, :]
                    
                    # Remove NaNs
                    valid = ~(np.isnan(r1) | np.isnan(r2))
                    if valid.sum() > 10:
                        corr = np.corrcoef(r1[valid], r2[valid])[0, 1]
                        corrs.append(corr)
                
                if len(corrs) > 0:
                    correlations.append({
                        'context1': ctx1,
                        'context2': ctx2,
                        'mean_correlation': np.mean(corrs),
                        'std_correlation': np.std(corrs)
                    })
        
        return pd.DataFrame(correlations)
    
    def plot_context_comparison(
        self,
        contexts: Optional[List[str]] = None,
        metric: str = 'mean',
        save: Optional[str] = None
    ):
        """
        Plot comparison of methylation levels across contexts
        
        Parameters:
            contexts: Contexts to compare (None = all)
            metric: 'mean', 'median', 'variance'
            save: Save path
        """
        import matplotlib.pyplot as plt
        import seaborn as sns
        
        if contexts is None:
            contexts = self.available_contexts
        
        # Calculate metric for each context
        data = []
        for ctx in contexts:
            rate = self.get_context(ctx, 'rate')
            
            if metric == 'mean':
                values = np.nanmean(rate, axis=1)
            elif metric == 'median':
                values = np.nanmedian(rate, axis=1)
            elif metric == 'variance':
                values = np.nanvar(rate, axis=1)
            else:
                raise ValueError(f"Unknown metric: {metric}")
            
            for val in values:
                data.append({
                    'context': ctx,
                    'value': val
                })
        
        df = pd.DataFrame(data)
        
        # Plot
        fig, ax = plt.subplots(figsize=(8, 6))
        
        sns.violinplot(
            data=df,
            x='context',
            y='value',
            ax=ax
        )
        
        ax.set_xlabel('Methylation Context')
        ax.set_ylabel(f'{metric.capitalize()} Methylation Level')
        ax.set_title(f'{metric.capitalize()} Methylation by Context')
        
        plt.tight_layout()
        
        if save:
            plt.savefig(save, dpi=300, bbox_inches='tight')
            print(f"Saved to {save}")
        
        return fig
    
    def __repr__(self) -> str:
        """Enhanced string representation"""
        descr = super().__repr__()
        
        if self.available_contexts:
            descr += f"\n    contexts: {', '.join(self.available_contexts)}"
            if 'active_context' in self.uns:
                descr += f" (active: {self.uns['active_context']})"
        
        return descr


def create_multi_context_data_from_allc(
    allc_files: Union[str, List[str]],
    regions: Optional[pd.DataFrame] = None,
    contexts: List[str] = ['CG', 'CH', 'all'],
    region_type: str = 'genome_wide',
    min_coverage: int = 1,
    n_jobs: int = -1,
    verbose: bool = True
) -> MultiContextMethylationData:
    """
    Create MultiContextMethylationData directly from ALLC files
    
    This function:
    1. Reads ALLC files
    2. Separates contexts (CG, CH, etc.)
    3. Quantifies methylation across regions
    4. Creates a multi-layer MethylationData object
    
    Parameters:
        allc_files: ALLC file path(s)
        regions: Genomic regions to quantify (DataFrame with chr, start, end)
                If None, uses genome-wide binning
        contexts: Methylation contexts to extract
        region_type: Name for region type
        min_coverage: Minimum coverage threshold
        n_jobs: Number of parallel jobs
        verbose: Show progress
        
    Returns:
        MultiContextMethylationData object with all contexts in separate layers
    """
    from glob import glob
    from joblib import Parallel, delayed
    from tqdm.auto import tqdm
    
    # Parse file paths
    if isinstance(allc_files, str):
        files = sorted(glob(allc_files))
    else:
        files = allc_files
    
    if len(files) == 0:
        raise ValueError(f"No ALLC files found: {allc_files}")
    
    if verbose:
        print(f"Creating multi-context data from {len(files)} ALLC files")
        print(f"Contexts: {contexts}")
    
    # If no regions specified, create genome-wide bins
    if regions is None:
        if verbose:
            print("No regions specified, using 100kb bins")
        regions = _create_genome_bins(bin_size=100000)
    
    # Process each file for each context
    results = Parallel(n_jobs=n_jobs)(
        delayed(_process_allc_multi_context)(
            filepath=f,
            regions=regions,
            contexts=contexts,
            min_coverage=min_coverage
        )
        for f in (tqdm(files, desc="Processing") if verbose else files)
    )
    
    # Combine results into matrices
    n_cells = len(files)
    n_regions = len(regions)
    
    # Initialize matrices for each context
    context_data = {}
    for ctx in contexts:
        context_data[ctx] = {
            'met': np.zeros((n_cells, n_regions), dtype=np.int32),
            'total': np.zeros((n_cells, n_regions), dtype=np.int32)
        }
    
    # Fill matrices
    sample_names = []
    for cell_idx, (cell_name, cell_data) in enumerate(results):
        sample_names.append(cell_name)
        
        for ctx in contexts:
            if ctx in cell_data:
                context_data[ctx]['met'][cell_idx, :] = cell_data[ctx]['met']
                context_data[ctx]['total'][cell_idx, :] = cell_data[ctx]['total']
    
    # Create obs (cell metadata)
    obs = pd.DataFrame({'sample_name': sample_names})
    obs.index = sample_names
    
    # Create var (region metadata)
    var = regions.copy()
    var.index = [f"{r.chr}:{r.start}-{r.end}" for _, r in regions.iterrows()]
    
    # Create MultiContextMethylationData
    # Use first context as default X
    first_ctx = contexts[0]
    with np.errstate(divide='ignore', invalid='ignore'):
        X = context_data[first_ctx]['met'].astype(float) / context_data[first_ctx]['total']
        X[~np.isfinite(X)] = np.nan
    
    mdata = MultiContextMethylationData(
        X=X,
        obs=obs,
        var=var,
        uns={'region_type': region_type, 'active_context': first_ctx}
    )
    
    # Add all contexts as layers
    for ctx in contexts:
        mdata.add_context(
            context=ctx,
            met=context_data[ctx]['met'],
            total=context_data[ctx]['total']
        )
    
    if verbose:
        print(f"\nCreated MultiContextMethylationData:")
        print(f"  Shape: {mdata.n_obs} cells × {mdata.n_vars} regions")
        print(f"  Contexts: {', '.join(mdata.available_contexts)}")
        print(f"  Active context: {mdata.uns['active_context']}")
    
    return mdata


def _process_allc_multi_context(
    filepath: str,
    regions: pd.DataFrame,
    contexts: List[str],
    min_coverage: int
) -> tuple:
    """Process single ALLC file for multiple contexts"""
    import gzip
    
    cell_name = Path(filepath).stem
    
    # Read ALLC
    if filepath.endswith('.gz'):
        with gzip.open(filepath, 'rt') as f:
            df = pd.read_csv(
                f, sep='\t', header=None,
                names=['chr', 'pos', 'strand', 'context', 'mc', 'cov', 'meth']
            )
    else:
        df = pd.read_csv(
            filepath, sep='\t', header=None,
            names=['chr', 'pos', 'strand', 'context', 'mc', 'cov', 'meth']
        )
    
    # Filter by coverage
    df = df[df['cov'] >= min_coverage]
    
    # Separate contexts
    context_dfs = {}
    for ctx in contexts:
        if ctx == 'CG':
            context_dfs[ctx] = df[df['context'].str.startswith('CG')]
        elif ctx == 'CH':
            context_dfs[ctx] = df[df['context'].str.match(r'^C[ATC]')]
        elif ctx == 'CHG':
            context_dfs[ctx] = df[df['context'].str.startswith('CHG')]
        elif ctx == 'CHH':
            context_dfs[ctx] = df[df['context'].str.startswith('CHH')]
        else:  # 'all'
            context_dfs[ctx] = df.copy()
    
    # Quantify each context across regions
    cell_data = {}
    
    for ctx, ctx_df in context_dfs.items():
        met_counts = np.zeros(len(regions), dtype=np.int32)
        total_counts = np.zeros(len(regions), dtype=np.int32)
        
        # For each region, sum overlapping sites
        for region_idx, region in regions.iterrows():
            overlaps = ctx_df[
                (ctx_df['chr'] == region['chr']) &
                (ctx_df['pos'] >= region['start']) &
                (ctx_df['pos'] < region['end'])
            ]
            
            if len(overlaps) > 0:
                met_counts[region_idx] = overlaps['mc'].sum()
                total_counts[region_idx] = overlaps['cov'].sum()
        
        cell_data[ctx] = {
            'met': met_counts,
            'total': total_counts
        }
    
    return (cell_name, cell_data)


def _create_genome_bins(bin_size: int = 100000, genome: str = 'hg38') -> pd.DataFrame:
    """Create genome-wide bins for quantification"""
    # Simplified chromosome sizes (hg38)
    chrom_sizes = {
        'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559,
        'chr4': 190214555, 'chr5': 181538259, 'chr6': 170805979,
        'chr7': 159345973, 'chr8': 145138636, 'chr9': 138394717,
        'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
        'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189,
        'chr16': 90338345, 'chr17': 83257441, 'chr18': 80373285,
        'chr19': 58617616, 'chr20': 64444167, 'chr21': 46709983,
        'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415
    }
    
    bins = []
    for chrom, size in chrom_sizes.items():
        for start in range(0, size, bin_size):
            end = min(start + bin_size, size)
            bins.append({'chr': chrom, 'start': start, 'end': end})
    
    return pd.DataFrame(bins)
