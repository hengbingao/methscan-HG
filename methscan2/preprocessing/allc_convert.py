"""
ALLC File Processing and Context Extraction

This module provides functions to:
1. Convert ALLC format to COV format
2. Extract specific methylation contexts (CG, CH, CHG, CHH)
3. Calculate global methylation levels
4. Process multiple files in parallel
"""

import numpy as np
import pandas as pd
import gzip
from pathlib import Path
from typing import Optional, Union, List, Dict
from joblib import Parallel, delayed
from tqdm.auto import tqdm
import subprocess
import os


def allc_to_cov(
    allc_files: Union[str, List[str]],
    output_dir: str = 'cov_files',
    contexts: Optional[List[str]] = None,
    separate_contexts: bool = True,
    compress: bool = False,
    n_jobs: int = -1,
    verbose: bool = True
) -> Dict[str, List[str]]:
    """
    Convert ALLC files to COV format with context separation
    
    ALLC format: chr pos strand context mc_count total_count methylation_level
    COV format:  chr start end meth_pct met_count unmet_count
    
    Parameters:
        allc_files: Path pattern or list of ALLC files
        output_dir: Output directory for COV files
        contexts: List of contexts to extract ['CG', 'CH', 'all']. 
                 None means extract all contexts separately
        separate_contexts: If True, create separate files for each context
        compress: Compress output files with gzip
        n_jobs: Number of parallel jobs
        verbose: Show progress
        
    Returns:
        Dictionary with context names as keys and output file lists as values
    """
    from glob import glob
    
    # Parse input files
    if isinstance(allc_files, str):
        files = sorted(glob(allc_files))
    else:
        files = allc_files
    
    if len(files) == 0:
        raise ValueError(f"No files found: {allc_files}")
    
    # Default contexts
    if contexts is None:
        contexts = ['CG', 'CH', 'all'] if separate_contexts else ['all']
    
    # Create output directory
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    if verbose:
        print(f"Converting {len(files)} ALLC files to COV format...")
        print(f"Contexts: {contexts}")
        print(f"Output directory: {output_dir}")
    
    # Process files in parallel
    results = Parallel(n_jobs=n_jobs)(
        delayed(_process_single_allc)(
            filepath=f,
            output_dir=output_dir,
            contexts=contexts,
            compress=compress,
            verbose=False
        )
        for f in (tqdm(files, desc="Converting") if verbose else files)
    )
    
    # Organize results by context
    output_files = {ctx: [] for ctx in contexts}
    for result in results:
        for ctx, filepath in result.items():
            output_files[ctx].append(filepath)
    
    if verbose:
        print("\nConversion complete!")
        for ctx, files in output_files.items():
            print(f"  {ctx}: {len(files)} files")
    
    return output_files


def _process_single_allc(
    filepath: str,
    output_dir: str,
    contexts: List[str],
    compress: bool,
    verbose: bool
) -> Dict[str, str]:
    """Process a single ALLC file and convert to COV format(s)"""
    
    basename = Path(filepath).stem
    if basename.endswith('.allc'):
        basename = basename[:-5]
    
    # Read ALLC file
    if filepath.endswith('.gz'):
        with gzip.open(filepath, 'rt') as f:
            df = pd.read_csv(
                f,
                sep='\t',
                header=None,
                names=['chr', 'pos', 'strand', 'context', 'mc', 'cov', 'meth_level']
            )
    else:
        df = pd.read_csv(
            filepath,
            sep='\t',
            header=None,
            names=['chr', 'pos', 'strand', 'context', 'mc', 'cov', 'meth_level']
        )
    
    output_files = {}
    
    # Process each context
    for ctx in contexts:
        # Filter by context
        if ctx == 'all':
            subset = df.copy()
            suffix = ''
        elif ctx == 'CG':
            subset = df[df['context'].str.startswith('CG')].copy()
            suffix = '_CG'
        elif ctx == 'CH':
            subset = df[df['context'].str.match(r'^C[ATC]')].copy()
            suffix = '_CH'
        elif ctx == 'CHG':
            subset = df[df['context'].str.startswith('CHG')].copy()
            suffix = '_CHG'
        elif ctx == 'CHH':
            subset = df[df['context'].str.startswith('CHH')].copy()
            suffix = '_CHH'
        else:
            # Custom context pattern
            subset = df[df['context'].str.startswith(ctx)].copy()
            suffix = f'_{ctx}'
        
        if len(subset) == 0:
            continue
        
        # Convert to COV format
        cov_df = pd.DataFrame({
            'chr': subset['chr'],
            'start': subset['pos'],
            'end': subset['pos'] + 1,
            'meth_pct': (subset['mc'] > 0).astype(int) * 100,  # 0 or 100
            'met_count': subset['mc'],
            'unmet_count': subset['cov'] - subset['mc']
        })
        
        # Output filename
        output_file = Path(output_dir) / f"{basename}{suffix}.cov"
        if compress:
            output_file = Path(str(output_file) + '.gz')
        
        # Write file
        if compress:
            with gzip.open(output_file, 'wt') as f:
                cov_df.to_csv(f, sep='\t', header=False, index=False)
        else:
            cov_df.to_csv(output_file, sep='\t', header=False, index=False)
        
        output_files[ctx] = str(output_file)
    
    return output_files


def extract_context_from_allc(
    allc_file: str,
    contexts: List[str] = ['CG', 'CH'],
    output_prefix: Optional[str] = None,
    chrom_sizes: Optional[str] = None,
    compress: bool = True,
    cpu: int = 1
) -> List[str]:
    """
    Extract specific methylation contexts from ALLC file
    
    This is a wrapper around allcools extract-allc command.
    Requires allcools to be installed.
    
    Parameters:
        allc_file: Path to ALLC file (.gz)
        contexts: List of contexts to extract (e.g., ['CG', 'CH'])
        output_prefix: Output file prefix
        chrom_sizes: Path to chromosome sizes file
        compress: Compress output
        cpu: Number of CPUs
        
    Returns:
        List of output files
    """
    if output_prefix is None:
        output_prefix = Path(allc_file).stem.replace('.allc', '')
    
    output_files = []
    
    for ctx in contexts:
        # Convert context format
        if ctx == 'CG':
            mc_context = 'CGN'
        elif ctx == 'CH':
            mc_context = 'CHN'
        elif ctx == 'CHG':
            mc_context = 'CHG'
        elif ctx == 'CHH':
            mc_context = 'CHH'
        else:
            mc_context = ctx
        
        output_file = f"{output_prefix}_{ctx}"
        
        # Build command
        cmd = [
            'allcools', 'extract-allc',
            '--allc_path', allc_file,
            '--output_prefix', output_file,
            '--mc_contexts', mc_context,
            '--cpu', str(cpu)
        ]
        
        if chrom_sizes:
            cmd.extend(['--chrom_size_path', chrom_sizes])
        
        # Run command
        try:
            subprocess.run(cmd, check=True, capture_output=True)
            output_files.append(f"{output_file}.allc.tsv.gz")
        except subprocess.CalledProcessError as e:
            print(f"Warning: Failed to extract {ctx} from {allc_file}")
            print(f"Error: {e.stderr.decode()}")
        except FileNotFoundError:
            raise ImportError(
                "allcools not found. Please install: pip install allcools"
            )
    
    return output_files


def calculate_global_methylation(
    cov_files: Union[str, List[str]],
    output_file: str = 'global_methylation_levels.tsv',
    context_name: str = 'mC',
    n_jobs: int = -1,
    verbose: bool = True
) -> pd.DataFrame:
    """
    Calculate global methylation levels for each cell
    
    Parameters:
        cov_files: Path pattern or list of COV files
        output_file: Output file path
        context_name: Context name (e.g., 'CG', 'CH', 'mC')
        n_jobs: Number of parallel jobs
        verbose: Show progress
        
    Returns:
        DataFrame with columns: cell, global_level, n_sites
    """
    from glob import glob
    
    # Parse files
    if isinstance(cov_files, str):
        files = sorted(glob(cov_files))
    else:
        files = cov_files
    
    if len(files) == 0:
        raise ValueError(f"No files found: {cov_files}")
    
    if verbose:
        print(f"Calculating global {context_name} methylation for {len(files)} cells...")
    
    # Process files in parallel
    results = Parallel(n_jobs=n_jobs)(
        delayed(_calculate_single_global)(f)
        for f in (tqdm(files, desc="Processing") if verbose else files)
    )
    
    # Create DataFrame
    df = pd.DataFrame(results, columns=['cell', f'global_{context_name}_level', f'{context_name}_sites'])
    
    # Save to file
    df.to_csv(output_file, sep='\t', index=False)
    
    if verbose:
        print(f"\nGlobal methylation statistics:")
        print(f"  Mean {context_name} level: {df[f'global_{context_name}_level'].mean():.3f}")
        print(f"  Median {context_name} level: {df[f'global_{context_name}_level'].median():.3f}")
        print(f"  Mean sites per cell: {df[f'{context_name}_sites'].mean():.0f}")
        print(f"\nResults saved to: {output_file}")
    
    return df


def _calculate_single_global(cov_file: str) -> tuple:
    """Calculate global methylation for a single COV file"""
    
    cell_name = Path(cov_file).stem
    if cell_name.endswith('.cov'):
        cell_name = cell_name[:-4]
    
    try:
        # Read COV file
        if cov_file.endswith('.gz'):
            df = pd.read_csv(
                cov_file,
                sep='\t',
                header=None,
                names=['chr', 'start', 'end', 'meth_pct', 'met', 'unmet']
            )
        else:
            df = pd.read_csv(
                cov_file,
                sep='\t',
                header=None,
                names=['chr', 'start', 'end', 'meth_pct', 'met', 'unmet']
            )
        
        # Calculate global methylation level
        # meth_pct is 0 or 100, so average gives global level
        global_level = df['meth_pct'].mean() / 100.0
        n_sites = len(df)
        
        return (cell_name, global_level, n_sites)
    
    except Exception as e:
        print(f"Warning: Failed to process {cov_file}: {e}")
        return (cell_name, np.nan, 0)


def batch_process_allc(
    allc_dir: str,
    output_base_dir: str,
    contexts: List[str] = ['CG', 'CH', 'all'],
    calculate_global: bool = True,
    n_jobs: int = -1,
    verbose: bool = True
) -> Dict[str, pd.DataFrame]:
    """
    Complete pipeline: ALLC → COV → Global methylation levels
    
    This replicates your bash workflow in Python with parallelization.
    
    Parameters:
        allc_dir: Directory containing ALLC files
        output_base_dir: Base output directory
        contexts: Contexts to extract
        calculate_global: Calculate global methylation levels
        n_jobs: Number of parallel jobs
        verbose: Show progress
        
    Returns:
        Dictionary of global methylation DataFrames for each context
    """
    from glob import glob
    
    allc_files = glob(os.path.join(allc_dir, '*.allc.tsv.gz'))
    if len(allc_files) == 0:
        allc_files = glob(os.path.join(allc_dir, '*.allc.gz'))
    if len(allc_files) == 0:
        raise ValueError(f"No ALLC files found in {allc_dir}")
    
    if verbose:
        print("=" * 60)
        print("ALLC Batch Processing Pipeline")
        print("=" * 60)
        print(f"Input directory: {allc_dir}")
        print(f"Found {len(allc_files)} ALLC files")
        print(f"Contexts: {contexts}")
        print(f"Output directory: {output_base_dir}")
    
    # Create output directories
    Path(output_base_dir).mkdir(parents=True, exist_ok=True)
    
    # Step 1: Convert ALLC to COV
    if verbose:
        print("\n[1/2] Converting ALLC to COV format...")
    
    output_files = allc_to_cov(
        allc_files=allc_files,
        output_dir=output_base_dir,
        contexts=contexts,
        separate_contexts=True,
        compress=False,
        n_jobs=n_jobs,
        verbose=verbose
    )
    
    # Step 2: Organize files into context-specific directories
    context_dirs = {}
    for ctx in contexts:
        ctx_dir = os.path.join(output_base_dir, f'{ctx}_cov')
        Path(ctx_dir).mkdir(exist_ok=True)
        context_dirs[ctx] = ctx_dir
        
        # Move files
        for f in output_files[ctx]:
            dest = os.path.join(ctx_dir, Path(f).name)
            if f != dest:
                os.rename(f, dest)
    
    if verbose:
        print("\n[2/2] Calculating global methylation levels...")
    
    # Step 3: Calculate global methylation for each context
    global_stats = {}
    
    if calculate_global:
        global_dir = os.path.join(output_base_dir, 'global_level')
        Path(global_dir).mkdir(exist_ok=True)
        
        for ctx in contexts:
            ctx_dir = context_dirs[ctx]
            cov_pattern = os.path.join(ctx_dir, '*.cov')
            
            output_file = os.path.join(
                global_dir,
                f'all_cell_global_average_m{ctx}.tsv'
            )
            
            if verbose:
                print(f"\n  Processing {ctx} context...")
            
            df = calculate_global_methylation(
                cov_files=cov_pattern,
                output_file=output_file,
                context_name=ctx,
                n_jobs=n_jobs,
                verbose=verbose
            )
            
            global_stats[ctx] = df
    
    if verbose:
        print("\n" + "=" * 60)
        print("Pipeline Complete!")
        print("=" * 60)
        print(f"\nOutput structure:")
        print(f"  {output_base_dir}/")
        for ctx in contexts:
            print(f"    ├── {ctx}_cov/        ({len(output_files[ctx])} files)")
        print(f"    └── global_level/")
        for ctx in contexts:
            print(f"          └── all_cell_global_average_m{ctx}.tsv")
    
    return global_stats


def merge_global_stats(
    stats_files: Union[List[str], Dict[str, str]],
    output_file: str = 'merged_global_stats.tsv'
) -> pd.DataFrame:
    """
    Merge global methylation statistics from different contexts
    
    Parameters:
        stats_files: List of or dict of global stats files
        output_file: Output merged file
        
    Returns:
        Merged DataFrame
    """
    if isinstance(stats_files, dict):
        dfs = []
        for ctx, filepath in stats_files.items():
            df = pd.read_csv(filepath, sep='\t')
            dfs.append(df)
    else:
        dfs = [pd.read_csv(f, sep='\t') for f in stats_files]
    
    # Merge on cell name
    merged = dfs[0]
    for df in dfs[1:]:
        merged = merged.merge(df, on='cell', how='outer')
    
    # Save
    merged.to_csv(output_file, sep='\t', index=False)
    
    return merged
