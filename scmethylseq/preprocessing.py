"""
Preprocessing module for allc files
Handles allc file processing, CG/CH extraction, and cov file generation
"""

import os
import gzip
import subprocess
from pathlib import Path
from typing import Optional, List, Tuple
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from functools import partial
import logging
import pandas as pd
import numpy as np
from tqdm import tqdm

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class AllcProcessor:
    """
    Process allc files to extract CG and CH methylation contexts
    """
    
    def __init__(self, 
                 allc_dir: str,
                 output_dir: str,
                 chrom_sizes: str,
                 n_jobs: int = 10):
        """
        Initialize AllcProcessor
        
        Parameters:
        -----------
        allc_dir : str
            Directory containing allc.gz files
        output_dir : str
            Output directory for processed files
        chrom_sizes : str
            Path to chromosome sizes file
        n_jobs : int
            Number of parallel processes
        """
        self.allc_dir = Path(allc_dir)
        self.output_dir = Path(output_dir)
        self.chrom_sizes = chrom_sizes
        self.n_jobs = n_jobs
        
        # Create output directories
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
    def _process_single_allc(self, allc_file: Path) -> Tuple[bool, str]:
        """
        Process a single allc file to extract CG and CH contexts
        
        Parameters:
        -----------
        allc_file : Path
            Path to allc.gz file
            
        Returns:
        --------
        Tuple[bool, str]
            Success status and message
        """
        try:
            cell_name = allc_file.stem.replace('.allc', '')
            
            # Tabix index the file
            tabix_cmd = f"tabix -p bed {allc_file}"
            subprocess.run(tabix_cmd, shell=True, check=True, 
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            
            # Extract CG context
            cg_output = self.output_dir / f"{cell_name}_CG.allc.gz"
            cg_cmd = (f"allcools extract-allc --allc_path {allc_file} "
                     f"--output_prefix {self.output_dir / cell_name}_CG "
                     f"--mc_contexts CGN "
                     f"--chrom_size_path {self.chrom_sizes} "
                     f"--cpu 1")
            subprocess.run(cg_cmd, shell=True, check=True,
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            
            # Extract CH context
            ch_output = self.output_dir / f"{cell_name}_CH.allc.gz"
            ch_cmd = (f"allcools extract-allc --allc_path {allc_file} "
                     f"--output_prefix {self.output_dir / cell_name}_CH "
                     f"--mc_contexts CHN "
                     f"--chrom_size_path {self.chrom_sizes} "
                     f"--cpu 1")
            subprocess.run(ch_cmd, shell=True, check=True,
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            
            return True, f"Successfully processed {cell_name}"
            
        except Exception as e:
            return False, f"Error processing {allc_file.name}: {str(e)}"
    
    def process_all(self) -> None:
        """
        Process all allc files in parallel
        """
        allc_files = list(self.allc_dir.glob("*.allc.gz"))
        logger.info(f"Found {len(allc_files)} allc files to process")
        
        with ProcessPoolExecutor(max_workers=self.n_jobs) as executor:
            futures = {executor.submit(self._process_single_allc, f): f 
                      for f in allc_files}
            
            with tqdm(total=len(allc_files), desc="Processing allc files") as pbar:
                for future in as_completed(futures):
                    success, message = future.result()
                    if success:
                        logger.info(message)
                    else:
                        logger.error(message)
                    pbar.update(1)


class CovFileGenerator:
    """
    Convert allc files to cov format
    """
    
    def __init__(self, 
                 allc_dir: str,
                 output_dir: str,
                 n_jobs: int = 20):
        """
        Initialize CovFileGenerator
        
        Parameters:
        -----------
        allc_dir : str
            Directory containing allc.gz files
        output_dir : str
            Output directory for cov files
        n_jobs : int
            Number of parallel processes
        """
        self.allc_dir = Path(allc_dir)
        self.output_dir = Path(output_dir)
        self.n_jobs = n_jobs
        
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def _allc_to_cov(self, allc_file: Path) -> Tuple[bool, str]:
        """
        Convert single allc file to cov format
        
        Format: chr  start  end  methylation_percent  mc_count  cov_count
        """
        try:
            cell_name = allc_file.stem.replace('.allc', '')
            output_file = self.output_dir / f"{cell_name}.cov"
            
            # Read allc file and convert
            with gzip.open(allc_file, 'rt') as fin, open(output_file, 'w') as fout:
                for line in fin:
                    if line.startswith('#'):
                        continue
                    
                    parts = line.strip().split('\t')
                    if len(parts) < 6:
                        continue
                    
                    chrom = parts[0]
                    pos = parts[1]
                    mc = int(parts[4])
                    cov = int(parts[5])
                    
                    # Calculate methylation percentage (0-100 scale)
                    if mc == 0:
                        meth_pct = 0
                    else:
                        meth_pct = 100
                    
                    # Write in cov format
                    fout.write(f"{chrom}\t{pos}\t{pos}\t{meth_pct}\t{mc}\t{cov}\n")
            
            return True, f"Converted {cell_name} to cov format"
            
        except Exception as e:
            return False, f"Error converting {allc_file.name}: {str(e)}"
    
    def convert_all(self, pattern: str = "*.allc.gz") -> None:
        """
        Convert all allc files to cov format in parallel
        
        Parameters:
        -----------
        pattern : str
            File pattern to match (e.g., "*_CG.allc.gz", "*_CH.allc.gz")
        """
        allc_files = list(self.allc_dir.glob(pattern))
        logger.info(f"Found {len(allc_files)} files to convert")
        
        with ThreadPoolExecutor(max_workers=self.n_jobs) as executor:
            futures = {executor.submit(self._allc_to_cov, f): f 
                      for f in allc_files}
            
            with tqdm(total=len(allc_files), desc="Converting to cov") as pbar:
                for future in as_completed(futures):
                    success, message = future.result()
                    if success:
                        logger.info(message)
                    else:
                        logger.error(message)
                    pbar.update(1)


class GlobalMethylationCalculator:
    """
    Calculate global CG and CH methylation levels for each cell
    """
    
    def __init__(self, 
                 cov_dir: str,
                 output_file: str,
                 context: str = "CG",
                 n_jobs: int = 20):
        """
        Initialize GlobalMethylationCalculator
        
        Parameters:
        -----------
        cov_dir : str
            Directory containing cov files
        output_file : str
            Output file path
        context : str
            Methylation context ("CG" or "CH")
        n_jobs : int
            Number of parallel processes
        """
        self.cov_dir = Path(cov_dir)
        self.output_file = output_file
        self.context = context
        self.n_jobs = n_jobs
    
    def _calculate_global_level(self, cov_file: Path) -> Tuple[str, float, int]:
        """
        Calculate global methylation level for a single cell
        
        Returns:
        --------
        Tuple[str, float, int]
            Cell name, methylation level (0-1), number of sites
        """
        try:
            cell_name = cov_file.stem.replace('.cov', '')
            
            # Read cov file
            df = pd.read_csv(cov_file, sep='\t', header=None,
                           names=['chr', 'start', 'end', 'meth_pct', 'mc', 'cov'])
            
            # Calculate global methylation level
            total_sites = len(df)
            avg_methylation = df['meth_pct'].mean() / 100  # Convert to 0-1 scale
            
            return cell_name, avg_methylation, total_sites
            
        except Exception as e:
            logger.error(f"Error processing {cov_file.name}: {str(e)}")
            return None, None, None
    
    def calculate_all(self) -> pd.DataFrame:
        """
        Calculate global methylation for all cells in parallel
        
        Returns:
        --------
        pd.DataFrame
            DataFrame with columns: cell, global_level, n_sites
        """
        cov_files = list(self.cov_dir.glob("*.cov"))
        logger.info(f"Calculating global {self.context} levels for {len(cov_files)} cells")
        
        results = []
        with ThreadPoolExecutor(max_workers=self.n_jobs) as executor:
            futures = {executor.submit(self._calculate_global_level, f): f 
                      for f in cov_files}
            
            with tqdm(total=len(cov_files), desc=f"Calculating {self.context}") as pbar:
                for future in as_completed(futures):
                    cell, level, sites = future.result()
                    if cell is not None:
                        results.append({
                            'cell': cell,
                            f'global_{self.context}_level': level,
                            f'{self.context}_sites': sites
                        })
                    pbar.update(1)
        
        # Create DataFrame
        df = pd.DataFrame(results)
        df = df.dropna()
        
        # Save results
        df.to_csv(self.output_file, sep='\t', index=False)
        logger.info(f"Saved global {self.context} levels to {self.output_file}")
        
        return df


class QualityControl:
    """
    Perform quality control on cells based on methylation metrics
    """
    
    def __init__(self, 
                 cg_file: str,
                 ch_file: str,
                 output_dir: str):
        """
        Initialize QualityControl
        
        Parameters:
        -----------
        cg_file : str
            File with global CG metrics
        ch_file : str
            File with global CH metrics
        output_dir : str
            Output directory for QC results
        """
        self.cg_file = cg_file
        self.ch_file = ch_file
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def filter_cells(self,
                    cg_min: float = 0.5,
                    cg_sites_min: int = 100000,
                    cg_sites_max: int = 2500000,
                    ch_max: float = 0.1,
                    ch_sites_min: int = 100000) -> pd.DataFrame:
        """
        Filter cells based on QC metrics
        
        Parameters:
        -----------
        cg_min : float
            Minimum global CG methylation level
        cg_sites_min : int
            Minimum number of CG sites
        cg_sites_max : int
            Maximum number of CG sites
        ch_max : float
            Maximum global CH methylation level
        ch_sites_min : int
            Minimum number of CH sites
            
        Returns:
        --------
        pd.DataFrame
            Filtered cell metadata
        """
        # Load data
        cg_df = pd.read_csv(self.cg_file, sep='\t')
        ch_df = pd.read_csv(self.ch_file, sep='\t')
        
        # Merge CG and CH data
        all_df = pd.merge(cg_df, ch_df, on='cell', how='inner')
        
        logger.info(f"Total cells before filtering: {len(all_df)}")
        
        # Apply filters
        filtered_df = all_df[
            (all_df['global_CG_level'] >= cg_min) &
            (all_df['CG_sites'] >= cg_sites_min) &
            (all_df['CG_sites'] <= cg_sites_max) &
            (all_df['global_CH_level'] <= ch_max) &
            (all_df['CH_sites'] >= ch_sites_min)
        ]
        
        logger.info(f"Cells after filtering: {len(filtered_df)}")
        logger.info(f"Filtered out: {len(all_df) - len(filtered_df)} cells")
        
        # Save results
        all_df.to_csv(self.output_dir / "cell_raw.csv", index=False)
        filtered_df.to_csv(self.output_dir / "cell_filtered.csv", index=False)
        
        return filtered_df
    
    def plot_qc_metrics(self, raw_df: pd.DataFrame, filtered_df: pd.DataFrame) -> None:
        """
        Generate QC plots
        
        Parameters:
        -----------
        raw_df : pd.DataFrame
            Raw cell metrics
        filtered_df : pd.DataFrame
            Filtered cell metrics
        """
        try:
            import matplotlib.pyplot as plt
            import seaborn as sns
            
            fig, axes = plt.subplots(2, 4, figsize=(16, 8))
            
            # CG plots - raw
            axes[0, 0].hist(raw_df['global_CG_level'], bins=50, color='black')
            axes[0, 0].axvline(x=0.5, color='red', linestyle='--')
            axes[0, 0].set_title('CG Level Distribution (Raw)')
            axes[0, 0].set_xlabel('Global CG Level')
            
            axes[0, 1].scatter(raw_df['global_CG_level'], 
                             np.log10(raw_df['CG_sites']), 
                             alpha=0.1, color='black', s=1)
            axes[0, 1].axvline(x=0.5, color='red', linestyle='--')
            axes[0, 1].axhline(y=np.log10(100000), color='red', linestyle='--')
            axes[0, 1].axhline(y=np.log10(2500000), color='red', linestyle='--')
            axes[0, 1].set_title('CG Sites vs Level (Raw)')
            axes[0, 1].set_xlabel('Global CG Level')
            axes[0, 1].set_ylabel('log10(CG Sites)')
            
            # CH plots - raw
            axes[0, 2].hist(raw_df['global_CH_level'], bins=50, color='black')
            axes[0, 2].axvline(x=0.1, color='red', linestyle='--')
            axes[0, 2].set_title('CH Level Distribution (Raw)')
            axes[0, 2].set_xlabel('Global CH Level')
            
            axes[0, 3].scatter(raw_df['global_CH_level'], 
                             np.log10(raw_df['CH_sites']), 
                             alpha=0.1, color='black', s=1)
            axes[0, 3].axvline(x=0.1, color='red', linestyle='--')
            axes[0, 3].axhline(y=np.log10(100000), color='red', linestyle='--')
            axes[0, 3].set_title('CH Sites vs Level (Raw)')
            axes[0, 3].set_xlabel('Global CH Level')
            axes[0, 3].set_ylabel('log10(CH Sites)')
            
            # Filtered plots
            axes[1, 0].hist(filtered_df['global_CG_level'], bins=50, color='black')
            axes[1, 0].set_title('CG Level Distribution (Filtered)')
            axes[1, 0].set_xlabel('Global CG Level')
            
            axes[1, 1].scatter(filtered_df['global_CG_level'], 
                             np.log10(filtered_df['CG_sites']), 
                             alpha=0.1, color='black', s=1)
            axes[1, 1].set_title('CG Sites vs Level (Filtered)')
            axes[1, 1].set_xlabel('Global CG Level')
            axes[1, 1].set_ylabel('log10(CG Sites)')
            
            axes[1, 2].hist(filtered_df['global_CH_level'], bins=50, color='black')
            axes[1, 2].set_title('CH Level Distribution (Filtered)')
            axes[1, 2].set_xlabel('Global CH Level')
            
            axes[1, 3].scatter(filtered_df['global_CH_level'], 
                             np.log10(filtered_df['CH_sites']), 
                             alpha=0.1, color='black', s=1)
            axes[1, 3].set_title('CH Sites vs Level (Filtered)')
            axes[1, 3].set_xlabel('Global CH Level')
            axes[1, 3].set_ylabel('log10(CH Sites)')
            
            plt.tight_layout()
            plt.savefig(self.output_dir / "qc_metrics.pdf", dpi=300)
            plt.close()
            
            logger.info(f"QC plots saved to {self.output_dir / 'qc_metrics.pdf'}")
            
        except ImportError:
            logger.warning("Matplotlib not available, skipping QC plots")
