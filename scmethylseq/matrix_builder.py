"""
Matrix builder module
Construct methylation matrices for gene bodies and genomic features
"""

import os
import subprocess
from pathlib import Path
from typing import Optional, List, Dict
from concurrent.futures import ProcessPoolExecutor, as_completed
import logging
import pandas as pd
import numpy as np
from tqdm import tqdm
import tempfile

logger = logging.getLogger(__name__)


class GeneBodyMethylation:
    """
    Calculate methylation levels for gene bodies
    """
    
    def __init__(self,
                 cov_dir: str,
                 gene_bed: str,
                 output_dir: str,
                 context: str = "CG",
                 n_jobs: int = 20):
        """
        Initialize GeneBodyMethylation
        
        Parameters:
        -----------
        cov_dir : str
            Directory containing cov files
        gene_bed : str
            BED file with gene coordinates (e.g., genebody_extend2kb.bed)
        output_dir : str
            Output directory
        context : str
            Methylation context ("CG" or "CH")
        n_jobs : int
            Number of parallel processes
        """
        self.cov_dir = Path(cov_dir)
        self.gene_bed = gene_bed
        self.output_dir = Path(output_dir)
        self.context = context
        self.n_jobs = n_jobs
        
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def _process_single_cell(self, cov_file: Path) -> Tuple[bool, str, Path]:
        """
        Calculate gene body methylation for a single cell
        
        Returns:
        --------
        Tuple[bool, str, Path]
            Success status, message, output file path
        """
        try:
            cell_name = cov_file.stem.replace('.cov', '')
            
            # Create temporary directory for intermediate files
            with tempfile.TemporaryDirectory() as tmpdir:
                tmpdir = Path(tmpdir)
                
                # Convert cov to BED format and sort
                bed_file = tmpdir / f"{cell_name}.bed"
                with open(cov_file, 'r') as fin, open(bed_file, 'w') as fout:
                    for line in fin:
                        parts = line.strip().split('\t')
                        if len(parts) >= 6:
                            chrom, start, end, meth_pct = parts[0], parts[1], parts[2], parts[3]
                            # Format: chr start end meth_pct count=1
                            fout.write(f"{chrom}\t{start}\t{int(start)+1}\t{meth_pct}\t1\n")
                
                # Sort the BED file
                sorted_bed = tmpdir / f"{cell_name}.sorted.bed"
                sort_cmd = f"bedtools sort -i {bed_file} > {sorted_bed}"
                subprocess.run(sort_cmd, shell=True, check=True,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                
                # Calculate methylation sum per gene
                sum_file = tmpdir / f"{cell_name}.sum"
                sum_cmd = f"bedtools map -a {self.gene_bed} -b {sorted_bed} -c 4 -o sum > {sum_file}"
                subprocess.run(sum_cmd, shell=True, check=True,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                
                # Calculate site count per gene
                count_file = tmpdir / f"{cell_name}.count"
                count_cmd = f"bedtools map -a {self.gene_bed} -b {sorted_bed} -c 5 -o sum > {count_file}"
                subprocess.run(count_cmd, shell=True, check=True,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                
                # Read results and calculate methylation levels
                sum_df = pd.read_csv(sum_file, sep='\t', header=None)
                count_df = pd.read_csv(count_file, sep='\t', header=None)
                
                # Get the last column (methylation sum and count)
                meth_sum = sum_df.iloc[:, -1]
                site_count = count_df.iloc[:, -1]
                
                # Calculate methylation level
                # Replace '.' with NaN
                meth_sum = pd.to_numeric(meth_sum, errors='coerce')
                site_count = pd.to_numeric(site_count, errors='coerce')
                
                # Calculate level (sum/count/100 to get 0-1 scale)
                meth_level = meth_sum / site_count / 100
                meth_level = meth_level.fillna(np.nan)
                
                # Save output
                output_file = self.output_dir / f"{cell_name}.gene_m{self.context}"
                with open(output_file, 'w') as fout:
                    fout.write(f"{cell_name}\n")
                    for val in meth_level:
                        if pd.isna(val):
                            fout.write("NA\n")
                        else:
                            fout.write(f"{val:.6f}\n")
                
                return True, f"Processed {cell_name}", output_file
                
        except Exception as e:
            logger.error(f"Error processing {cov_file.name}: {str(e)}")
            return False, f"Error: {str(e)}", None
    
    def process_all_cells(self) -> List[Path]:
        """
        Process all cells in parallel
        
        Returns:
        --------
        List[Path]
            List of output file paths
        """
        cov_files = list(self.cov_dir.glob("*.cov"))
        logger.info(f"Processing {len(cov_files)} cells for {self.context} gene body methylation")
        
        output_files = []
        with ProcessPoolExecutor(max_workers=self.n_jobs) as executor:
            futures = {executor.submit(self._process_single_cell, f): f 
                      for f in cov_files}
            
            with tqdm(total=len(cov_files), desc=f"Processing {self.context}") as pbar:
                for future in as_completed(futures):
                    success, message, output_file = future.result()
                    if success:
                        output_files.append(output_file)
                    pbar.update(1)
        
        logger.info(f"Successfully processed {len(output_files)} cells")
        return output_files


class MethylationMatrixBuilder:
    """
    Build cell x gene methylation matrix
    """
    
    def __init__(self,
                 gene_meth_dir: str,
                 gene_names_file: str,
                 output_file: str,
                 context: str = "CG"):
        """
        Initialize MethylationMatrixBuilder
        
        Parameters:
        -----------
        gene_meth_dir : str
            Directory containing gene methylation files
        gene_names_file : str
            File containing gene names (one per line)
        output_file : str
            Output matrix file
        context : str
            Methylation context ("CG" or "CH")
        """
        self.gene_meth_dir = Path(gene_meth_dir)
        self.gene_names_file = gene_names_file
        self.output_file = output_file
        self.context = context
    
    def build_matrix(self) -> pd.DataFrame:
        """
        Build methylation matrix by concatenating all cell files
        
        Returns:
        --------
        pd.DataFrame
            Cell x gene methylation matrix
        """
        logger.info(f"Building {self.context} methylation matrix")
        
        # Get all gene methylation files
        pattern = f"*.gene_m{self.context}"
        meth_files = sorted(list(self.gene_meth_dir.glob(pattern)))
        
        if len(meth_files) == 0:
            raise ValueError(f"No files matching pattern {pattern} found in {self.gene_meth_dir}")
        
        logger.info(f"Found {len(meth_files)} cell files")
        
        # Read gene names
        with open(self.gene_names_file, 'r') as f:
            gene_names = [line.strip() for line in f if line.strip()]
        
        logger.info(f"Loaded {len(gene_names)} gene names")
        
        # Paste all files together using system paste command (more efficient)
        logger.info("Concatenating cell files...")
        
        # Create output with gene names
        temp_output = self.output_file + ".tmp"
        
        # Use paste command to combine all files
        file_list = " ".join([str(f) for f in meth_files])
        paste_cmd = f"paste {self.gene_names_file} {file_list} > {temp_output}"
        
        try:
            subprocess.run(paste_cmd, shell=True, check=True,
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            
            # Move temp file to final output
            subprocess.run(f"mv {temp_output} {self.output_file}", shell=True, check=True)
            
            logger.info(f"Matrix saved to {self.output_file}")
            
            # Also load and return as DataFrame
            logger.info("Loading matrix into memory...")
            df = pd.read_csv(self.output_file, sep='\t', header=None)
            
            # Set gene names as first column
            df.iloc[:, 0] = gene_names
            
            return df
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Error building matrix: {str(e)}")
            raise
    
    def clean_intermediate_files(self) -> None:
        """
        Remove intermediate gene methylation files
        """
        pattern = f"*.gene_m{self.context}"
        files_to_remove = list(self.gene_meth_dir.glob(pattern))
        
        logger.info(f"Removing {len(files_to_remove)} intermediate files")
        for f in tqdm(files_to_remove, desc="Cleaning up"):
            f.unlink()


class VMRMethylationExtractor:
    """
    Extract methylation levels for Variable Methylated Regions (VMRs)
    """
    
    def __init__(self,
                 cov_dir: str,
                 vmr_bed: str,
                 output_dir: str,
                 context: str = "CG",
                 n_jobs: int = 20):
        """
        Initialize VMRMethylationExtractor
        
        Parameters:
        -----------
        cov_dir : str
            Directory containing cov files
        vmr_bed : str
            BED file with VMR coordinates
        output_dir : str
            Output directory
        context : str
            Methylation context
        n_jobs : int
            Number of parallel processes
        """
        self.cov_dir = Path(cov_dir)
        self.vmr_bed = vmr_bed
        self.output_dir = Path(output_dir)
        self.context = context
        self.n_jobs = n_jobs
        
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def extract_vmr_methylation(self) -> pd.DataFrame:
        """
        Extract VMR methylation for all cells
        
        Returns:
        --------
        pd.DataFrame
            Cell x VMR methylation matrix
        """
        # Similar implementation to GeneBodyMethylation
        # but for VMR regions
        logger.info("VMR extraction not yet implemented")
        raise NotImplementedError("VMR extraction will be implemented in future version")


class TSSMethylationExtractor:
    """
    Extract methylation levels around Transcription Start Sites (TSS)
    """
    
    def __init__(self,
                 cov_dir: str,
                 tss_bed: str,
                 output_dir: str,
                 window_size: int = 2000,
                 context: str = "CG",
                 n_jobs: int = 20):
        """
        Initialize TSSMethylationExtractor
        
        Parameters:
        -----------
        cov_dir : str
            Directory containing cov files
        tss_bed : str
            BED file with TSS coordinates
        output_dir : str
            Output directory
        window_size : int
            Window size around TSS (bp)
        context : str
            Methylation context
        n_jobs : int
            Number of parallel processes
        """
        self.cov_dir = Path(cov_dir)
        self.tss_bed = tss_bed
        self.output_dir = Path(output_dir)
        self.window_size = window_size
        self.context = context
        self.n_jobs = n_jobs
        
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def extract_tss_methylation(self) -> pd.DataFrame:
        """
        Extract TSS methylation for all cells
        
        Returns:
        --------
        pd.DataFrame
            Cell x TSS methylation matrix
        """
        logger.info("TSS extraction not yet implemented")
        raise NotImplementedError("TSS extraction will be implemented in future version")
