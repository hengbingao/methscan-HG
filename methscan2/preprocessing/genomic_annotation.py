"""
Genomic Region Annotation Module

Calculate methylation levels across genomic features:
- TSS (Transcription Start Sites)
- Gene body
- Upstream/Downstream regions
- Promoters
- Enhancers
- Custom regions
"""

import numpy as np
import pandas as pd
import pyranges as pr
from pathlib import Path
from typing import Optional, Union, List, Dict, Tuple
from joblib import Parallel, delayed
from tqdm.auto import tqdm
import warnings


class GenomicRegionAnnotator:
    """
    Annotate and quantify methylation across genomic regions
    """
    
    def __init__(
        self,
        gtf_file: Optional[str] = None,
        bed_file: Optional[str] = None,
        genome: str = 'hg38'
    ):
        """
        Initialize with genomic annotation
        
        Parameters:
            gtf_file: GTF/GFF annotation file
            bed_file: BED file with custom regions
            genome: Genome assembly name
        """
        self.genome = genome
        self.genes = None
        self.custom_regions = None
        
        if gtf_file:
            self.load_gtf(gtf_file)
        elif bed_file:
            self.load_bed(bed_file)
    
    def load_gtf(self, gtf_file: str):
        """Load gene annotations from GTF file"""
        print(f"Loading GTF: {gtf_file}")
        
        # Read GTF using pyranges
        self.genes = pr.read_gtf(gtf_file)
        
        # Filter for genes only
        if 'Feature' in self.genes.columns:
            self.genes = self.genes[self.genes.Feature == 'gene']
        
        print(f"Loaded {len(self.genes)} genes")
    
    def load_bed(self, bed_file: str):
        """Load custom regions from BED file"""
        print(f"Loading BED: {bed_file}")
        self.custom_regions = pr.read_bed(bed_file)
        print(f"Loaded {len(self.custom_regions)} regions")
    
    def define_tss_regions(
        self,
        upstream: int = 2000,
        downstream: int = 500
    ) -> pr.PyRanges:
        """
        Define TSS regions (upstream + downstream of TSS)
        
        Parameters:
            upstream: Distance upstream of TSS
            downstream: Distance downstream of TSS
            
        Returns:
            PyRanges object with TSS regions
        """
        if self.genes is None:
            raise ValueError("No gene annotations loaded")
        
        # Convert to DataFrame for easier manipulation
        df = self.genes.df
        
        # Calculate TSS positions based on strand
        tss_regions = []
        
        for _, gene in df.iterrows():
            if gene['Strand'] == '+':
                # Forward strand: TSS is start
                tss = gene['Start']
                start = max(0, tss - upstream)
                end = tss + downstream
            else:
                # Reverse strand: TSS is end
                tss = gene['End']
                start = max(0, tss - downstream)
                end = tss + upstream
            
            tss_regions.append({
                'Chromosome': gene['Chromosome'],
                'Start': start,
                'End': end,
                'Strand': gene['Strand'],
                'gene_id': gene.get('gene_id', ''),
                'gene_name': gene.get('gene_name', ''),
                'region_type': 'TSS'
            })
        
        return pr.PyRanges(pd.DataFrame(tss_regions))
    
    def define_gene_body_regions(self) -> pr.PyRanges:
        """Define gene body regions (from TSS to TES)"""
        if self.genes is None:
            raise ValueError("No gene annotations loaded")
        
        df = self.genes.df.copy()
        df['region_type'] = 'gene_body'
        
        return pr.PyRanges(df)
    
    def define_upstream_downstream_regions(
        self,
        upstream: int = 5000,
        downstream: int = 5000
    ) -> Dict[str, pr.PyRanges]:
        """
        Define upstream and downstream regions
        
        Returns:
            Dictionary with 'upstream' and 'downstream' PyRanges
        """
        if self.genes is None:
            raise ValueError("No gene annotations loaded")
        
        df = self.genes.df
        
        upstream_regions = []
        downstream_regions = []
        
        for _, gene in df.iterrows():
            if gene['Strand'] == '+':
                # Upstream is before start
                up_start = max(0, gene['Start'] - upstream)
                up_end = gene['Start']
                # Downstream is after end
                down_start = gene['End']
                down_end = gene['End'] + downstream
            else:
                # Reverse strand
                up_start = gene['End']
                up_end = gene['End'] + upstream
                down_start = max(0, gene['Start'] - downstream)
                down_end = gene['Start']
            
            upstream_regions.append({
                'Chromosome': gene['Chromosome'],
                'Start': up_start,
                'End': up_end,
                'Strand': gene['Strand'],
                'gene_id': gene.get('gene_id', ''),
                'gene_name': gene.get('gene_name', ''),
                'region_type': 'upstream'
            })
            
            downstream_regions.append({
                'Chromosome': gene['Chromosome'],
                'Start': down_start,
                'End': down_end,
                'Strand': gene['Strand'],
                'gene_id': gene.get('gene_id', ''),
                'gene_name': gene.get('gene_name', ''),
                'region_type': 'downstream'
            })
        
        return {
            'upstream': pr.PyRanges(pd.DataFrame(upstream_regions)),
            'downstream': pr.PyRanges(pd.DataFrame(downstream_regions))
        }
    
    def define_promoter_regions(
        self,
        upstream: int = 2000,
        downstream: int = 500
    ) -> pr.PyRanges:
        """
        Define promoter regions (similar to TSS but commonly used term)
        
        This is an alias for define_tss_regions with standard promoter distances
        """
        return self.define_tss_regions(upstream, downstream)


def calculate_region_methylation(
    allc_file: str,
    regions: pr.PyRanges,
    context: str = 'CG',
    min_coverage: int = 1,
    region_name: str = 'region'
) -> pd.DataFrame:
    """
    Calculate methylation levels within specified genomic regions
    
    Parameters:
        allc_file: Path to ALLC file
        regions: PyRanges object with genomic regions
        context: Methylation context ('CG', 'CH', 'CHG', 'CHH', 'all')
        min_coverage: Minimum coverage threshold
        region_name: Name for this region type
        
    Returns:
        DataFrame with region methylation statistics
    """
    import gzip
    
    # Read ALLC file
    if allc_file.endswith('.gz'):
        with gzip.open(allc_file, 'rt') as f:
            allc_df = pd.read_csv(
                f,
                sep='\t',
                header=None,
                names=['chr', 'pos', 'strand', 'mc_context', 'mc', 'cov', 'meth_level']
            )
    else:
        allc_df = pd.read_csv(
            allc_file,
            sep='\t',
            header=None,
            names=['chr', 'pos', 'strand', 'mc_context', 'mc', 'cov', 'meth_level']
        )
    
    # Filter by context
    if context == 'CG':
        allc_df = allc_df[allc_df['mc_context'].str.startswith('CG')]
    elif context == 'CH':
        allc_df = allc_df[allc_df['mc_context'].str.match(r'^C[ATC]')]
    elif context == 'CHG':
        allc_df = allc_df[allc_df['mc_context'].str.startswith('CHG')]
    elif context == 'CHH':
        allc_df = allc_df[allc_df['mc_context'].str.startswith('CHH')]
    # 'all' means use all contexts
    
    # Filter by coverage
    allc_df = allc_df[allc_df['cov'] >= min_coverage]
    
    # Convert ALLC to PyRanges
    allc_pr = pr.PyRanges(
        chromosomes=allc_df['chr'],
        starts=allc_df['pos'],
        ends=allc_df['pos'] + 1,
        mc=allc_df['mc'],
        cov=allc_df['cov']
    )
    
    # Find overlaps with regions
    overlaps = regions.join(allc_pr)
    
    # Calculate methylation for each region
    if len(overlaps) == 0:
        return pd.DataFrame()
    
    overlap_df = overlaps.df
    
    # Group by region and calculate statistics
    results = []
    
    for region_id in overlap_df['gene_id'].unique() if 'gene_id' in overlap_df.columns else range(len(regions)):
        if 'gene_id' in overlap_df.columns:
            region_sites = overlap_df[overlap_df['gene_id'] == region_id]
            gene_name = region_sites['gene_name'].iloc[0] if 'gene_name' in region_sites.columns else region_id
        else:
            region_sites = overlap_df
            gene_name = f'region_{region_id}'
        
        if len(region_sites) == 0:
            continue
        
        total_mc = region_sites['mc'].sum()
        total_cov = region_sites['cov'].sum()
        
        meth_level = total_mc / total_cov if total_cov > 0 else np.nan
        n_sites = len(region_sites)
        
        results.append({
            'gene_id': region_id if 'gene_id' in overlap_df.columns else gene_name,
            'gene_name': gene_name,
            'region_type': region_name,
            'context': context,
            'methylation_level': meth_level,
            'n_sites': n_sites,
            'total_mc': total_mc,
            'total_cov': total_cov
        })
    
    return pd.DataFrame(results)


def batch_calculate_region_methylation(
    allc_files: List[str],
    regions_dict: Dict[str, pr.PyRanges],
    contexts: List[str] = ['CG', 'CH'],
    min_coverage: int = 1,
    n_jobs: int = -1,
    verbose: bool = True
) -> Dict[str, pd.DataFrame]:
    """
    Calculate region methylation for multiple cells and contexts
    
    Parameters:
        allc_files: List of ALLC files
        regions_dict: Dictionary of region types to PyRanges objects
            e.g., {'tss': tss_regions, 'gene_body': gb_regions}
        contexts: List of methylation contexts
        min_coverage: Minimum coverage threshold
        n_jobs: Number of parallel jobs
        verbose: Show progress
        
    Returns:
        Dictionary with results for each region type and context
    """
    if verbose:
        print(f"Calculating region methylation for {len(allc_files)} cells")
        print(f"Regions: {list(regions_dict.keys())}")
        print(f"Contexts: {contexts}")
    
    all_results = {f"{region}_{ctx}": [] for region in regions_dict for ctx in contexts}
    
    # Process each cell
    for allc_file in (tqdm(allc_files, desc="Processing cells") if verbose else allc_files):
        cell_name = Path(allc_file).stem
        
        # Process each region type and context
        for region_name, regions in regions_dict.items():
            for context in contexts:
                result_df = calculate_region_methylation(
                    allc_file=allc_file,
                    regions=regions,
                    context=context,
                    min_coverage=min_coverage,
                    region_name=region_name
                )
                
                if len(result_df) > 0:
                    result_df['cell'] = cell_name
                    all_results[f"{region_name}_{context}"].append(result_df)
    
    # Concatenate results
    final_results = {}
    for key, dfs in all_results.items():
        if len(dfs) > 0:
            final_results[key] = pd.concat(dfs, ignore_index=True)
    
    return final_results
