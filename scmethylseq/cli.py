"""
Command-line interface for scMethylSeq
"""

import click
import logging
from pathlib import Path
import yaml
from typing import Optional

from .preprocessing import (
    AllcProcessor,
    CovFileGenerator,
    GlobalMethylationCalculator,
    QualityControl
)
from .matrix_builder import (
    GeneBodyMethylation,
    MethylationMatrixBuilder
)
from .clustering import (
    PCAAnalyzer,
    UMAPAnalyzer,
    LeidenClustering
)
from .adata_builder import AnnDataBuilder

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


@click.group()
@click.version_option()
def main():
    """
    scMethylSeq - Single-cell DNA Methylation Analysis Toolkit
    
    A comprehensive toolkit for analyzing single-cell DNA methylation data
    from allc files to AnnData objects.
    """
    pass


@main.command()
@click.option('--allc-dir', required=True, type=click.Path(exists=True),
              help='Directory containing allc.gz files')
@click.option('--output-dir', required=True, type=click.Path(),
              help='Output directory')
@click.option('--chrom-sizes', required=True, type=click.Path(exists=True),
              help='Chromosome sizes file')
@click.option('--n-jobs', default=10, type=int,
              help='Number of parallel jobs')
def process_allc(allc_dir, output_dir, chrom_sizes, n_jobs):
    """
    Process allc files to extract CG and CH contexts
    """
    logger.info("Starting allc processing...")
    
    processor = AllcProcessor(
        allc_dir=allc_dir,
        output_dir=output_dir,
        chrom_sizes=chrom_sizes,
        n_jobs=n_jobs
    )
    
    processor.process_all()
    logger.info("Allc processing complete!")


@main.command()
@click.option('--allc-dir', required=True, type=click.Path(exists=True),
              help='Directory containing allc.gz files')
@click.option('--output-dir', required=True, type=click.Path(),
              help='Output directory for cov files')
@click.option('--context', type=click.Choice(['CG', 'CH', 'all']), default='all',
              help='Methylation context to convert')
@click.option('--n-jobs', default=20, type=int,
              help='Number of parallel jobs')
def convert_to_cov(allc_dir, output_dir, context, n_jobs):
    """
    Convert allc files to cov format
    """
    logger.info("Starting allc to cov conversion...")
    
    converter = CovFileGenerator(
        allc_dir=allc_dir,
        output_dir=output_dir,
        n_jobs=n_jobs
    )
    
    if context == 'all':
        converter.convert_all(pattern="*.allc.gz")
    elif context == 'CG':
        converter.convert_all(pattern="*_CG.allc.gz")
    else:
        converter.convert_all(pattern="*_CH.allc.gz")
    
    logger.info("Conversion complete!")


@main.command()
@click.option('--cov-dir', required=True, type=click.Path(exists=True),
              help='Directory containing cov files')
@click.option('--output-file', required=True, type=click.Path(),
              help='Output file path')
@click.option('--context', type=click.Choice(['CG', 'CH']), required=True,
              help='Methylation context')
@click.option('--n-jobs', default=20, type=int,
              help='Number of parallel jobs')
def calc_global_meth(cov_dir, output_file, context, n_jobs):
    """
    Calculate global methylation levels
    """
    logger.info(f"Calculating global {context} methylation...")
    
    calculator = GlobalMethylationCalculator(
        cov_dir=cov_dir,
        output_file=output_file,
        context=context,
        n_jobs=n_jobs
    )
    
    df = calculator.calculate_all()
    logger.info(f"Results saved to {output_file}")


@main.command()
@click.option('--cg-file', required=True, type=click.Path(exists=True),
              help='Global CG metrics file')
@click.option('--ch-file', required=True, type=click.Path(exists=True),
              help='Global CH metrics file')
@click.option('--output-dir', required=True, type=click.Path(),
              help='Output directory')
@click.option('--cg-min', default=0.5, type=float,
              help='Minimum CG methylation level')
@click.option('--cg-sites-min', default=100000, type=int,
              help='Minimum CG sites')
@click.option('--cg-sites-max', default=2500000, type=int,
              help='Maximum CG sites')
@click.option('--ch-max', default=0.1, type=float,
              help='Maximum CH methylation level')
@click.option('--ch-sites-min', default=100000, type=int,
              help='Minimum CH sites')
def filter_cells(cg_file, ch_file, output_dir, 
                cg_min, cg_sites_min, cg_sites_max, 
                ch_max, ch_sites_min):
    """
    Filter cells based on QC metrics
    """
    logger.info("Filtering cells...")
    
    qc = QualityControl(
        cg_file=cg_file,
        ch_file=ch_file,
        output_dir=output_dir
    )
    
    filtered_df = qc.filter_cells(
        cg_min=cg_min,
        cg_sites_min=cg_sites_min,
        cg_sites_max=cg_sites_max,
        ch_max=ch_max,
        ch_sites_min=ch_sites_min
    )
    
    logger.info(f"Filtered {len(filtered_df)} cells")


@main.command()
@click.option('--cov-dir', required=True, type=click.Path(exists=True),
              help='Directory containing cov files')
@click.option('--gene-bed', required=True, type=click.Path(exists=True),
              help='Gene body BED file')
@click.option('--output-dir', required=True, type=click.Path(),
              help='Output directory')
@click.option('--context', type=click.Choice(['CG', 'CH']), required=True,
              help='Methylation context')
@click.option('--n-jobs', default=20, type=int,
              help='Number of parallel jobs')
def calc_gene_meth(cov_dir, gene_bed, output_dir, context, n_jobs):
    """
    Calculate gene body methylation levels
    """
    logger.info(f"Calculating {context} gene body methylation...")
    
    gene_meth = GeneBodyMethylation(
        cov_dir=cov_dir,
        gene_bed=gene_bed,
        output_dir=output_dir,
        context=context,
        n_jobs=n_jobs
    )
    
    output_files = gene_meth.process_all_cells()
    logger.info(f"Processed {len(output_files)} cells")


@main.command()
@click.option('--gene-meth-dir', required=True, type=click.Path(exists=True),
              help='Directory containing gene methylation files')
@click.option('--gene-names-file', required=True, type=click.Path(exists=True),
              help='File with gene names')
@click.option('--output-file', required=True, type=click.Path(),
              help='Output matrix file')
@click.option('--context', type=click.Choice(['CG', 'CH']), required=True,
              help='Methylation context')
@click.option('--cleanup/--no-cleanup', default=False,
              help='Remove intermediate files after building matrix')
def build_matrix(gene_meth_dir, gene_names_file, output_file, context, cleanup):
    """
    Build methylation matrix from gene methylation files
    """
    logger.info(f"Building {context} methylation matrix...")
    
    builder = MethylationMatrixBuilder(
        gene_meth_dir=gene_meth_dir,
        gene_names_file=gene_names_file,
        output_file=output_file,
        context=context
    )
    
    matrix = builder.build_matrix()
    
    if cleanup:
        logger.info("Cleaning up intermediate files...")
        builder.clean_intermediate_files()
    
    logger.info(f"Matrix saved to {output_file}")


@main.command()
@click.option('--cg-matrix', required=True, type=click.Path(exists=True),
              help='CG methylation matrix file')
@click.option('--ch-matrix', required=True, type=click.Path(exists=True),
              help='CH methylation matrix file')
@click.option('--metadata', required=True, type=click.Path(exists=True),
              help='Cell metadata file (CSV)')
@click.option('--output', required=True, type=click.Path(),
              help='Output h5ad file')
def build_adata(cg_matrix, ch_matrix, metadata, output):
    """
    Build AnnData object from methylation matrices
    """
    logger.info("Building AnnData object...")
    
    builder = AnnDataBuilder(
        cg_matrix_file=cg_matrix,
        ch_matrix_file=ch_matrix,
        metadata_file=metadata,
        output_file=output
    )
    
    adata = builder.build_and_save()
    logger.info(f"AnnData saved to {output}")


@main.command()
@click.option('--config', required=True, type=click.Path(exists=True),
              help='Configuration YAML file')
def run_pipeline(config):
    """
    Run complete pipeline from configuration file
    """
    logger.info(f"Loading configuration from {config}")
    
    with open(config, 'r') as f:
        cfg = yaml.safe_load(f)
    
    logger.info("Starting complete pipeline...")
    
    # Step 1: Process allc files if specified
    if cfg.get('run_allc_processing', False):
        logger.info("=" * 60)
        logger.info("STEP 1: Processing allc files")
        logger.info("=" * 60)
        
        processor = AllcProcessor(
            allc_dir=cfg['allc_dir'],
            output_dir=cfg['processed_allc_dir'],
            chrom_sizes=cfg['chrom_sizes'],
            n_jobs=cfg.get('n_jobs', 10)
        )
        processor.process_all()
    
    # Step 2: Convert to cov
    if cfg.get('run_cov_conversion', False):
        logger.info("=" * 60)
        logger.info("STEP 2: Converting to cov format")
        logger.info("=" * 60)
        
        for context in ['CG', 'CH']:
            converter = CovFileGenerator(
                allc_dir=cfg['processed_allc_dir'],
                output_dir=cfg[f'{context.lower()}_cov_dir'],
                n_jobs=cfg.get('n_jobs', 20)
            )
            converter.convert_all(pattern=f"*_{context}.allc.gz")
    
    # Step 3: Calculate global methylation
    if cfg.get('run_global_meth', False):
        logger.info("=" * 60)
        logger.info("STEP 3: Calculating global methylation")
        logger.info("=" * 60)
        
        for context in ['CG', 'CH']:
            calculator = GlobalMethylationCalculator(
                cov_dir=cfg[f'{context.lower()}_cov_dir'],
                output_file=cfg[f'global_{context.lower()}_file'],
                context=context,
                n_jobs=cfg.get('n_jobs', 20)
            )
            calculator.calculate_all()
    
    # Step 4: QC and filter cells
    if cfg.get('run_qc', False):
        logger.info("=" * 60)
        logger.info("STEP 4: Quality control and filtering")
        logger.info("=" * 60)
        
        qc = QualityControl(
            cg_file=cfg['global_cg_file'],
            ch_file=cfg['global_ch_file'],
            output_dir=cfg['qc_output_dir']
        )
        
        qc_params = cfg.get('qc_params', {})
        filtered_df = qc.filter_cells(**qc_params)
    
    # Step 5: Calculate gene body methylation
    if cfg.get('run_gene_meth', False):
        logger.info("=" * 60)
        logger.info("STEP 5: Calculating gene body methylation")
        logger.info("=" * 60)
        
        for context in ['CG', 'CH']:
            gene_meth = GeneBodyMethylation(
                cov_dir=cfg[f'{context.lower()}_cov_dir'],
                gene_bed=cfg['gene_bed'],
                output_dir=cfg[f'{context.lower()}_gene_meth_dir'],
                context=context,
                n_jobs=cfg.get('n_jobs', 20)
            )
            gene_meth.process_all_cells()
    
    # Step 6: Build matrices
    if cfg.get('run_matrix_building', False):
        logger.info("=" * 60)
        logger.info("STEP 6: Building methylation matrices")
        logger.info("=" * 60)
        
        for context in ['CG', 'CH']:
            builder = MethylationMatrixBuilder(
                gene_meth_dir=cfg[f'{context.lower()}_gene_meth_dir'],
                gene_names_file=cfg['gene_names_file'],
                output_file=cfg[f'{context.lower()}_matrix_file'],
                context=context
            )
            builder.build_matrix()
            
            if cfg.get('cleanup_intermediate', False):
                builder.clean_intermediate_files()
    
    # Step 7: Build AnnData
    if cfg.get('run_adata_building', False):
        logger.info("=" * 60)
        logger.info("STEP 7: Building AnnData object")
        logger.info("=" * 60)
        
        # Metadata should include UMAP and clustering results
        # These need to be generated separately or provided
        
        adata_builder = AnnDataBuilder(
            cg_matrix_file=cfg['cg_matrix_file'],
            ch_matrix_file=cfg['ch_matrix_file'],
            metadata_file=cfg['metadata_file'],
            output_file=cfg['output_h5ad']
        )
        adata_builder.build_and_save()
    
    logger.info("=" * 60)
    logger.info("Pipeline complete!")
    logger.info("=" * 60)


if __name__ == '__main__':
    main()
