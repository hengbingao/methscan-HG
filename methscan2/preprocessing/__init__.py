"""
MethSCAn2 Preprocessing Module
"""

from .io import read_cov_files, read_allc_files, read_bedgraph
from .qc import calculate_qc_metrics, plot_qc_metrics
from .filter import filter_cells, filter_sites
from .smooth import smooth_methylation
from .allc_convert import (
    allc_to_cov,
    extract_context_from_allc,
    calculate_global_methylation,
    batch_process_allc,
    merge_global_stats
)
from .genomic_annotation import (
    GenomicRegionAnnotator,
    calculate_region_methylation,
    batch_calculate_region_methylation
)

__all__ = [
    'read_cov_files',
    'read_allc_files', 
    'read_bedgraph',
    'calculate_qc_metrics',
    'plot_qc_metrics',
    'filter_cells',
    'filter_sites',
    'smooth_methylation',
    'allc_to_cov',
    'extract_context_from_allc',
    'calculate_global_methylation',
    'batch_process_allc',
    'merge_global_stats',
    'GenomicRegionAnnotator',
    'calculate_region_methylation',
    'batch_calculate_region_methylation',
]
