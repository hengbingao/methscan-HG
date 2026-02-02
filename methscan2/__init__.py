"""
MethSCAn2 - Single-cell DNA Methylation Analysis Toolkit

A comprehensive Python package for analyzing single-cell DNA methylation data.
"""

__version__ = "0.2.0"
__author__ = "MethSCAn2 Development Team"

# Import main modules
from . import preprocessing as pp
from . import tools as tl
from . import plotting as pl

# Import data readers
from .preprocessing.io import (
    read_cov_files,
    read_allc_files,
    read_bedgraph,
)

# Import core data structures
from .core.methylation_data import MethylationData
from .core.multi_context_data import (
    MultiContextMethylationData,
    create_multi_context_data_from_allc
)

# Set up logging
import logging
logging.basicConfig(
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    level=logging.INFO
)

__all__ = [
    'pp',
    'tl', 
    'pl',
    'read_cov_files',
    'read_allc_files',
    'read_bedgraph',
    'MethylationData',
    'MultiContextMethylationData',
    'create_multi_context_data_from_allc',
]
