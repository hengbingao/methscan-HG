"""
scMethylSeq - Single-cell DNA Methylation Sequencing Analysis Toolkit

A comprehensive Python toolkit for analyzing single-cell DNA methylation data,
from allc files to AnnData objects with clustering and visualization.
"""

__version__ = "0.1.0"
__author__ = "Your Name"

from .preprocessing import (
    AllcProcessor,
    CovFileGenerator,
    GlobalMethylationCalculator,
    QualityControl,
)
from .matrix_builder import (
    MethylationMatrixBuilder,
    GeneBodyMethylation,
)
from .clustering import (
    PCAAnalyzer,
    UMAPAnalyzer,
    LeidenClustering,
)
from .adata_builder import (
    AnnDataBuilder,
)

__all__ = [
    "AllcProcessor",
    "CovFileGenerator",
    "GlobalMethylationCalculator",
    "QualityControl",
    "MethylationMatrixBuilder",
    "GeneBodyMethylation",
    "PCAAnalyzer",
    "UMAPAnalyzer",
    "LeidenClustering",
    "AnnDataBuilder",
]
