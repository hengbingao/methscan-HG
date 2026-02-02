"""
MethSCAn2 Core Module
"""

from .methylation_data import MethylationData, concat
from .multi_context_data import (
    MultiContextMethylationData,
    create_multi_context_data_from_allc
)

__all__ = [
    'MethylationData',
    'concat',
    'MultiContextMethylationData',
    'create_multi_context_data_from_allc'
]
