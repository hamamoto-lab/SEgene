"""
SEgene_analyzer_erna - eRNAbase Analysis Tool

This package provides tools for analyzing genomic regions against the eRNAbase database.
It includes functionality for identifying overlapping eRNAs, analyzing tissue specificity,
and generating comprehensive reports.
"""

from .erna_analyzer import ERNARegionAnalyzer
from .sedb_analyzer import SEdbRegionAnalyzer

__version__ = '1.0.0'
__all__ = ['ERNARegionAnalyzer', 'SEdbRegionAnalyzer']

# Export CLI entry point
from . import cli