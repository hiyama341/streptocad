"""Top-level package for StreptoCAD."""

__author__ = """Lucas Levassor"""
__email__ = "luclev@biosustain.dtu.dk"
__version__ = "0.1.0"

from streptocad.biopython_compat import patch_dseqrecord_linear_argument

patch_dseqrecord_linear_argument()
