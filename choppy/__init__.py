"""
Choppy: Functionalities for chopping up DNA for TAR cloning.

This package provides tools to chop DNA sequences into smaller fragments
for in vivo DNA assembly by homologous recombination in yeast S. cerevisiae.
"""

from .fragment_annotator import annotate_fragments, extract_no_homology_regions, fragment_from_file
from .homology_finder import annotate_homology, file_process_homology, store_trie_from_file

__version__ = "0.1.0"

__all__ = [
	"annotate_fragments",
	"extract_no_homology_regions",
	"fragment_from_file",
	"annotate_homology",
	"file_process_homology",
	"store_trie_from_file",
]
