# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import sys
from pathlib import Path

sys.path.insert(0, str(Path('..', '..', 'src').resolve()))
#sys.path.insert(0, os.path.abspath('.'))

project = 'The PPOLs Model'
copyright = '2025, S. McCloat'
author = 'S. McCloat'
release = '1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
'sphinx.ext.autodoc',
#'sphinx.ext.autosummary',
'numpydoc']

templates_path = ['_templates']
exclude_patterns = []

autodoc_mock_imports = [
	"astropy",
	"matplotlib"
	]



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
