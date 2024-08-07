# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys
sys.path.insert(0,os.path.abspath("../../"))
import chemdiff as cd

project = 'newcandy'
copyright = '2024, Eric Van Clepper'
author = 'Eric Van Clepper'
release = cd.__version__

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.viewcode',      # add links to code in documentation
    'sphinx.ext.autodoc',       # create autodocs for api from docstrings
    'sphinx.ext.napoleon',      # read google style docstrings
    ]

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_theme_path = ['_themes']
html_static_path = ['_static']
