# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'CelLink'
copyright = '2024, Xin Luo'
author = 'Xin Luo'
release = '0.1.3'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

templates_path = ['_templates']
exclude_patterns = []

language = 'English'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']

extensions = [
    'nbsphinx',
    'sphinx.ext.mathjax',  # For rendering equations
    'sphinx.ext.viewcode', # Optional, for linking code
]

# Exclude build directory
exclude_patterns = ['_build', '**.ipynb_checkpoints']
source_suffix = ['.rst', '.ipynb']

# Set notebooks to render on a specific kernel (if needed)
nbsphinx_kernel_name = 'python3'

