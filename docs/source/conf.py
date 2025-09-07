# Configuration file for the Sphinx documentation builder.

# -- Project information

project = 'CelLink'
copyright = '2024, Xin Luo'
author = 'Xin Luo'

release = '1.0'
version = '0.1.4'

# -- General configuration

extensions = [
    'nbsphinx',
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
]

# nbsphinx configuration: do not execute notebooks during documentation build on ReadTheDocs
nbsphinx_execute = 'never'
# Alternatively, to avoid execution but still show outputs, use 'never'. If you want execution set to 'auto' or 'always'.

nbsphinx_allow_errors = False
nbsphinx_timeout = 120

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'

# -- Options for EPUB output
epub_show_urls = 'footnote'

def setup(app):
    app.add_css_file('my_theme.css')

html_static_path = ['_static']
