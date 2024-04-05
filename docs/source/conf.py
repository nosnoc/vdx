# Configuration file for the Sphinx documentation builder.

# -- Project information

project = 'VDX'
copyright = '2024, Anton Pozharskiy, Armin Nurkanovic, Jonathan Frey, Moritz Diehl'
author = 'Anton Pozharskiy'

release = '0.1'
version = '0.1.0'

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinxcontrib.matlab',
]

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

matlab_src_dir = '.'
matlab_short_links = True
matlab_auto_link = "basic"
