# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Pharmacophore-Toolkit'
copyright = '2026, Tony E. Lin'
author = 'Tony E. Lin'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'myst_nb',
    'sphinx_book_theme',
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'sphinx.ext.intersphinx',
    'sphinx.ext.coverage',
    'IPython.sphinxext.ipython_console_highlighting',
    'sphinx_design',
]

nb_execution_mode = "off"
nb_execution_allow_errors = False
nb_execution_raise_on_error = True
nb_execution_timeout = 480

myst_enable_extensions = [
    "colon_fence",
    "html_image",
    "html_admonition",
]

templates_path = ['_templates']
exclude_patterns = []

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "myst-nb",
    ".ipynb": "myst-nb",
    ".myst": "myst-nb",
}

html_theme_options = {
    "repository_url": "https://github.com/tlint101/Pharmacophore-Toolkit",
    "show_toc_level": 3,
    "path_to_docs": "docs",
    "use_source_button": True,
    "use_download_button": True,
    "use_repository_button": True,
    "use_issues_button": True,
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/tlint101/Pharmacophore-Toolkit",
            "icon": "fa-brands fa-square-github",
            "type": "fontawesome",
        },
    ],
}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_book_theme'
html_static_path = ['img']
