# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation: # https://www.sphinx-doc.org/en/master/usage/configuration.html
# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('./..'))


# -- Project information -----------------------------------------------------

project = 'Pipe Stress Infinity'
copyright = '2021, Denis Gomes'
author = 'Denis Gomes'

# The full version, including alpha/beta/rc tags
release = '0.0.1'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.napoleon',
              'sphinx.ext.mathjax',
              'jupyter_sphinx.execute',
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '**.ipynb_checkpoints']

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'alabaster'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']


html_sidebars = {
    '**': [
        'about.html',
        'navigation.html',
        # 'relations.html',
        'searchbox.html',
        # 'donate.html',
    ]
}


html_theme_options = {
    'logo': 'logo.png',
    # 'logo_name': True,
    # 'logo_text_align': 'left',
    'github_user': 'denisgomes',
    'github_repo': 'psi',
    'description': 'The pipe stress analysis and design software.',
    # 'fixed_sidebar': True,
    # 'github_banner': True,
    # 'github_button': True,
    'show_relbars': True,
    # 'font_size': 8,
    'extra_nav_links':
    {'PSI @ PyPI': 'https://pypi.org',
     'PSI @ GitHub': 'https://github.com/denisgomes/psi',
     'PSI @ Discord': 'https://discord.gg/xnHnwbD',
     'Issue Tracker': 'https://github.com/denisgomes/psi/issues',
     'Mailing List': 'https://groups.google.com/group/pipestressinfinity-users',
     'PDF Documentation': 'https://readthedocs.org/projects/pipe-stress-infinity/downloads/pdf/latest',
    }
}


# for readthedocs build
master_doc = 'index'

# for figure numbering
# numfig = True

# mathjax_config = {
#     "jax": ["input/TeX","output/HTML-CSS"],
#     "displayAlign": "left"
# }
