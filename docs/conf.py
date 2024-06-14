# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.append ( os.path.abspath ('../HandyTools') )
sys.path.append ( os.path.abspath ('../DataTools') )
sys.path.append ( os.path.abspath ('../PYtoCPP/DataTools') )

from HandyTools import HandyTools

autodoc_mock_imports = ['DataWranglingToolsPYtoCPP', 'FilterToolsPYtoCPP']
from DataTools import DataTools


# -- Project information -----------------------------------------------------

project = 'GeneralTools for Scientists - Documentation'
author = 'Maarten Roos-Serote'

# The full version, including alpha/beta/rc tags
release = '1.0'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx_rtd_theme',
    'sphinx.ext.autodoc',
]

autodoc_default_options = {
    'member-order': 'bysource',
}

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# Automatic numbering of figures and tables.
numfig = False
numtab = False


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    "stickysidebar": "true",
    "description": "General Tools for Scientitsts"
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# I added this based on the entry of 4 August 2022 in https://github.com/readthedocs/sphinx_rtd_theme/issues/301
# The issue was that the sphinx_rtd_theme has a bug of placing the equation numbers above te equations, and not to the right of them.
# This seems to solve it, or rather, it is a work-around.
# html_css_files = ['css/custom.css']

html_last_updated_fmt = '%b %d, %Y, %X'
