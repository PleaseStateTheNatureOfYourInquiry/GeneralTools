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
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))

import sphinx_rtd_theme



# -- Project information -----------------------------------------------------

project = 'GeneralTools Documentation'
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

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


math_eqref_format = 'Eq.{number}'


# COMMENTAAR MAARTEN: met onderstaande 'boolean' variabelen kan je de automatische nummering van
#  figuren en tabellen aan (True) of uit (False) zetten.
#  Een 'boolean' variabele werkt precies als een schakelaartje, denk aan een lichtschakelaar bijvoorbeeld: het staan aan (True of 1) of uit (False of 0).
#  Ik heb het nu even uit gezet. Is dit handig/mooier?? 

# In een Python script, wat dit bestand 'conf.py' is (daarom heeft het ook de uitgang .py), moet commentaar worden voorafgegaan door een #-teken. 
#  Dat is dus anders dan in een .rst bestand!

numfig = False
numtab = False


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'
# html_theme_options = {
#     "stickysidebar": "true",
#     "description": "The reference documentation for all EMC-EFO software."
# }

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# I added this based on the entry of 4 August 2022 in https://github.com/readthedocs/sphinx_rtd_theme/issues/301
# The issue was that the sphinx_rtd_theme has a bug of placing the equation numbers above te equations, and not to the right of them.
# This seems to solve it, or rather, it is a work-around.
html_css_files = ['css/custom.css']

html_last_updated_fmt = '%b %d, %Y, %X'
