# -*- coding: utf-8 -*-
##########################################################################
# pySAP - Copyright (C) CEA, 2017 - 2018
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# Module current version
version_major = 0
version_minor = 0
version_micro = 1

# Expected by setup.py: string of form "X.Y.Z"
__version__ = "{0}.{1}.{2}".format(version_major, version_minor, version_micro)

# Expected by setup.py: the status of the project
CLASSIFIERS = ["Development Status :: 1 - Planning",
               "Environment :: Console",
               "Operating System :: OS Independent",
               "Programming Language :: Python",
               "Topic :: Scientific/Engineering"]

# Project descriptions
description = """
PYthon CosmoStat Package
"""
SUMMARY = """
.. container:: summary-carousel

    pycs is a Python module by **CosmoStat lab **

"""
long_description = (
    "CosmoStat development package"
)
# Main setup parameters
NAME = "python-pycs"
ORGANISATION = "CEA"
MAINTAINER = "Sam Farrens"
MAINTAINER_EMAIL = "sam.farrens@cea.fr"
DESCRIPTION = description
LONG_DESCRIPTION = long_description
EXTRANAME = "CosmoStat webPage"
EXTRAURL = "http://cosmic.cosmostat.org/"
URL = "https://github.com/CosmoStat/cosmostat"
DOWNLOAD_URL = "https://github.com/CosmoStat/cosmostat"
LICENSE = "CeCILL-B"
CLASSIFIERS = CLASSIFIERS
AUTHOR = """
Francois Lanusse <francois.lanusse@cea.fr>
Samuel Farrens <samuel.farrens@cea.fr>
Jean-Luc Starck <jl.starck@cea.fr>
"""
AUTHOR_EMAIL = "sam.farrens@cea.fr"
PLATFORMS = "Linux,OSX"
ISRELEASE = True
VERSION = __version__
PROVIDES = ["pycs"]
REQUIRES = [
    "numpy>=1.16.4",
]

PREINSTALL_REQUIRES = [
    "pybind11>=2.3.0"
]
