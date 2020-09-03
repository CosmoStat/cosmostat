# -*- coding: utf-8 -*-

# Module current version
version_major = 0
version_minor = 0
version_micro = 1

# Expected by setup.py: string of form "X.Y.Z"
__version__ = "{0}.{0}.{1}".format(version_major, version_minor, version_micro)

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

    pycs is a Python package developed by **CosmoStat**

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
Samuel Farrens <samuel.farrens@cea.fr>
Francois Lanusse <francois.lanusse@cea.fr>
Jean-Luc Starck <jl.starck@cea.fr>
"""
AUTHOR_EMAIL = "sam.farrens@cea.fr"
PLATFORMS = "Linux,OSX"
ISRELEASE = True
VERSION = __version__
PROVIDES = ["pycs"]
CMAKE_VERSION = "3.12.0"
REQUIRES = [
    "astropy>=4.0",
    "healpy>=1.14",
    # "lenspack" // wait for lenspack release on PyPi
    "matplotlib>=3.3"
    "modopt>=1.4.4",
    "numpy>=1.19",
    "pyqt5>=5.12.2",
    "pyqtgraph>=0.11",
    "scipy>=1.5",
    "seaborn>=0.10"
]

PREINSTALL_REQUIRES = [
    "pybind11>=2.5"
]
