from setuptools import setup, find_packages
import os

setup(
    name="pycs",
    author="CosmoStat Laboratory",
    author_email="",
    version="0.0.1rc1",
    packages=find_packages(),
    install_requires = [
      'lenspack @ git+https://github.com/CosmoStat/lenspack.git@master#egg=lenspack',
    ]
)
