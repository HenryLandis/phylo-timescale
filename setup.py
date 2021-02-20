#!/usr/bin/env python

"""
Call `pip install -e .` to install package locally for testing.
"""

from setuptools import setup

# Build setup.
setup(
    name="phylotimescale",
    version="0.0.1",
    author="Henry Landis and Deren Eaton",
    author_email="hnl2109@columbia.edu",
    description="A package to evaluate node dating methods for the inference of accurate phylogenetic trees.",
    entry_points={
        "console_scripts": [
            "pht = phylotimescale.__main__:main", # command = package:file:function
            ]
    },
)