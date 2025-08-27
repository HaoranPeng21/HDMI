#!/usr/bin/env python3
"""
Setup script for HDMI (HGT Detection from MAGs in Individual)
"""

from setuptools import setup, find_packages
import os

# Read the README file
def read_readme():
    with open("README.md", "r", encoding="utf-8") as fh:
        return fh.read()

# Get version
def get_version():
    return "1.0.0"

setup(
    name="HDMI",
    version=get_version(),
    author="Haoran Peng",
    author_email="penghr21@gmail.com",
    description="HGT Detection from MAGs in Individual - A comprehensive pipeline for detecting horizontal gene transfer events",
    long_description=read_readme(),
    long_description_content_type="text/markdown",
    url="https://github.com/HaoranPeng21/HDMI",
    packages=find_packages(),
    py_modules=[],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
    python_requires=">=3.7",
    install_requires=[
        "pandas>=1.3.0",
        "numpy>=1.21.0",
        "biopython>=1.79",
        "pysam>=0.19.0",
    ],
    entry_points={
        "console_scripts": [
            "HDMI=HDMI:main",
        ],
    },
    include_package_data=True,
    package_data={
        "": ["*.py"],
    },
    scripts=[
        "bin/HDMI",
    ],
    keywords="bioinformatics, metagenomics, horizontal gene transfer, HGT, MAGs",
    project_urls={
        "Bug Reports": "https://github.com/HaoranPeng21/HDMI/issues",
        "Source": "https://github.com/HaoranPeng21/HDMI",
        "Documentation": "https://github.com/HaoranPeng21/HDMI#readme",
    },
)
