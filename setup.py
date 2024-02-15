#!/usr/bin/env python3
# coding: utf-8

import io
import os
import glob

from setuptools import find_packages
from setuptools import setup

# name of package ElliptiCBn

#Package meta-data
DESCRIPTION= \
""" Python package for analyzing ellipticity in CBn crystal structures automatically. """
URL = "https://github.com/harmslab/ElliptiCBn"  # temp URL
EMAIL = "mshavlik@uoregon.edu"
AUTHOR = "Michael Shavlik"
REQUIRES_PYTHON = ">=3.7.0"
VERSION = "1.2.0"

here = os.path.abspath(os.path.dirname(__file__))

# Import README for description
with io.open(os.path.join(here,'README.md'),encoding='utf-8') as f:
    full_description = '\n' + f.read()

packages = find_packages()
package_data = {"":["ElliptiCBn/data/*.csv",
                    "ElliptiCBn/data/*.txt",
                    "notebooks/*.ipynb"]}
    
# Now the part where we do setup
setup(
    name='ElliptiCBn',
    version='1.2.0',
    author=AUTHOR,
    author_email=EMAIL,
    description=DESCRIPTION,
    long_description=full_description,
    url=URL,
    license='MIT',
    packages=packages,
    package_data=package_data,
    include_package_data=True,
    scripts=glob.glob("bin/*"),
    keywords='CBn; cucurbituril; host; guest; chemistry; molecule; science; analysis; crystal structure ',
    classifiers = ["Development Status :: 3 - Alpha",
                  'Intended Audience :: Science/Research',
                  'Programming Language :: Python :: 3.7'
                  'Programming Language :: Python :: 3.8',
                  'Programming Language :: Python :: 3.9',
                  'Programming Language :: Python :: 3.10',
                  'Programming Language :: Python :: 3.11'],
    zip_safe=False)