#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from setuptools import setup, find_packages

# NOTE: disease is intended to be replaced here at some point in the future with a build script that will 
# dynamically create and move code into a subfolder within optima/ for a specific disease area. 
disease = 'tb'
with open("optima/%s/_version.py"%disease, "r") as f:
    version_file = {}
    exec(f.read(), version_file)
    version = version_file["__version__"]

try:
    from pypandoc import convert
except ImportError:
    import io
    def convert(filename, fmt):
        with io.open(filename, encoding='utf-8') as fd:
            return fd.read()

CLASSIFIERS = [
    'Environment :: Console',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: GPLv3',
    'Operating System :: OS Independent',
    'Programming Language :: Python',
    'Topic :: Software Development :: Libraries :: Python Modules',
    'Development Status :: 3 - Alpha',
    'Programming Language :: Python :: 2.7',
]

setup(
    name='optima.tb',
    version=version,
    author='David Kedziora, Sarah Jarvis, Azfar Hussain',
    author_email='info@optimamodel.com',
    description='Software package for modeling H2H infectious disease epidemics',
    #long_description=convert('README.md', 'md'),
    url='http://github.com/optimamodel/tb-ucl',
    keywords=['optima','disease'],
    platforms=['OS Independent'],
    classifiers=CLASSIFIERS,
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'matplotlib>=1.4.2',
        'numpy>=1.10.1',
    ],
)
