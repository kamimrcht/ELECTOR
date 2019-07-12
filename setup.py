#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

from itertools import tee

setup(
    name = "elector",
    version = "0.9",
    packages = find_packages(),

    author = "Camile Marchet",
    author_email = "marchetcamille@gmail.com",
    description = "A tools to evalute long-reads correction tools",
    long_description = open('README.md').read(),
    url = "https://github.com/kamimrcht/ELECTOR",
    
    install_requires = ['biopython'],
    include_package_data = True,
    
    classifiers = [
        "Programming Language :: Python :: 3",
        "Development Status :: 2 - Pre-Alpha"
    ],

    entry_points = {
        'console_scripts': [
            'elector = elector.__main__:main',
        ]
    }
)
