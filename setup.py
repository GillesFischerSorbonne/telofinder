# -*- coding: utf-8 -*-
# License: 3-clause BSD
__revision__ = "$Id: $"  # for the SVN Id
import sys
import os
from setuptools import setup, find_packages
from setuptools import setup
import glob

_MAJOR = 0
_MINOR = 1
_MICRO = 0
version = "%d.%d.%d" % (_MAJOR, _MINOR, _MICRO)
release = "%d.%d" % (_MAJOR, _MINOR)


metainfo = {
    "authors": {"main": ("XX", "XX")},
    "maintainer": {"main": ("XX", "XX")},
    "version": version,
    "license": "new BSD",
    "download_url": "https://github.com/GillesFischerSorbonne/dubii_project".format(
        version
    ),
    "url": "http://github.com/GillesFisherSoronne/dubii_project",
    "description": "A library for telomere prediction.",
    "platforms": ["Linux", "Unix"],
    "keywords": ["NGS", "telomers"],
    "classifiers": [
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Information Analysis",
        "Topic :: Scientific/Engineering :: Mathematics",
        "Topic :: Scientific/Engineering :: Physics",
    ],
}


setup(
    name="telofinder",
    version=version,
    maintainer=metainfo["authors"]["main"][0],
    maintainer_email=metainfo["authors"]["main"][1],
    author=metainfo["authors"]["main"][0],
    author_email=metainfo["authors"]["main"][1],
    long_description=open("README.md").read(),
    keywords=metainfo["keywords"],
    description=metainfo["description"],
    license=metainfo["license"],
    platforms=metainfo["platforms"],
    url=metainfo["url"],
    download_url=metainfo["download_url"],
    classifiers=metainfo["classifiers"],
    # package installation
    # package_dir = {'telofinder': 'src/python_script'},
    packages=find_packages(),
    install_requires=["matplotlib", "Bio", "easydev"],
    package_data={"telofinder.data": ["*.*"],},
    zip_safe=False,
    # entry_points = {
    #    'console_scripts':[
    #       'telofinder=src.telofinder.analyze_telom_length:main',
    #    ]
    # },
)
