# -*- coding: utf-8 -*-
# License: 3-clause BSD
__revision__ = "$Id: $"  # for the SVN Id

from setuptools import setup, find_packages

_MAJOR = 1
_MINOR = 0
_MICRO = 0
version = "%d.%d.%d" % (_MAJOR, _MINOR, _MICRO)
release = "%d.%d" % (_MAJOR, _MINOR)

requirements = open("requirements.txt").read().split()

metainfo = {
    "authors": {
        "main": (
            "Gilles Fischer",
            "gilles.fischer@sorbonne-universite.fr",
            "Etienne Kornobis",
            "Thomas Cokelaer",
        )
    },
    "maintainer": {"main": ("Gilles Fischer")},
    "version": version,
    "license": "new BSD",
    "download_url": "https://github.com/GillesFischerSorbonne/telofinder".format(
        version
    ),
    "url": "http://github.com/GillesFisherSoronne/telofinder",
    "description": "A library for telomere prediction.",
    "platforms": ["Linux", "Unix"],
    "keywords": ["genome assembly", "telomeres", "telomeric repeats"],
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
    install_requires=requirements,
    package_data={
        "telofinder.data": ["*.*"],
    },
    zip_safe=False,
    entry_points={
        "console_scripts": [
            "telofinder=telofinder.main:main",
        ]
    },
)
