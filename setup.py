# vim: fdm=indent
'''
author:     Taylor Kessinger, Richard Neher
date:       14/11/14
content:    Setup script for betatree.
'''
import os
from setuptools import setup


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


setup(
    name = "betatree",
    version = "0.1",
    author = "Taylor Kessinger and Richard Neher",
    author_email = "richard.neher@tuebingen.mpg.de",
    description = ("A collection of python scripts to generate and analyze "
                   "trees from the beta coalescent ensemble"),
    license = "MIT",
    keywords = "coalescent evolution",
    url = "https://github.com/neherlab/betatree",
    packages=['betatree'],
    long_description=read('README'),
    install_requires = [
        'biopython>=1.66',
        'numpy>=1.10.4',
        'scipy>=0.16.1'
    ],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Science",
        "License :: OSI Approved :: MIT License",
    ],
)
