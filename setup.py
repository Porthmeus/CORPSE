# Porthmeus
# 26.08.21

from setuptools import find_packages, setup
# To use a consistent encoding
from codecs import open
from os import path

HERE = path.abspath(path.dirname(__file__))
# Get the long description from the README file
with open(path.join(HERE, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='corpse',
    version='0.1.0',
    description='CObRaPy tool SEt - auxiliary functions for cobrapy',
    url = 'https://github.com/Porthmeus/CORPSE',
    author='Jan Taubenheim',
    author_email = "j.taubenheim@iem.uni-kiel.de",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license='GPL3',
    packages=["corpse"],
    install_requires = [
        "cobra>=0.20.0",
        "pandas>=1.1.4",
        "numpy>=1.19.4",
        "joblib>=0.17.0",
        "scipy>=1.5.4",
        "troppo",
        "cobamp"
        ]
)
