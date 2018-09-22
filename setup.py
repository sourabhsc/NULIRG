from distutils.core import setup
from io import open
from setuptools import setup, find_packages
from os import path
here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='NULIRG',
    version='0.1dev',
    description='NULIRG documentation',
    author='Sourabh Chauhan',
    packages=['NULIRG'],
    include_package_data=True,
    author_email=['sourabh@astro.umn.edu'],
    package_data={'NULIRG': ['data/*', 'data_dark/*']},
    license='Creative Commons Attribution-Noncommercial-Share Alike license',
    long_description=open('README.rst').read(),
)
