"""
Repertoire Zoo

A Collection of animal-themed repertoire analysis tools.

See README.md for details on usage.
"""

from setuptools import setup, find_packages

setup(
    name='repertoire_zoo',
    version='0.0.1',
    description='A set of python modules for repertoire analysis.',
    author='Samuel Wollenburg',
    author_email='Samuel.Wollenburg@UTSouthwestern.edu',
    packages=find_packages(),
    python_requires='>=3.6',
    license='GPLv3',
    install_requires=[
        'airr',
        'matplotlib',
        'numpy',
        'pandas',
        'pycirclize',
        'seaborn',
    ],
    classifiers=[
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent',
    ],
)
