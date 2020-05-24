"""
Decline Curve Models
Originally developed for David S. Fulford's thesis research

Author
------
David S. Fulford
Derrick W. Turk

Notes
-----
Created on August 5, 2019
"""

import os
import sys

from petbox.dca import __version__

try:
    from setuptools import setup  # type: ignore
except ImportError:
    from distutils.core import setup

if sys.argv[-1] == 'build':
    print('\nBuilding...')
    os.system('rm -r dist\\')
    os.system('python setup.py sdist bdist_wheel')


setup(
    name='petbox-dca',
    version=__version__,
    description='Decline Curve Library',
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    url='https://github.com/petbox-dev/dca',
    author='David S. Fulford',
    author_email='dsfulford@gmail.com',
    install_requires=['numpy', 'scipy'],
    zip_safe=False,
    packages=['petbox.dca'],
    package_data={
        'petbox.dca': ['py.typed']
    },
)
