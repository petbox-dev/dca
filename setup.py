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
from re import sub

from petbox.dca import __version__

try:
    from setuptools import setup  # type: ignore
except ImportError:
    from distutils.core import setup


def get_long_description() -> str:
    # Fix display issues on PyPI caused by RST markup
    readme = open('README.rst').read()

    version_lines = []
    with open('docs/versions.rst') as infile:
        next(infile)
        for line in infile:
            line = line.rstrip().replace('.. automodule:: more_itertools', '')
            version_lines.append(line)
    version_history = '\n'.join(version_lines)
    version_history = sub(r':func:`([a-zA-Z0-9._]+)`', r'\1', version_history)

    ret = readme + '\n\n' + version_history
    return ret


if sys.argv[-1] == 'build':
    print('\nBuilding...')
    os.system('rm -r dist\\')
    os.system('python setup.py sdist bdist_wheel')


setup(
    name='petbox-dca',
    version=__version__,
    description='Decline Curve Library',
    long_description=get_long_description(),
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
