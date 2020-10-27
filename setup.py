# -*- coding: utf-8 -*-
#!/usr/bin/env python

# Learn more: https://github.com/kennethreitz/setup.py

from setuptools import setup, find_packages


with open('README.rst') as f:
    readme = f.read()

# with open('LICENSE') as f:
#     license = f.read()

setup(
    name='oci_l1a',
    version='0.1.0',
    description="Package to read Level 1A data from PACE OCI",
    long_description=readme,
    author="Christopher Field",
    author_email="Christopher.T.Field@nasa.gov",
    url='Not yet defined',
    # license=license,
    packages=find_packages(exclude=('tests', 'docs'))
)
