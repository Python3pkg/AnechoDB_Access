# -*- coding: utf-8 -*-
"""
@author: Jacopo Martelli
"""

import setuptools

setuptools.setup(
    name='AnechoDB_Access',
    version='1.0',
    description='A library to comunicate with the database Belen to search, download and compute beam stored in it',
    author='Jacopo Martelli',
    author_email='martelli.jacopo@gmail.com',
    url='https://github.com/JacopoMartelli/AnechoDB_Access',
    requires=['Python (>3.3)'],
    license='MIT',
    long_description=open('README.rst').read(),
    packages=['share_belen',],
    #test_suite='test_sb',
    install_requires=['h5py'],
)
