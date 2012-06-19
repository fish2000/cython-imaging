#!/usr/bin/env python
# encoding: utf-8
"""
setup.py

Created by FI$H 2000 on 2012-06-19.
Copyright (c) 2012 Objects In Space And Time, LLC. All rights reserved.

"""

__author__ = 'Alexander Bohn'
__contact__ = 'fish2000@gmail.com'
__version__ = (0, 2, 0)


from setuptools import find_packages
from os.path import join, dirname, abspath
from distutils.core import setup, Extension
from distutils.sysconfig import get_python_inc


try:
    import numpy
except ImportError:
    print("NumPy failed to import!")
    import sys
    sys.exit(1)

try:
    pth = abspath(__file__)
except (ValueError, AttributeError):
    import os
    pth = os.getcwd()

ext_pathex = lambda *pth: join('CIL', 'ext', 'src', *pth)
ext_one = Extension('CIL.ext.rastersystem',
            sources=[ext_pathex('rastersystem.cpp')],
            language='c++',

            include_dirs=[
                join(dirname(abspath(__file__)), 'CIL', 'ext', 'include'),
                numpy.get_include(),
                get_python_inc(plat_specific=1),
                '/usr/local/include/gsl',
                '/usr/local/include',
                '/usr/X11/include',
                '/usr/include'],

            library_dirs=[
                '/usr/local/lib',
                '/usr/X11/lib',
                '/usr/lib'],

            libraries=['stdc++',
                'pHash', 'hdf5', 'hdf5_cpp',
                'gsl', 'adolc', 'pthread',
                'boost_system-mt', 'boost_regex-mt', 'boost_thread-mt'],

            extra_compile_args=[
                '-std=c++0x',
                "-PIC", "-dynamiclib", "-march=core2",
                "-msse4.1", "-pipe", "-frounding-math"],

            depends=[
                ext_pathex('rastersystem.cpp')])



setup(
    name='cython-imaging',
    version='%s.%s.%s' % __version__,
    description='C++ imaging power with Pythonic aplomb',
    long_description=""" C++ imaging power, with Pythonic aplomb. """,

    author=__author__,
    author_email=__contact__,
    maintainer=__author__,
    maintainer_email=__contact__,

    license='BSD',
    url='http://github.com/fish2000/django-signalqueue/',
    keywords=[
        'imaging',
        'CImg', 'C++',
        'CBIR', 'HDF5',
        'image analysis',
        'image comparison',
        'image processing',
        'perceptual hash',
        'Cython',
        'Ctypes',
    ],
    
    packages=find_packages(),
    
    ext_modules=[ext_one],

    install_requires=[
        'numpy', 'scipy', 'h5py', 'ipython',
    ],
    
    tests_require=[
        'nose', 'rednose',
    ],

    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Web Environment',
        'Framework :: Django',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 2.7',
        'Topic :: Utilities'
    ],
)
