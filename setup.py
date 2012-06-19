#!/usr/bin/env python
# encoding: utf-8
"""
setup.py

Created by FI$H 2000 on 2012-06-19.
Copyright (c) 2012 Objects In Space And Time, LLC. All rights reserved.
"""

from distutils.core import setup, Extension

try:
    import numpy
except ImportError:
    print("NumPy failed to import!")
    import sys
    sys.exit(1)

ext_modules =   [
                    Extension('CIL.ext.rastersystem', sources=['CIL/ext/rastersystem.cpp'],
                        language='c++',
                        include_dirs=[
                            numpy.get_include(),
                            '/usr/local/include/gsl',
                            '/usr/local/include',
                            '/usr/X11/include',
                            '/usr/include'],
                        library_dirs=[
                            '/usr/local/lib',
                            '/usr/X11/lib',
                            '/usr/lib'],
                        libraries=[
                            'pHash','hdf5','hdf5_cpp','gsl',
                            'boost','adolc','pthread'],
                        extra_compile_args=[
                            '-std=c++0x',
                            "-PIC", "-dynamiclib", "-march=core2",
                            "-msse4.1", "-pipe", "-frounding-math"],
                        ),
                ]

setup(

    name='pHash',
    version= '0.1.0',
    description='C++ imaging power with Pythonic aplomb.',
    ext_modules=ext_modules)
