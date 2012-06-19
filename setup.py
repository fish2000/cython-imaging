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

# You have to import stuff from setuptools before
# you import from distutils -- and by "before" I mean
# "from a line in the file nearer to the top" --
# or the act of doing the import will somehow cause
# setuptools to sabotage the distutilitarian imports
# ... notably the Extension class won't work. WTF.
# that is srsly some unpythonic shit, you guys.
# get it together.

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
    septh = dirname(abspath(__file__))
except (ValueError, AttributeError):
    # not reliably this directory
    import os
    septh = os.getcwd() 


pathex =        lambda *pth: join(septh, *pth)
ext_pathex =    lambda *pth: join('CIL', 'ext', 'src', *pth)

ext_modules = [
        
        Extension('rastersystem',
            sources=[ext_pathex('rastersystem.cpp')],
            language='c++',

            include_dirs=[
                pathex('CIL', 'ext', 'include'),
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
                'boost_system-mt',
                'boost_thread-mt',
                'boost_regex-mt'],

            extra_compile_args=[
                '-std=c++0x', "-O4", "-DNDEBUG"
                "-PIC", "-dynamiclib", "-march=core2",
                "-msse4.1", "-pipe", "-frounding-math"],

            depends=[
                pathex('setup.py'),
                ext_pathex('rastersystem.cpp')]
            
            ),
        
        ]


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
    namespace_packages=['CIL'],
    include_package_data=True,
    package_data={
        'CIL': [
            'ext/*.cpp',
            'ext/*.h',
            'ext/src/*',
            'ext/include/*'],
    },

    ext_package='CIL.ext',
    ext_modules=ext_modules,

    setup_requires=[
        'cython', 'ctypes'
    ],

    install_requires=[
        'numpy', 'scipy', 'h5py', 'ipython',
        'PIL', 'ctypes'
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
