#!/usr/bin/env python
# encoding: utf-8
"""
setup.py

Created by FI$H 2000 on 2012-06-19.
Copyright (c) 2012 Objects In Space And Time, LLC. All rights reserved.
"""

__author__ = 'Alexander Bohn'
__version__ = (0, 2, 0)


from distutils.core import setup, Extension
from setuptools import find_packages


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
    name='django-signalqueue',
    version='%s.%s.%s' % __version__,
    description='Truly asynchronous signal dispatch for Django!',

    author=__author__,
    author_email='fish2000@gmail.com',
    maintainer=__author__,
    maintainer_email='fish2000@gmail.com',

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

    namespace_packages=['CIL'],
    packages=find_packages(),
    include_package_data=True,
    package_data={
        'CIL': [
            'ext/*.cpp',
            'ext/*.h',
            'ext/include/*'],
    },
    install_requires=[
        'django-delegate>=0.2.2', 'tornado', 'tornadio2', 'redis',
    ],
    tests_require=[
        'nose', 'rednose', 'django-nose',
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
    ]
)

