# -*- coding: utf-8 -*-
from __future__ import print_function
from os import sys

try:
    from skbuild import setup
except ImportError:
    print('scikit-build is required to build from source.', file=sys.stderr)
    print('Please run:', file=sys.stderr)
    print('', file=sys.stderr)
    print('  python -m pip install scikit-build')
    sys.exit(1)

setup(
    name='itk-minimalpathextraction',
    version='1.0.2',
    author='Insight Software Consortium',
    author_email='itk+community@discourse.itk.org',
    packages=['itk'],
    package_dir={'itk': 'itk'},
    download_url=r'https://github.com/InsightSoftwareConsortium/ITKMinimalPathExtraction',
    description=r'A minimal path extraction framework based on Fast Marching arrival functions.',
    long_description='itk-minimalpathextraction provides a minimal path '
                     'extraction framework based on Fast Marching arrival '
                     'functions.\n'
                     'Please refer to:\n'
                     'Mueller, D. "Fast Marching Minimal Path Extraction in ITK", '
                     'Insight Journal, January-June 2008, http://hdl.handle.net/1926/1332.',
    classifiers=[
        "License :: OSI Approved :: Apache Software License",
        "Programming Language :: Python",
        "Programming Language :: C++",
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "Intended Audience :: Education",
        "Intended Audience :: Healthcare Industry",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Medical Science Apps.",
        "Topic :: Scientific/Engineering :: Information Analysis",
        "Topic :: Software Development :: Libraries",
        "Operating System :: Android",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: POSIX",
        "Operating System :: Unix",
        "Operating System :: MacOS"
        ],
    license='Apache',
    keywords='ITK InsightToolkit',
    url=r'https://github.com/InsightSoftwareConsortium/ITKMinimalPathExtraction',
    install_requires=[
        r'itk>=5.0.0.post1'
    ]
    )
