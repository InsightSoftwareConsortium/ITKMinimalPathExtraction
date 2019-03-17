ITKMinimalPathExtraction
========================

.. image::  https://dev.azure.com/InsightSoftwareConsortium/ITKModules/_apis/build/status/InsightSoftwareConsortium.ITKMinimalPathExtraction?branchName=master
    :target: https://dev.azure.com/InsightSoftwareConsortium/ITKModules/_build/latest?definitionId=16&branchName=master
    :alt: Build Status


Overview
--------

This is an `ITK <http://itk.org>`_ module that implements a minimal path
extraction framework based on Fast Marching arrival functions.

A more detailed description can be found in
`the Insight Journal article <http://hdl.handle.net/1926/1332>`_::

  Mueller, D.
  "Fast Marching Minimal Path Extraction in ITK"
  The Insight Journal. January-June, 2008.
  http://hdl.handle.net/1926/1332
  http://www.insight-journal.org/browse/publication/213


Installation
------------

Since ITK 4.8.0, this module is available in the ITK source tree as a remote
module. To enable it, set::

  Module_MinimalPathExtraction:BOOL=ON

in ITK's CMake build configuration.


License
-------

This software is distributed under the Apache 2.0 license. Please see
the *LICENSE* file for details.
