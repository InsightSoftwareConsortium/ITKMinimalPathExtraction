ITKMinimalPathExtraction
========================

.. image:: https://circleci.com/gh/InsightSoftwareConsortium/ITKMinimalPathExtraction.svg?style=svg
    :target: https://circleci.com/gh/InsightSoftwareConsortium/ITKMinimalPathExtraction

This is an `ITK <http://itk.org>`_ module that implements a minimal path
extraction framework based on Fast Marching arrival functions.

A more detailed description can be found in
`the Insight Journal article <http://hdl.handle.net/1926/1332>`_::

  Mueller, D. "Fast Marching Minimal Path Extraction in ITK"
  http://hdl.handle.net/1926/1332
  http://www.insight-journal.org/browse/publication/213
  March, 2008.

Since ITK 4.8.0, this module is available in the ITK source tree as a Remote
module.  To enable it, set::

  Module_MinimalPathExtraction:BOOL=ON

in ITK's CMake build configuration.
