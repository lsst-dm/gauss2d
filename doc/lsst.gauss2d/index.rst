.. py:currentmodule:: lsst.gauss2d

.. _lsst.gauss2d:

############
lsst.gauss2d
############

gauss2d provides classes and methods for defining 2D Gaussian mixtures and
evaluating their (approximate) integrals over square pixels. Several ellipse
parameterizations are implemented but evaluations are done strictly with the
(sigma_x, sigma_y, rho) variant, since this keeps rho bounded between -1 and 1
and not periodic (like a position angle would be).

gauss2d can also evaluate the first derivatives of a model (i.e. the Jacobian)
or its likelihood analytically.

.. _lsst.gauss2d-using:

Using lsst.gauss2d
==============================

Example usage can be found in the unit tests and also in dependent packages,
particularly `gauss2dfit <https://github.com/lsst-dm/gauss2dfit>`_.

.. toctree::
   :maxdepth: 2

.. _lsst.gauss2d-contributing:

Contributing
============

``lsst.gauss2d`` is developed at https://github.com/lsst-dm/gauss2d.
You can find Jira issues for this module under the
`gauss2d <https://rubinobs.atlassian.net/issues/?jql=project%20%3D%20DM%20AND%20component%20%3D%20gauss2d>`_
component.

.. If there are topics related to developing this module (rather than using it), link to this from a toctree placed here.

.. .. toctree::
..    :maxdepth: 2

.. _lsst.gauss2d-pyapi:

Python API reference
====================

``lsst.gauss2d`` has Python bindings for classes using numpy-based single
and double precision arrays. Support for GSL arrays is forthcoming with
`DM-38617 <https://jira.lsstcorp.org/browse/DM-38617>`_.

.. automodapi:: lsst.gauss2d
   :no-main-docstr:
   :no-inheritance-diagram:
