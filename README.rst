Gauss2D
###########

.. todo image:: https://img.shields.io/pypi/v/gauss2d.svg
   .. todo   :target: https://pypi.python.org/pypi/gauss2d

.. todo image:: https://img.shields.io/pypi/pyversions/gauss2d.svg
   .. todo   :target: https://pypi.python.org/pypi/gauss2d

*gauss2d* is a package for defining and evaluating 2D Gaussian mixtures and images thereof. It is being
developed primarily for use in astronomy - specifically, by
`Vera C. Rubin Observatory Data Management <https://www.lsst.org/about/dm>`_ for the
`Legacy Survey of Space and Time <https://www.lsst.org/about>`_, and the
`MultiProFit <https://github.com/lsst-dm/multiprofit/>`_ source modelling package - but it can be used for
any kind of image. Future versions will attempt to be entirely generic and avoid astronomy-specific terms.

*gauss2d* requires Python 3, along with `pybind11 <https://github.com/pybind/pybind11>`_ for C++ bindings.
It can be installed using setup.py like so:

python3 setup.py install --user

.. todo *gauss2d* is available in `PyPI <https://pypi.python.org/pypi/multiprofit>`_
   .. and thus can be easily installed via::

.. pip install multiprofit
