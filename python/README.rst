Gauss2D
#######

.. todo image:: https://img.shields.io/pypi/v/gauss2d.svg
.. todo   :target: https://pypi.python.org/pypi/gauss2d

.. todo image:: https://img.shields.io/pypi/pyversions/gauss2d.svg
.. todo   :target: https://pypi.python.org/pypi/gauss2d

*gauss2d* provides python bindings for the C++ libgauss2d library.
Python 3.10+, `pybind11 <https://github.com/pybind/pybind11>`_, and numpy are
requirements, as numpy arrays are used to store images. The library may build
on Python version as old as 3.8 but these are supported on a best-effort basis.

With libgauss2d installed, these bindings can be built using meson:

Build scripts are provided in the base directory of gauss2d.

.. todo *gauss2d* is available in `PyPI <https://pypi.python.org/pypi/multiprofit>`_
.. and thus can be easily installed via::

.. pip install multiprofit
