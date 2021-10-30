Gauss2D
#######

.. todo image:: https://img.shields.io/pypi/v/gauss2d.svg
   .. todo   :target: https://pypi.python.org/pypi/gauss2d

.. todo image:: https://img.shields.io/pypi/pyversions/gauss2d.svg
   .. todo   :target: https://pypi.python.org/pypi/gauss2d

*gauss2d* provides python bindings for the C++ libgauss2d library.
Python 3.7+, `pybind11 <https://github.com/pybind/pybind11>`_, and numpy are 
requirements, as numpy arrays are used to store images.

With libgauss2d installed, these bindings can be built using meson:

meson --prefix=~/.local build && meson compile -C build && meson install -C build

... where the prefix argument is your desired installation path.

You may also need to configure pkg-config, e.g. prepend 
PKG_CONFIG_PATH=~/.local/lib64/pkgconfig to the build command for a local
and/or non-root installation.

Alternatively, it can be built for local development using poetry like so:

poetry build && python3 -m pip install .

... using the include pyproject.toml.

.. todo *gauss2d* is available in `PyPI <https://pypi.python.org/pypi/multiprofit>`_
   .. and thus can be easily installed via::

.. pip install multiprofit
