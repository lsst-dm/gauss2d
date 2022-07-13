
Gauss2D
#######

*gauss2d* is a package for defining and evaluating 2D Gaussian mixtures and 
images thereof. It is being developed primarily for use in astronomy - 
specifically, by
`Vera C. Rubin Observatory Data Management <https://www.lsst.org/about/dm>`_ 
for the `Legacy Survey of Space and Time <https://www.lsst.org/about>`_, and the
`MultiProFit <https://github.com/lsst-dm/multiprofit/>`_ source modelling 
package - but it can be used for any kind of image or domain.

*gauss2d* is provided in two parts; a C++ shared library (libgauss2d) and 
optional Python bindings.

*gauss2d* is built with `Meson <https://github.com/mesonbuild>`_ and uses
`pkg-config <https://github.com/pkgconf/pkgconf>`_ to store package information.
You will likely want to configure it to install in a local directory. For example,
 if you are using a conda environment:

``PKG_CONFIG_PATH=$CONDA_PREFIX/.local/lib64/pkgconfig meson --prefix=$CONDA_PREFIX/.local build``

You will also need to run ``meson build`` once in the ``python/lib`` folder,
run ``meson compile -C build && meson install -C build`` once, and then
``meson build`` in the  ``python`` folder to build the Python bindings.

A full example setup script is provided in ``setup_conda.sh``.

Once the build command is run once to create the build directories, subsequent
rebuilds can use the provided ``build*.sh`` scripts.
