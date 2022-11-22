
Gauss2D
#######

*gauss2d* is a package for defining and evaluating 2D Gaussian mixtures and 
images thereof. It is being developed primarily for use in astronomy - 
specifically, by
`Vera C. Rubin Observatory Data Management <https://www.lsst.org/about/dm>`_ 
for the `Legacy Survey of Space and Time <https://www.lsst.org/about>`_,
and the
`MultiProFit <https://github.com/lsst-dm/multiprofit/>`_ source modelling 
package - but it can be used for any kind of image or domain.

*gauss2d* is provided in two parts; a C++ shared library (libgauss2d) and 
optional Python bindings. The Python bindings include a shared library for
the convenience of developers who wish to generate bindings for their own
templated Image class (e.g. based on something other than numpy arrays).

*gauss2d* is built with `Meson <https://github.com/mesonbuild>`_ and uses
`pkg-config <https://github.com/pkgconf/pkgconf>`_ to store package
information. You will likely want to configure it to install in a local
directory. For example, if you are using a conda environment:

``PKG_CONFIG_PATH=$CONDA_PREFIX/.local/lib64/pkgconfig meson 
--prefix=$CONDA_PREFIX/.local build``

You will also need to run ``meson build`` once, followed by
``meson compile -C build && meson install -C build`` to (re)compile.
This must be done in the base directory and then the ``python/lib``
and ``python`` subfolders (in that order) if bindings are desired.

A full example setup script is provided in ``setup_conda-release.sh``.
This defaults to  using ``$CONDA_PREFIX``, but can be configured to output
elsewhere (e.g. ``~/.local``):

``CONDA_PREFIX=~ sh setup-conda-release.sh``

Once the build command is run once to create the build directories, subsequent
rebuilds can use the provided ``build*.sh`` scripts.

One can also build a configuration for debugging purposes: see
``setup-debug.sh`` for an example. Note that meson defaults to a build
directory named ``builddir``, to disambiguate the  common ``meson build``
command. The example scripts configure separate ``build-debug`` and 
``build-release`` directories. If you are using an IDE that expects to find a 
 ``build``, simply create a symbolic link (e.g. ``ln -s build-debug build``).

Note: You may need to change include and/or LD paths to debug the Python
bindings with ``build-debug``; this has not been tested yet.
