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
optional Python bindings. The Python bindings include a sub-package for
the convenience of developers who wish to generate bindings for their own
templated Image class (e.g. based on something other than numpy arrays).

Building the Package
####################

*gauss2d* is built with `Meson <https://github.com/mesonbuild>`_. It can
be built as a standalone package or as an eups-managed Science Pipelines
package.

EUPS build
##########

To build with `eups <https://github.com/RobertLuptonTheGood/eups>`_ for use
with the Rubin/LSST `Science Pipelines <https://pipelines.lsst.io/>`_,
call:

.. code-block:: sh
   :name: build

   setup -jr .
   eupspkg -e -v 1 config
   eupspkg -e -v 1 build

If testing a ticket with Jenkins, see full instructions with more steps in the
`developer guide <https://developer.lsst.io/stack/packaging-third-party-eups-dependencies.html#testing-the-package>`_.

Doc build
#########

Sphinx package documentation is built by meson using
`documenteer <https://github.com/lsst-sqre/documenteer/>`_, if available.
Unfortunately, C++ doctrings are not passed through to pybind11 objects, so
the Doxygen docs are generally more useful.

C++ Doxygen docs are built with scons using
`sconsUtils <https://github.com/lsst/sconsUtils>`_, if available, either by
calling eupspkg install or manually running scons.
sconsUtils' ``tickets/DM-44144`` branch can be pip-installed outside of the
Science Pipelines (the next section).

Standalone builds
#################

A full example setup script is provided in ``setup_conda.sh``.
This defaults to  using ``$CONDA_PREFIX``, but can be configured to output
elsewhere (e.g. ``~/.local``) like so:

``CONDA_PREFIX=~ sh setup-conda.sh``

Once the build command is run once to create the build directories, subsequent
rebuilds can use the provided ``build.sh`` script.

Otherwise, to manually create a build directory, call:

``meson builddir/``

Standalone builds require `pkg-config <https://github.com/pkgconf/pkgconf>`_
to manage package configuration metadata.
If not using the provided `setup-conda.sh`, you will likely want to configure
meson to install in a local directory. For example, with conda:

``PKG_CONFIG_PATH=$CONDA_PREFIX/.local/lib64/pkgconfig meson
--prefix=$CONDA_PREFIX/.local builddir/``

Note: the default meson build directory is ``builddir/`` (build-release is
used in ``setup-conda.sh``), to disambiguate with the command ``meson build``.
However, many IDEs expect a build directory in ``build/``, as is typical with
``cmake``. It may be convenient to create a symbolic link between them, e.g. by
``ln -s build-release build``. Alternatively, some IDEs may support opening the
built ``compile_commands.json``, which you may also want to symlink to the
root directory.
