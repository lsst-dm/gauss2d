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

The shared library can be built using Meson as follows:

meson --prefix=~/.local build && meson compile -C build && meson install -C build

... where the prefix argument is your desired installation path.

You may also need to configure pkg-config, e.g. prepend 
PKG_CONFIG_PATH=~/.local/lib64/pkgconfig to the build command for a local 
and/or non-root installation.

