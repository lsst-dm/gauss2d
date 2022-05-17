export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$CONDA_PREFIX/.local/lib64/pkgconfig
prefix="--prefix=$CONDA_PREFIX/.local"
meson $prefix builddir
sh build-cc.sh
meson $prefix python/lib/builddir
CXXFLAGS="$CXXFLAGS -O3" meson compile -C python/lib/builddir && meson install -C python/lib/builddir
meson $prefix python/builddir
CXXFLAGS="$CXXFLAGS -O3" meson compile -C python/builddir && meson install -C python/builddir

