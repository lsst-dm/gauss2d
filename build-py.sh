CXXFLAGS="$CXXFLAGS -O3" meson compile -C python/lib/builddir && meson install -C python/lib/builddir
CXXFLAGS="$CXXFLAGS -O3" meson compile -C python/builddir && meson install -C python/builddir
