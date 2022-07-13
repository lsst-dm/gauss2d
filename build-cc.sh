CXXFLAGS="$CXXFLAGS -O3" meson compile -C builddir && meson test -C builddir && meson install -C builddir
