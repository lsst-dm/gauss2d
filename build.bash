CXXFLAGS="$CXXFLAGS -O3" meson compile -C build && meson test -C build && meson install -C build && meson compile -C python/build && meson install -C python/build
