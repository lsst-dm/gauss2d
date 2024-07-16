# Compile and test C++ debug library
meson compile -C build-debug && meson test -C build-debug && meson install -C build-debug
