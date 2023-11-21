# Compile, test and install C++ library
meson compile -C build-release && meson test -C build-release && meson install -C build-release

