# Compile python debug libraries
meson compile -C python/lib/build-debug && meson install -C python/lib/build-debug && meson test -C python/lib/build-debug
meson compile -C python/build-debug && meson install -C python/build-debug && meson test -C python/build-debug
