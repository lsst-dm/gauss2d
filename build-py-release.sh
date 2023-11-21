# Compile, test and install python release libraries
# Unfortunately, the library must be installed before any python tests are run
# (otherwise the import will fail)
meson compile -C python/lib/build-release && meson install -C python/lib/build-release && meson test -C python/lib/build-release
meson compile -C python/build-release && meson install -C python/build-release && meson test -C python/build-release

