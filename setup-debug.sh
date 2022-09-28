export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$CONDA_PREFIX/.local/lib64/pkgconfig
meson --prefix=$PWD/debug --buildtype debug build-debug
./build-cc-debug.sh
cd python/lib
meson --prefix=$PWD/debug --buildtype debug build-debug
cd ..
meson --prefix=$PWD/debug --buildtype debug build-debug
cd ..
./build-py-debug.sh

