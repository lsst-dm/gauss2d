export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$CONDA_PREFIX/.local/lib64/pkgconfig
prefix="--prefix=$CONDA_PREFIX/.local"
meson setup $prefix --buildtype release build-release
./build-cc-release.sh
cd python/lib
meson setup $prefix --buildtype release build-release
cd ..
meson setup $prefix --buildtype release build-release
cd ..
./build-py-release.sh

