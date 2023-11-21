# Configure package and build debug libraries
meson --prefix=$PWD/debug --buildtype debug build-debug
./build-cc-debug.sh
cd python/lib
meson --prefix=$PWD/debug --buildtype debug build-debug
cd ..
meson --prefix=$PWD/debug --buildtype debug build-debug
cd ..
./build-py-debug.sh

