# Configure package and build release libraries
meson setup "$@" --buildtype release build-release
./build-cc-release.sh
cd python/lib
meson setup "$@" --buildtype release build-release
cd ..
meson setup "$@" --buildtype release build-release
cd ..
./build-py-release.sh
