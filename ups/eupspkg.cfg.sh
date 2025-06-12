# EupsPkg config file. Sourced by 'eupspkg'

build_cc()
{
    ./build-cc-release.sh
}

build_py()
{
    ./build-py-release.sh
}

build()
{
    (build_cc && build_py)
}

clean()
{
    ./clean.sh
}

config_cc()
{
    ([ -d "$GAUSS2D_DIR" ] && ./clean-cc.sh && meson setup --prefix="$GAUSS2D_DIR/build-release" \
     --buildtype release --libdir "lib" build-release)
}

config_py()
{
    ([ -d "$GAUSS2D_DIR" ] && ./clean-py.sh \
     && cd python/lib \
     && meson setup --prefix="$GAUSS2D_DIR/python/lib/build-release" --buildtype release build-release \
     && cd .. \
     && meson setup --prefix="$GAUSS2D_DIR/python/build-release" --buildtype release build-release)
}

config()
{
    config_cc && build_cc && config_py
}
