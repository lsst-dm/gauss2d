# Configure package to install to a local conda folder 
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$CONDA_PREFIX/.local/lib64/pkgconfig

mkdir -p $CONDA_PREFIX/etc/conda/activate.d/
mkdir -p $CONDA_PREFIX/etc/conda/deactivate.d/

# To automatically add the C++ library paths to $LD_LIBRARY_PATH, see:
# https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#macos-and-linux
# (adapted from https://stackoverflow.com/a/49238956)

if [ ! -f $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh ]; then
	echo "export LD_LIBRARY_PATH_OLD=\${LD_LIBRARY_PATH};export LD_LIBRARY_PATH=\${CONDA_PREFIX}/.local/lib64:\${LD_LIBRARY_PATH}" | sed "s/;/\n/g" > $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh
fi
if [ ! -f $CONDA_PREFIX/etc/conda/deactivate.d/env_vars.sh ]; then
	echo "export LD_LIBRARY_PATH=\${LD_LIBRARY_PATH_OLD};unset LD_LIBRARY_PATH_OLD" | sed "s/;/\n/g" > $CONDA_PREFIX/etc/conda/deactivate.d/env_vars.sh
fi

./setup-release.sh "--prefix=$CONDA_PREFIX/.local"

