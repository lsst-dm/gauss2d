sources = [
    'centroid.cc',
    'ellipse.cc',
    'gaussian.cc',
    'image.cc',
    'object.cc',
    'pybind11.cc',
]
module = python.extension_module(
    '_gauss2d', # Name of the module
    sources,  # List of sources
    install: true,  # Install it
    dependencies: [ # List of dependencies
        python.dependency(), # Add Python.h as dependency
        gauss2d_dep,
        gauss2d_python_dep,
        pybind11_dep,
    ],
    subdir: 'gauss2d', # Folder relative to site-packages to install to
)

# Pure Python sources
python_sources = [
    '__init__.py',
    'utils.py',
]

# Install pure Python
python.install_sources(
    python_sources,
    pure: false, # Will be installed next to binaries
    subdir: 'gauss2d'  # Folder relative to site-packages to install to
)