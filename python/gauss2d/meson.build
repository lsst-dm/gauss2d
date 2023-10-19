sources = [
    'centroid.cc',
    'ellipse.cc',
    'gaussian.cc',
    'image.cc',
    'object.cc',
    'pybind11.cc',
]

kwargs = {}
eups = find_program('eups', required: false)
use_eups = eups.found()

if use_eups
	use_eups = run_command('eups', 'list', '-d', 'gauss2d', check: false).returncode() == 0
endif

if use_eups
    # meson.current_build_dir() would be ideal but it seems to delete the file?
    kwargs = {'install_dir': meson.build_root() + '/lib/gauss2d'}
else
    # Folder relative to site-packages to install to
    kwargs = {'subdir': 'gauss2d'}
endif

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
    kwargs: kwargs
)

# Pure Python sources
python_sources = [
    '__init__.py',
    'utils.py',
]

if use_eups
    install_data(python_sources, kwargs: kwargs)
else
    # Install pure Python
    python.install_sources(
        python_sources,
        pure: false, # Will be installed next to binaries
        subdir: 'gauss2d'  # Folder relative to site-packages to install to
    )
endif
