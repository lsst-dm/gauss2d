subdir_gauss2d = 'lsst' / 'gauss2d'
sources = [
    subdir_gauss2d / 'centroid.cc',
    subdir_gauss2d / 'coordinatesystem.cc',
    subdir_gauss2d / 'ellipse.cc',
    subdir_gauss2d / 'gaussian.cc',
    subdir_gauss2d / 'image.cc',
    subdir_gauss2d / 'object.cc',
    subdir_gauss2d / 'pybind11.cc',
]

add_project_arguments('-DVERSION=@0@'.format(meson.project_version()), language : 'cpp')

eups = find_program('eups', required: false)
use_eups = eups.found()

if use_eups
	use_eups = run_command('eups', 'list', '-d', 'gauss2d', check: false).returncode() == 0
endif

kwargs = {}
if use_eups
    # meson.current_build_dir() would be ideal but it seems to delete the file?
    kwargs = {'install_dir': meson.build_root() / 'lib' / subdir_gauss2d}
else
    # Folder relative to site-packages to install to
    kwargs = {'subdir': subdir_gauss2d}
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
py_sources = [
    'lsst' / '__init__.py',
    subdir_gauss2d / '__init__.py',
    subdir_gauss2d / 'utils.py',
]

if use_eups
    install_data(py_sources, kwargs: kwargs)
else
    # Install pure Python
    python.install_sources(
        py_sources,
        pure: false, # Will be installed next to binaries
        preserve_path: true,
    )
endif