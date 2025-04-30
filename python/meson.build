project(
    'gauss2d_py',
    'cpp',
    default_options : ['cpp_std=c++20',],
    version : files('VERSION'),
)

add_project_arguments('-DVERSION=@0@'.format(meson.project_version()), language : 'cpp')

eups = find_program('eups', required: false)
use_eups = eups.found()

if use_eups
	eupsdir = run_command('sh', '-c', 'echo $GAUSS2D_DIR', check: false, capture: true)
	use_eups = eupsdir.returncode() == 0
endif

if use_eups
    eupsdir = eupsdir.stdout().strip() / 'build-release/'
    message('eupsdir=' + eupsdir)
    cpp = meson.get_compiler('cpp')
    gauss2d_dep = declare_dependency(
        dependencies: [cpp.find_library('gauss2d', dirs: [eupsdir / 'lib/lsst']),],
        include_directories: eupsdir / 'include',
    )
    gauss2d_python_dep = declare_dependency(
        # Relative path because it's in a subdir of this file
        include_directories: 'lib/build-release/include',
    )
else
    gauss2d_dep = dependency('lsst_gauss2d')
    gauss2d_python_dep = dependency('lsst_gauss2d_python')
endif

pybind11_dep = dependency('pybind11')
pymod = import('python')
python = pymod.find_installation(
    'python3',
    required : true,
)
subdir('src')

if not use_eups
    # A hacky way to pip install (to say the least) the Python module
    # this is "necessary" to allow its use as a dependency in other packages
    # (e.g. in requirements.txt, pyproject.toml)
    # because meson doesn't generate a dist-info
    # the downside is there are now two ways to manage/install the package
    # pip doesn't seem to find the pybind extension module automagically either
    # (perhaps some config is missing in pyproject.toml)
    subdir_module = meson.source_root() / 'src/lsst/gauss2d'
	pip_cmd = [
	    'python', '-m', 'pip', 'install', meson.source_root(),
	     '-t', python.get_install_dir(),
	     '--no-dependencies',
	     '--no-build-isolation',
	     '--force-reinstall',
     ]
	message('Will run on install:')
	message(pip_cmd)
	add_link_cmd = ['ln', '-fs', module.full_path(), subdir_module]
	meson.add_install_script(add_link_cmd)
	meson.add_install_script(pip_cmd)
	remove_link_cmd = ['rm', '-f', subdir_module / module.full_path().split('/')[-1]]
	meson.add_install_script(remove_link_cmd)
endif

# Run tests - not sure why workdir is required, but it is
test('pytest', python, args : ['-m', 'pytest', 'tests'], workdir : meson.project_source_root())

# Build package docs (optional)
subdir('doc')
