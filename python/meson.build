project(
    'gauss2d_py',
    'cpp',
    default_options : ['cpp_std=c++17',],
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
        dependencies: [cpp.find_library('gauss2d', dirs: [eupsdir + 'lib64/lsst/gauss2d']),],
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
subdir('lsst')
# Run tests - not sure why workdir is required, but it is
test('pytest', python, args : ['-m', 'pytest', 'tests'], workdir : meson.source_root())

if not use_eups
	pip_cmd = [
	    'python', '-m', 'pip', 'install', meson.source_root(),
	     '-t', python.get_install_dir(),
	     '--no-dependencies',
	     '--no-build-isolation',
	     '--upgrade',
     ]
	message('Will run on install:')
	message(pip_cmd)
	meson.add_install_script(pip_cmd)
endif
