project(
    'gauss2d',
    'cpp',
    default_options : ['cpp_std=c++17',],
    version : run_command('cat', meson.source_root() + '/VERSION', check: true).stdout().strip(),
)

add_project_arguments('-DVERSION=@0@'.format(meson.project_version()), language : 'cpp')

eups = find_program('eups', required: false)
use_eups = eups.found()

if use_eups
	eupsdir = run_command('eups', 'list', '-d', 'gauss2d', check: false)
	use_eups = eupsdir.returncode() == 0
endif

if use_eups
    eupsdir = eupsdir.stdout().strip()
    cpp = meson.get_compiler('cpp')
    gauss2d_dep = declare_dependency(
        dependencies: [cpp.find_library('gauss2d', dirs: [eupsdir + '/build-release/lib64']),],
        include_directories: eupsdir + '/build-release/include',
    )
    gauss2d_python_dep = declare_dependency(
        # Relative path because it's in a subdir of this file
        include_directories: 'lib/build-release/include',
    )
else
    gauss2d_dep = dependency('gauss2d')
    gauss2d_python_dep = dependency('gauss2d_python')
endif

pybind11_dep = dependency('pybind11')
pymod = import('python')
python = pymod.find_installation(
    'python3',
    required : true,
)
subdir('gauss2d')

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
