project(
  'libgauss2d_python',
  'cpp',
  default_options : ['cpp_std=c++17',],
  license : 'GPL-3.0-or-later',
  version : run_command('cat', meson.source_root() + '/VERSION', check: true).stdout().strip(),
)

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
        dependencies: [cpp.find_library('gauss2d', dirs: [eupsdir / 'lib64/lsst']),],
        include_directories: eupsdir / 'include',
    )
else
    gauss2d_dep = dependency('lsst_gauss2d')
endif

py = import('python').find_installation('python3', modules : ['numpy'])
pybind11_dep = dependency('pybind11')

inc = include_directories('include')

subdir('include')
subdir('tests')

if not use_eups
    pkg_mod = import('pkgconfig')
    pkg_mod.generate(version : meson.project_version(),
                     name : 'libgauss2d_python',
                     filebase : 'lsst_gauss2d_python',
                     description : 'Python-specific headers for gauss2d.')
endif
