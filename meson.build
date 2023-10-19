project(
  'libgauss2d',
  'cpp',
  version : run_command('cat', meson.source_root() + '/VERSION', check: true).stdout().strip(),
  license : 'GPL-3.0-or-later',
  default_options : ['cpp_std=c++17',],
)

eups = find_program('eups', required: false)
use_eups = eups.found()
if use_eups
    use_eups = run_command('eups', 'list', 'gauss2d', check: false).returncode() == 0
endif

inc = include_directories('include')

subdir('include')
subdir('src')
subdir('tests')

if not use_eups
    pkg_mod = import('pkgconfig', required: false)
    if pkg_mod.found()
        pkg_mod.generate(
            libraries : gauss2dlib,
            version : meson.project_version(),
            name : 'libgauss2d',
            filebase : 'gauss2d',
            description : 'Create, manipulate and evaluate 2D Gaussian mixtures and images thereof.'
        )
    endif
endif

subdir('doc')
