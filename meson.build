project(
  'libgauss2d',
  'cpp',
  version : files('VERSION'),
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

# Unit Tests
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

# Build package docs (optional)
subdir('doc')

# Build doxygen docs (done by eupspkg install if use_eups)
if not use_eups
    pymod = import('python')
    python = pymod.find_installation(
        'python3',
        required : false,
        modules: ['lsst.sconsUtils'],
    )
    has_py3 = python.found()
    scons = find_program('scons', required: false)
    has_scons = scons.found()
    if has_scons and has_py3
        message('Running SCons to build docs; see ' + meson.build_root() + '/meson-logs/meson-log.txt for output')
        run_command('sh', '-c', 'cd $MESON_SOURCE_ROOT && echo $PWD && scons', check: false, capture: true)
    endif
endif
