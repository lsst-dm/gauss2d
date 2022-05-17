project(
  'libgauss2d_python',
  'cpp',
  version : '0.1.0',
  license : 'GPL-3.0-or-later',
  default_options : ['cpp_std=c++17',],
)

gauss2d_dep = dependency('gauss2d')
pybind11_dep = dependency('pybind11')

inc = include_directories('include')

depends = declare_dependency(
  include_directories: inc,
  dependencies: [gauss2d_dep, pybind11_dep],
)

subdir('include')

pkg_mod = import('pkgconfig')
pkg_mod.generate(version : '0.1',
                 name : 'libgauss2d_python',
                 filebase : 'gauss2d_python',
                 description : 'Python-specific headers for gauss2d.')
