project('libgauss2d', 'cpp',
  version : '0.1.0',
  license : 'GPL-3.0-or-later')

inc = include_directories('include')

subdir('include')
subdir('src')
#subdir('test')

pkg_mod = import('pkgconfig')
pkg_mod.generate(libraries : gauss2dlib,
                 version : '0.1',
                 name : 'libgauss2d',
                 filebase : 'gauss2d',
                 description : 'Create, manipulate and evaluate 2D Gaussian mixtures and images thereof.')
