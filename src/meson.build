gauss2d_sources = [
    'centroid.cc',
    'coordinatesystem.cc',
    'ellipse.cc',
    'evaluate.cc',
    'gaussian.cc',
]
gauss2dlib = shared_library('gauss2d',
                            gauss2d_sources,
                            include_directories : inc,
                            install : true,
                            install_dir : get_option('libdir') / 'lsst',
                            soversion : meson.project_version(),
                           )
