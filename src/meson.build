gauss2d_sources = ['gaussian.cc', 'evaluate.cc']
gauss2dlib = shared_library('gauss2d',
                            gauss2d_sources,
                            include_directories : inc,
                            dependencies : pybind11_dep,
                            install : true,
                            soversion : '0.1.0')
