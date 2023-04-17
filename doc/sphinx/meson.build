sphinx = find_program(['sphinx-build'], required: false)

if sphinx.found()
  command = [
    sphinx, '-b', 'html',
    '-Dbreathe_projects.gauss2d=' + meson.build_root() + '/doc/xml',
    meson.current_source_dir() + '/source',
    'doc/sphinx/html'
  ]

  sphinx_target = custom_target('sphinx',
                                output : 'html',
                                command : command,
                                install : true,
                                install_dir : meson.current_build_dir() + '/install')
endif
