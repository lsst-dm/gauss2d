doxygen = find_program(['doxygen'], required: false)

if doxygen.found()
    build = meson.current_build_dir()

    conf_data = configuration_data()
    conf_data.set('DOXYGEN_OUTPUT_DIR', build)
    conf_data.set('DOXYGEN_INPUT_DIR',  meson.source_root() + '/include')


    doxygen_file = configure_file(output: 'doxyfile',
                            input: 'doxyfile.in',
                            configuration: conf_data,
                            install: false)

    # TODO: Determine why install: true is needed
    # Right now this builds the docs in build + '/doc'
    # Then "installs" them in build + '/doc/build'
    # With install: false it does not build them at all
    # (or perhaps does and then deletes them?)
    doxygen_target = custom_target('doxygen',
                                input: doxygen_file,
                                output: 'xml',
                                command: [doxygen, doxygen_file],
                                install: true,
                                install_dir: build + '/build-release',
                                )

    subdir('sphinx')
endif
