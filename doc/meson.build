documenteer = find_program('package-docs', required: false)

if documenteer.found()
    message('Running documenteer to build docs; see ' + meson.build_root() + '/meson-logs/meson-log.txt for output')
    run_command('sh', '-c', 'cd $MESON_SOURCE_ROOT/$MESON_SUBDIR && package-docs build', check: false, capture: true)
endif
