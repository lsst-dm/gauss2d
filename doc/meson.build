documenteer = find_program('package-docs', required: false)

if documenteer.found()
    run_command('sh', '-c', 'cd $MESON_SOURCE_ROOT/$MESON_SUBDIR && package-docs build', check: false, capture: true)
endif
