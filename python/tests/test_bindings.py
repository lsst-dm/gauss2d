import lsst.gauss2d as g2


def test_bindings():
    errors = []
    bad_args_msg = "incompatible function arguments"

    for name_class in dir(g2):
        klass = getattr(g2, name_class)

        if callable(klass):
            try:
                obj = klass()
            except ValueError as e:
                # Likely non-trivial __init__
                continue
            except TypeError as e:
                msg = str(e)
                if not (
                    msg.endswith("No constructor defined!")
                    or msg.startswith(f"__init__(): incompatible constructor arguments.")
                    or msg.startswith("extend_path() missing")
                    or (msg.startswith("make_gaussians_pixel") and (bad_args_msg in msg))
                ):
                    errors.append(e)
                    # str and repr should not raise
            try:
                repr(obj)
            except Exception as e:
                errors.append(e)
            try:
                str(obj)
            except Exception as e:
                errors.append(e)

    if errors:
        print('\n'.join((str(e) for e in errors)))
        assert not 'errors not empty; see stdout'
