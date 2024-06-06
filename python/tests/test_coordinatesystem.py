import lsst.gauss2d as g2d

import math
import pytest


def test_CoordinateSystem():
    for value_bad in (0., math.inf, math.nan):
        with pytest.raises(ValueError):
            g2d.CoordinateSystem(value_bad)
        with pytest.raises(ValueError):
            g2d.CoordinateSystem(1.0, value_bad)

    for value_bad in (math.inf, math.nan):
        with pytest.raises(ValueError):
            print(value_bad)
            g2d.CoordinateSystem(1.0, 1.0, value_bad)
            g2d.CoordinateSystem(1.0, 1.0, 1.0, value_bad)

    coordsys1 = g2d.CoordinateSystem(1.0, 3.0, -10., 23.)
    coordsys2 = g2d.CoordinateSystem(1.0, 3.0, -10., 23.)

    assert coordsys1 == coordsys2
    str_cs1 = "CoordinateSystem(dx1=1.000000e+00, dy2=3.000000e+00, x_min=-1.000000e+01, y_min=2.300000e+01)"
    assert str(coordsys1) == str_cs1
    assert repr(coordsys1) == "lsst.gauss2d." + str_cs1
