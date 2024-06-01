import lsst.gauss2d as g2

import math
import pytest


def test_CoordinateSystem():
    for value_bad in (0., math.inf, math.nan):
        with pytest.raises(ValueError):
            g2.CoordinateSystem(value_bad)
        with pytest.raises(ValueError):
            g2.CoordinateSystem(1.0, value_bad)

    for value_bad in (math.inf, math.nan):
        with pytest.raises(ValueError):
            print(value_bad)
            g2.CoordinateSystem(1.0, 1.0, value_bad)
            g2.CoordinateSystem(1.0, 1.0, 1.0, value_bad)

    coordsys1 = g2.CoordinateSystem(1.0, 3.0, -10., 23.)
    coordsys2 = g2.CoordinateSystem(1.0, 3.0, -10., 23.)

    assert coordsys1 == coordsys2
    assert str(coordsys1) == "CoordinateSystem(dx1=1.000000, dy2=3.000000, x_min=-10.000000, y_min=23.000000)"
