# This file is part of gauss2d.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

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
