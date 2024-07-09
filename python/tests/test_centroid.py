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


def test_CentroidValues():
    values = g2d.CentroidValues()
    assert values.xy == [0., 0.]

    values.x, values.y = -2., 3.
    assert values.xy == [-2., 3.]

    str_values = "CentroidValues(x=-2.000000e+00, y=3.000000e+00)"
    assert str(values) == str_values
    assert repr(values) == f"lsst.gauss2d.{str_values}"


def test_Centroid():
    values = g2d.CentroidValues(-1.0, 2.0)
    cen1 = g2d.Centroid(values)
    cen2 = g2d.Centroid(values)
    cen3 = g2d.Centroid(-1.0, 2.0)

    assert cen1 == cen2
    assert not (cen1 != cen3)

    cen1.x = cen1.y
    assert cen2.x == cen1.x
    assert cen1 != cen3

    str_values = str(values)
    assert str(cen1) == f"Centroid(data={str_values})"
    assert repr(cen1) == f"lsst.gauss2d.Centroid(data=lsst.gauss2d.{str_values})"
