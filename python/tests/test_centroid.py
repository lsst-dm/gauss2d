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
