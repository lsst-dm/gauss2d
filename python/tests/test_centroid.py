import lsst.gauss2d as g2


def test_CentroidValues():
    values = g2.CentroidValues()
    assert values.xy == [0., 0.]

    values.x, values.y = -2., 3.
    assert values.xy == [-2., 3.]

    assert str(values) == "CentroidValues(x=-2.000000, y=3.000000)"


def test_Centroid():
    values = g2.CentroidValues(-1.0, 2.0)
    cen1 = g2.Centroid(values)
    cen2 = g2.Centroid(values)
    cen3 = g2.Centroid(-1.0, 2.0)

    assert cen1 == cen2
    assert not (cen1 != cen3)

    cen1.x = cen1.y
    assert cen2.x == cen1.x
    assert cen1 != cen3
