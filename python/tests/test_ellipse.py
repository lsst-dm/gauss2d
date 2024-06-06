import lsst.gauss2d as g2d

import math
import pytest

rho_min = math.nextafter(-1, -2)
rho_max = math.nextafter(1, 2)
pos_min = math.nextafter(0, 1)
prefix_namespace = "lsst.gauss2d."


def test_Covariance():
    with pytest.raises(ValueError):
        g2d.Covariance(-1)
    with pytest.raises(ValueError):
        g2d.Covariance(0, -1)
    for rho_bad in (rho_min, rho_max):
        with pytest.raises(ValueError):
            g2d.Covariance(0, 0, rho_bad)

    covar_0 = g2d.Covariance()
    assert (covar_0.sigma_x_sq, covar_0.sigma_y_sq, covar_0.cov_xy) == (0, 0, 0)
    assert covar_0.xyc == [0, 0, 0]
    assert covar_0 != g2d.Covariance(pos_min, 0, 0)
    assert covar_0 != g2d.Covariance(0, pos_min, 0)
    assert g2d.Covariance(1, 1, 0) != g2d.Covariance(1, 1, pos_min)

    covar_conv = g2d.Covariance(9., 9., 0).make_convolution(g2d.Covariance(16., 16., 0))
    assert covar_conv == g2d.Covariance(25., 25., 0)

    str_covar_conv = "Covariance(sigma_x_sq=2.500000e+01, sigma_y_sq=2.500000e+01, cov_xy=0.000000e+00)"
    assert str(str_covar_conv) == str_covar_conv
    assert repr(covar_conv) == f"{prefix_namespace}{str_covar_conv}"


def test_Ellipse():
    with pytest.raises(ValueError):
        g2d.Ellipse(-1)
    with pytest.raises(ValueError):
        g2d.Ellipse(0, -1)
    for rho_bad in (rho_min, rho_max):
        with pytest.raises(ValueError):
            g2d.Ellipse(0, 0, rho_bad)

    ell_0 = g2d.Ellipse()
    assert ell_0.xyr == [0, 0, 0]
    assert ell_0 != g2d.Ellipse(pos_min, 0, 0)
    assert ell_0 != g2d.Ellipse(0, pos_min, 0)
    ell_1 = g2d.Ellipse(1, 1, 0)
    assert ell_1 != g2d.Ellipse(1, 1, pos_min)

    ell_conv = g2d.Ellipse(3., 3., 0).make_convolution(g2d.Ellipse(4., 4., 0))
    assert ell_conv == g2d.Ellipse(5., 5., 0)
    assert ell_conv.get_radius_trace() == pytest.approx(5*math.sqrt(2.), rel=1e-10, abs=1e-10)

    str_data = "EllipseValues(sigma_x=5.000000e+00, sigma_y=5.000000e+00, rho=0.000000e+00)"
    assert str(ell_conv) == f"Ellipse(data={str_data})"
    assert repr(ell_conv) == f"{prefix_namespace}Ellipse(data={prefix_namespace}{str_data})"

    ell_conv.set(g2d.Covariance(ell_0))
    print(g2d.Covariance(ell_0))
    assert ell_conv == ell_0
    ell_1_maj = g2d.EllipseMajor(ell_1)
    ell_conv.set(ell_1_maj)
    assert ell_conv == g2d.Ellipse(ell_1_maj)

"""    
    .def("set_h", &gauss2d::Ellipse::set_h)
    .def_property("rho", &gauss2d::Ellipse::get_rho, &gauss2d::Ellipse::set_rho)
    .def_property("hwhm_x", &gauss2d::Ellipse::get_hwhm_x, &gauss2d::Ellipse::set_hwhm_x)
    .def_property("hwhm_y", &gauss2d::Ellipse::get_hwhm_y, &gauss2d::Ellipse::set_hwhm_y)
    .def_property("sigma_x", &gauss2d::Ellipse::get_sigma_x, &gauss2d::Ellipse::set_sigma_x)
    .def_property("sigma_y", &gauss2d::Ellipse::get_sigma_y, &gauss2d::Ellipse::set_sigma_y)
    .def_property("hxyr", &gauss2d::Ellipse::get_hxyr, &gauss2d::Ellipse::set_hxyr)
    .def_property("xyr", &gauss2d::Ellipse::get_xyr, &gauss2d::Ellipse::set_xyr)
"""

def test_EllipseMajor():
    covar = g2d.Covariance(0.08333332098858685, 0.08333332098858683, 1.337355953645e-13)
    ellipse_maj = g2d.EllipseMajor(covar)
    assert ellipse_maj.r_major > 0
