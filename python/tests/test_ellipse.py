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
    ell_1 = g2d.Ellipse(sigma_x=1, sigma_y=1, rho=-0.1)
    assert (ell_1.sigma_x, ell_1.sigma_y, ell_1.rho) == (1, 1, -0.1)
    assert [ell_1.hwhm_x, ell_1.hwhm_y, ell_1.rho] == ell_1.hxyr
    ell_1.set_h(hwhm_x=1, hwhm_y=1, rho=0.1)
    assert (ell_1.hwhm_x, ell_1.hwhm_y, ell_1.rho) == (1, 1, 0.1)
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


def test_EllipseMajor():
    covar = g2d.Covariance(0.08333332098858685, 0.08333332098858683, 1.337355953645e-13)
    ellipse_maj = g2d.EllipseMajor(covar)
    assert ellipse_maj.r_major > 0
