// -*- LSST-C++ -*-
/*
 * This file is part of gauss2d.
 *
 * Developed for the LSST Data Management System.
 * This product includes software developed by the LSST Project
 * (https://www.lsst.org).
 * See the COPYRIGHT file at the top-level directory of this distribution
 * for details of code ownership.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"

#include "lsst/gauss2d/ellipse.h"

namespace g2d = lsst::gauss2d;

const double RHO_MAX = std::nextafter(1., 2.);
const double RHO_MIN = std::nextafter(-1., -2.);

TEST_CASE("Covariance") {
    CHECK_THROWS_AS(g2d::Covariance(-1), std::invalid_argument);
    CHECK_THROWS_AS(g2d::Covariance(0, -1), std::invalid_argument);
    CHECK_THROWS_AS(g2d::Covariance(0, 0, RHO_MAX), std::invalid_argument);
    CHECK_THROWS_AS(g2d::Covariance(0, 0, RHO_MIN), std::invalid_argument);
    CHECK_EQ(g2d::Covariance(0, 0, 0).get_xyc(), std::array<double, 3>({0., 0., 0.}));
    CHECK_EQ(*(g2d::Covariance(9, 9, 0).make_convolution(g2d::Covariance(16, 16, 0))),
             g2d::Covariance(25, 25, 0));
}

TEST_CASE("Ellipse") {
    CHECK_THROWS_AS(g2d::Ellipse(-1), std::invalid_argument);
    CHECK_THROWS_AS(g2d::Ellipse(0, -1), std::invalid_argument);
    CHECK_THROWS_AS(g2d::Ellipse(0, 0, RHO_MAX), std::invalid_argument);
    CHECK_THROWS_AS(g2d::Ellipse(0, 0, RHO_MIN), std::invalid_argument);
    CHECK_EQ(*(g2d::Ellipse(3, 3, 0).make_convolution(g2d::Ellipse(4, 4, 0))), g2d::Ellipse(5, 5, 0));
    auto ell = g2d::Ellipse(3, 0, 0);
    CHECK_EQ(g2d::Ellipse(g2d::Covariance(ell)), ell);
    auto ellmaj = g2d::EllipseMajor(ell);
    CHECK_EQ(g2d::Ellipse(ellmaj), ell);

    auto ell2 = g2d::Ellipse(1, 1, 0.3);
    CHECK_EQ(ell2.get_xyr(), std::array<double, 3>{1, 1, 0.3});
    CHECK_EQ(ell2.get_hxyr(), std::array<double, 3>{g2d::M_SIGMA_HWHM, g2d::M_SIGMA_HWHM, 0.3});
    ell2.set_hxyr({1, 1, -0.3});
    CHECK_EQ(ell2.get_xyr(), std::array<double, 3>{g2d::M_HWHM_SIGMA, g2d::M_HWHM_SIGMA, -0.3});
}

TEST_CASE("Marginal Major Ellipse") {
    auto covar = g2d::Covariance(0.08333332098858685, 0.08333332098858683, 1.337355953645e-13);
    auto ellipse_maj = g2d::EllipseMajor(covar);
    CHECK_GT(ellipse_maj.get_r_major(), 0);
}
