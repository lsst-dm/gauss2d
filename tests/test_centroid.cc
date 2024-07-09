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

#include <memory>

#include "lsst/gauss2d/centroid.h"

namespace g2d = lsst::gauss2d;

TEST_CASE("Centroid") {
    double x = 1, y = -1;
    auto values = std::make_shared<g2d::CentroidValues>(x, y);
    CHECK_EQ(*values, g2d::CentroidValues(x, y));
    auto c1 = std::make_shared<g2d::Centroid>(values);
    CHECK_EQ(c1->get_x(), x);
    auto c2 = std::make_shared<g2d::Centroid>(values);
    CHECK_EQ(*c1, *c2);
    c2->set_xy({0, -2});
    CHECK_EQ(c1->get_x(), 0);
    CHECK_EQ(c1->get_y(), -2);

    CHECK_EQ(c1->repr(false),
             "lsst::gauss2d::Centroid(lsst::gauss2d::CentroidValues(0.000000e+00, -2.000000e+00))");
    CHECK_EQ(c1->repr(true, g2d::Object::PY_NAMESPACE_SEPARATOR),
             "lsst.gauss2d.Centroid(data=lsst.gauss2d.CentroidValues(x=0.000000e+00, y=-2.000000e+00))");
}