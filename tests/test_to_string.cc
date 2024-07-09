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

#include <list>
#include <set>
#include <vector>

#include "lsst/gauss2d/to_string.h"

namespace g2d = lsst::gauss2d;

TEST_CASE("to_string_float") {
    CHECK_EQ(g2d::to_string_float(1.0), "1.000000e+00");
    CHECK_EQ(g2d::to_string_float(-1.99999001, 6, false), "-1.999990");
    CHECK_EQ(g2d::to_string_float(42.0, 1, false), "42.0");
    CHECK_EQ(g2d::to_string_float(99.338367, 3), "9.934e+01");
}

TEST_CASE("to_string_float_iter") {
    std::vector<double> t1{1.0, 1.9999901, 42.0, 0.0099338367};
    CHECK_EQ(g2d::to_string_float_iter(t1), "[1.000000e+00, 1.999990e+00, 4.200000e+01, 9.933837e-03]");
}

TEST_CASE("to_string_iter") {
    std::list<int> t1{1, -5, 24};
    CHECK_EQ(g2d::to_string_iter(t1), "[1, -5, 24]");
    std::set<size_t> t2{0, 1, 2, 4};
    CHECK_EQ(g2d::to_string_iter(t2), "[0, 1, 2, 4]");
}
