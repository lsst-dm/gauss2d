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

#include "lsst/gauss2d/type_name.h"

namespace g2d = lsst::gauss2d;

namespace lsst {
class TestLsst {
public:
    double y;
};
namespace gauss2d {
class Test {
public:
    int x;
};
}  // namespace gauss2d
}  // namespace lsst

TEST_CASE("type_name: Test") {
    CHECK_EQ(g2d::type_name_str<lsst::TestLsst>(false), "lsst::TestLsst");
    CHECK_EQ(g2d::type_name_str<g2d::Test>(false), "lsst::gauss2d::Test");
    CHECK_EQ(g2d::type_name_str<g2d::Test>(true), "Test");
    CHECK_EQ(g2d::type_name_str<g2d::Test>(false, "."), "lsst.gauss2d.Test");
    CHECK_EQ(g2d::type_name_str<double>(false), "double");
}
