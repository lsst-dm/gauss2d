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

#include "lsst/gauss2d/string_utils.h"

namespace g2d = lsst::gauss2d;

TEST_CASE("replace_all_str_str") {
    CHECK_EQ(g2d::replace_all("testWHOwhatWHOSE", std::string("WHO"), std::string("WHERE")),
             "testWHEREwhatWHERESE");
}

TEST_CASE("replace_all_mixed") {
    std::string_view token = "3x";
    CHECK_EQ(g2d::replace_all("ABC123x3xyz", token, std::string("")), "ABC12yz");
    CHECK_EQ(g2d::replace_all("ABC123x3xyz", std::string("x3"), token), "ABC1233xxyz");
}

TEST_CASE("replace_all_none") {
    std::string_view token = "3x";
    CHECK_EQ(g2d::replace_all_none("ABC123x3xyz", token), "ABC12yz");
    CHECK_EQ(g2d::replace_all_none("ABC123x3xyz", std::string("3")), "ABC12xxyz");
}
