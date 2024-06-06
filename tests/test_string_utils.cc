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
