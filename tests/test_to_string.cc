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
