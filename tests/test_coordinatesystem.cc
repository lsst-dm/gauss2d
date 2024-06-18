#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"

#include <limits>
#include <memory>
#include <stdexcept>

#include "lsst/gauss2d/coordinatesystem.h"

namespace g2d = lsst::gauss2d;

TEST_CASE("CoordinateSystem") {
    auto inf = std::numeric_limits<double>::infinity();
    auto nan = std::numeric_limits<double>::quiet_NaN();
    CHECK_THROWS_AS(g2d::CoordinateSystem{inf}, std::invalid_argument);
    CHECK_THROWS_AS(g2d::CoordinateSystem{0.}, std::invalid_argument);
    CHECK_THROWS_AS(g2d::CoordinateSystem(5.0, nan), std::invalid_argument);
    CHECK_THROWS_AS(g2d::CoordinateSystem(1.0, 1.0, -0.1, -inf), std::invalid_argument);
    CHECK_THROWS_AS(g2d::CoordinateSystem(1.0, 1.0, nan, 1.0), std::invalid_argument);
    double dx1 = 0.1, dy2 = 0.8, x_min = 18.0, y_min = -5.5;
    auto coord = g2d::CoordinateSystem(dx1, dy2, x_min, y_min);
    CHECK_EQ(dx1, coord.get_dx1());
    CHECK_EQ(dy2, coord.get_dy2());
    CHECK_EQ(x_min, coord.get_x_min());
    CHECK_EQ(y_min, coord.get_y_min());
    CHECK_EQ(true, coord.is_xy_aligned());
    CHECK_EQ(coord, g2d::CoordinateSystem(dx1, dy2, x_min, y_min));
    CHECK_NE(coord, g2d::CoordinateSystem(dx1 + 1e-5, dy2, x_min, y_min));
    CHECK_NE(coord, g2d::CoordinateSystem(dx1, dy2 + 1e-5, x_min, y_min));
    CHECK_NE(coord, g2d::CoordinateSystem(dx1, dy2, x_min - 1e-5, y_min));
    CHECK_NE(coord, g2d::CoordinateSystem(dx1, dy2, x_min, y_min - 1e-5));

    CHECK_EQ(coord.repr(),
             "lsst::gauss2d::CoordinateSystem(1.000000e-01, 8.000000e-01, 1.800000e+01, -5.500000e+00)");
}