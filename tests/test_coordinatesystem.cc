#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"

#include <limits>
#include <memory>
#include <stdexcept>

#include "lsst/gauss2d/coordinatesystem.h"

namespace g2 = lsst::gauss2d;

TEST_CASE("Centroid") {
    CHECK_THROWS_AS(g2::CoordinateSystem(std::numeric_limits<double>::infinity()), std::invalid_argument);
    CHECK_THROWS_AS(g2::CoordinateSystem(0.), std::invalid_argument);
    CHECK_THROWS_AS(g2::CoordinateSystem(5.0, std::numeric_limits<double>::quiet_NaN()),
                    std::invalid_argument);
    double dx1 = 0.1, dy2 = 0.8, x_min = 18.0, y_min = -5.5;
    auto coord = g2::CoordinateSystem(dx1, dy2, x_min, y_min);
    CHECK_EQ(dx1, coord.get_dx1());
    CHECK_EQ(dy2, coord.get_dy2());
    CHECK_EQ(x_min, coord.get_x_min());
    CHECK_EQ(y_min, coord.get_y_min());
    CHECK_EQ(true, coord.is_xy_aligned());
    CHECK_EQ(coord, g2::CoordinateSystem(dx1, dy2, x_min, y_min));
    CHECK_NE(coord, g2::CoordinateSystem(dx1 + 1e-5, dy2, x_min, y_min));
    CHECK_NE(coord, g2::CoordinateSystem(dx1, dy2 + 1e-5, x_min, y_min));
    CHECK_NE(coord, g2::CoordinateSystem(dx1, dy2, x_min - 1e-5, y_min));
    CHECK_NE(coord, g2::CoordinateSystem(dx1, dy2, x_min, y_min - 1e-5));
}