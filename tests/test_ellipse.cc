#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"

#include <memory>

#include "ellipse.h"

namespace g2 = gauss2d;

TEST_CASE("Ellipse")
{
    CHECK_THROWS_AS(g2::Covariance(-1), std::invalid_argument);
    CHECK_THROWS_AS(g2::Covariance(0, -1), std::invalid_argument);
    CHECK_THROWS_AS(g2::Covariance(0, 0, 1), std::invalid_argument);
    CHECK(g2::Covariance(0, 0, 0).get_xyc() == std::array<double, 3>({0., 0., 0.}));
    CHECK(*(g2::Covariance(9, 9, 0).make_convolution(g2::Covariance(16, 16, 0)))
        == g2::Covariance(25, 25, 0));

    CHECK_THROWS_AS(g2::Ellipse(-1), std::invalid_argument);
    CHECK_THROWS_AS(g2::Ellipse(0, -1), std::invalid_argument);
    CHECK_THROWS_AS(g2::Ellipse(0, 0, 1), std::invalid_argument);
    CHECK(*(g2::Ellipse(3, 3, 0).make_convolution(g2::Ellipse(4, 4, 0)))
        == g2::Ellipse(5, 5, 0));
}