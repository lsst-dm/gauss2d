#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"

#include <cmath>
#include <memory>

#include "ellipse.h"

namespace g2 = gauss2d;

TEST_CASE("Ellipse") {
    double rho_max = std::nextafter(1., 2.);
    double rho_min = std::nextafter(-1., -2.);
    CHECK_THROWS_AS(g2::Covariance(-1), std::invalid_argument);
    CHECK_THROWS_AS(g2::Covariance(0, -1), std::invalid_argument);
    CHECK_THROWS_AS(g2::Covariance(0, 0, rho_max), std::invalid_argument);
    CHECK_THROWS_AS(g2::Covariance(0, 0, rho_min), std::invalid_argument);
    CHECK(g2::Covariance(0, 0, 0).get_xyc() == std::array<double, 3>({0., 0., 0.}));
    CHECK(*(g2::Covariance(9, 9, 0).make_convolution(g2::Covariance(16, 16, 0)))
          == g2::Covariance(25, 25, 0));

    CHECK_THROWS_AS(g2::Ellipse(-1), std::invalid_argument);
    CHECK_THROWS_AS(g2::Ellipse(0, -1), std::invalid_argument);
    CHECK_THROWS_AS(g2::Ellipse(0, 0, rho_max), std::invalid_argument);
    CHECK_THROWS_AS(g2::Ellipse(0, 0, rho_min), std::invalid_argument);
    CHECK(*(g2::Ellipse(3, 3, 0).make_convolution(g2::Ellipse(4, 4, 0))) == g2::Ellipse(5, 5, 0));
    auto ell = g2::Ellipse(3, 0, 0);
    CHECK(g2::Ellipse(g2::Covariance(ell)) == ell);
    auto ellmaj = g2::EllipseMajor(ell);
    std::cout << ellmaj << std::endl;
    CHECK(g2::Ellipse(ellmaj) == ell);

    auto ell2 = g2::Ellipse(1, 1, 0.3);
    CHECK(ell2.get_xyr() == std::array<double, 3>{1, 1, 0.3});
    CHECK(ell2.get_hxyr() == std::array<double, 3>{g2::M_SIGMA_HWHM, g2::M_SIGMA_HWHM, 0.3});
    ell2.set_hxyr({1, 1, -0.3});
    CHECK(ell2.get_xyr() == std::array<double, 3>{g2::M_HWHM_SIGMA, g2::M_HWHM_SIGMA, -0.3});
}
