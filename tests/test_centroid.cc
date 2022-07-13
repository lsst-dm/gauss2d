#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"

#include <memory>

#include "centroid.h"

namespace g2 = gauss2d;

TEST_CASE("Centroid")
{
    double x = 1, y = -1;
    auto values = std::make_shared<g2::CentroidValues>(x, y);
    CHECK(*values == g2::CentroidValues(x, y));
    auto c1 = std::make_shared<g2::Centroid>(values);
    CHECK(c1->get_x() == x);
    auto c2 = std::make_shared<g2::Centroid>(values);
    CHECK(*c1 == *c2);
    c2->set_xy({0, -2});
    CHECK(c1->get_x() == 0);
    CHECK(c1->get_y() == -2);
}