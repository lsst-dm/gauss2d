#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"

#include <memory>

#include "lsst/gauss2d/centroid.h"

namespace g2d = lsst::gauss2d;

TEST_CASE("Centroid") {
    double x = 1, y = -1;
    auto values = std::make_shared<g2d::CentroidValues>(x, y);
    CHECK_EQ(*values, g2d::CentroidValues(x, y));
    auto c1 = std::make_shared<g2d::Centroid>(values);
    CHECK_EQ(c1->get_x(), x);
    auto c2 = std::make_shared<g2d::Centroid>(values);
    CHECK_EQ(*c1, *c2);
    c2->set_xy({0, -2});
    CHECK_EQ(c1->get_x(), 0);
    CHECK_EQ(c1->get_y(), -2);

    CHECK_EQ(c1->repr(false),
             "lsst::gauss2d::Centroid(lsst::gauss2d::CentroidValues(0.000000e+00, -2.000000e+00))");
    CHECK_EQ(c1->repr(true, g2d::Object::PY_NAMESPACE_SEPARATOR),
             "lsst.gauss2d.Centroid(data=lsst.gauss2d.CentroidValues(x=0.000000e+00, y=-2.000000e+00))");
}