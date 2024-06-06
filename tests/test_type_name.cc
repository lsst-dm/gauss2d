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
