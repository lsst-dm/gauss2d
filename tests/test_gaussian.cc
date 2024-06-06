#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"

#include <memory>
#include <stdexcept>

#include "lsst/gauss2d/centroid.h"
#include "lsst/gauss2d/gaussian.h"

namespace g2d = lsst::gauss2d;

TEST_CASE("Gaussian") {
    auto g1 = g2d::Gaussian{};
    auto g2 = std::make_shared<g2d::Gaussian>(g1.get_centroid_ptr(), g1.get_ellipse_ptr(),
                                              std::make_shared<g2d::GaussianIntegralValue>(2));
    CHECK_NE(g1, *g2);
    g2->set_integral_value(g1.get_integral_value());
    CHECK_EQ(g1, *g2);
    auto cen = std::make_shared<g2d::Centroid>(0.5, 0.5);
    g2->set_centroid_ptr(cen);
    CHECK_EQ(g2->get_centroid().get_x(), 0.5);
    CHECK_EQ(g2->get_centroid().get_y(), 0.5);
    CHECK_NE(g1, *g2);
}

TEST_CASE("Gaussians") {
    auto g1 = std::make_shared<g2d::Gaussian>();
    g2d::Gaussians::Data data{g1, g1};
    auto gs = g2d::Gaussians(data);
    CHECK_EQ(&(gs[0]), &(gs[1]));
    CHECK_THROWS_AS(gs.at(2), std::out_of_range);
    for (g2d::Gaussians::const_iterator it = gs.cbegin(); it != gs.cend(); ++it) {
        CHECK_EQ(*it, g1);
    }
    for (const auto& gauss : gs) {
        CHECK_EQ(gauss, g1);
    }

    std::vector<std::optional<const g2d::Gaussians::Data>> data_vec;
    data_vec.emplace_back(g2d::Gaussians::Data());
    data_vec.push_back(data);

    auto gs2 = g2d::Gaussians(data_vec);
    auto gs23 = g2d::Gaussians({data, data});
}

TEST_CASE("ConvolvedGaussian") {
    auto g1 = std::make_shared<g2d::Gaussian>(nullptr, std::make_shared<g2d::Ellipse>(3, 6, 0));
    auto g2 = std::make_shared<g2d::Gaussian>(nullptr, std::make_shared<g2d::Ellipse>(4, 8, 0),
                                              std::make_shared<g2d::GaussianIntegralValue>(2));
    auto gc = g2d::ConvolvedGaussian(g1, g2);
    CHECK_EQ(&(gc.get_source()), &(*g1));
    CHECK_EQ(&(gc.get_kernel()), &(*g2));
}

TEST_CASE("ConvolvedGaussians") {
    auto g1 = std::make_shared<g2d::Gaussian>(nullptr, std::make_shared<g2d::Ellipse>(3, 6, 0));
    auto g2 = std::make_shared<g2d::Gaussian>(nullptr, std::make_shared<g2d::Ellipse>(4, 8, 0),
                                              std::make_shared<g2d::GaussianIntegralValue>(2));
    auto g3 = std::make_shared<g2d::Gaussian>();

    auto gc1 = std::make_shared<g2d::ConvolvedGaussian>(g1, g2);
    auto gc2 = std::make_shared<g2d::ConvolvedGaussian>(g1, g3);
    g2d::ConvolvedGaussians::Data data{gc1, gc2};

    std::vector<std::optional<const g2d::ConvolvedGaussians::Data>> data_vec;
    data_vec.emplace_back(g2d::ConvolvedGaussians::Data());
    data_vec.push_back(data);

    auto gcs = g2d::ConvolvedGaussians(data);
    size_t i = 0;
    for (g2d::ConvolvedGaussians::const_iterator it = gcs.cbegin(); it != gcs.cend(); ++it) {
        CHECK_EQ(*it, data.at(i++));
    }
    i = 0;
    for (const auto& gauss : gcs) {
        CHECK_EQ(gauss, data.at(i++));
    }
}
