#include <initializer_list>
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"

#include <memory>
#include <stdexcept>

#include "centroid.h"
#include "gaussian.h"

namespace g2 = gauss2d;

TEST_CASE("Gaussian")
{
    auto g1 = g2::Gaussian{};
    auto g2 = std::make_shared<g2::Gaussian>(
        g1.get_centroid_ptr(), g1.get_ellipse_ptr(),
        std::make_shared<g2::GaussianIntegralValue>(2)
    );
    CHECK(g1 != *g2);
    g2->set_integral_value(g1.get_integral_value());
    CHECK(g1 == *g2);
    auto cen = std::make_shared<g2::Centroid>(0.5, 0.5);
    g2->set_centroid_ptr(cen);
    CHECK(g2->get_centroid().get_x() == 0.5);
    CHECK(g2->get_centroid().get_y() == 0.5);
    CHECK(g1 != *g2);
}

TEST_CASE("Gaussians")
{
    auto g1 = std::make_shared<g2::Gaussian>();
    g2::Gaussians::Data data{g1, g1};
    auto gs = g2::Gaussians(data);
    CHECK(&(gs[0]) == &(gs[1]));
    CHECK_THROWS_AS(gs.at(2), std::out_of_range);
    for(g2::Gaussians::const_iterator it = gs.cbegin(); it != gs.cend(); ++it) {
        CHECK(*it == g1);
    }
    for(const auto & gauss : gs) {
        CHECK(gauss == g1);
    }

    std::vector<std::optional<const g2::Gaussians::Data>> data_vec;
    data_vec.emplace_back(g2::Gaussians::Data());
    data_vec.push_back(data);

    auto gs2 = g2::Gaussians(data_vec);
    auto gs23 = g2::Gaussians({data, data});
}

TEST_CASE("ConvolvedGaussian")
{
    auto g1 = std::make_shared<g2::Gaussian>(
        nullptr, std::make_shared<g2::Ellipse>(3, 6, 0));
    auto g2 = std::make_shared<g2::Gaussian>(
        nullptr, std::make_shared<g2::Ellipse>(4, 8, 0),
        std::make_shared<g2::GaussianIntegralValue>(2)
    );
    auto gc = g2::ConvolvedGaussian(g1, g2);
    CHECK(&(gc.get_source_const()) == &(*g1));
        CHECK(&(gc.get_kernel_const()) == &(*g2));
}

TEST_CASE("ConvolvedGaussians")
{
    auto g1 = std::make_shared<g2::Gaussian>(
        nullptr, std::make_shared<g2::Ellipse>(3, 6, 0));
    auto g2 = std::make_shared<g2::Gaussian>(
        nullptr, std::make_shared<g2::Ellipse>(4, 8, 0),
        std::make_shared<g2::GaussianIntegralValue>(2)
    );
    auto g3 = std::make_shared<g2::Gaussian>();

    auto gc1 = std::make_shared<g2::ConvolvedGaussian>(g1, g2);
    auto gc2 = std::make_shared<g2::ConvolvedGaussian>(g1, g3);
    g2::ConvolvedGaussians::Data data{gc1, gc2};

    std::vector<std::optional<const g2::ConvolvedGaussians::Data>> data_vec;
    data_vec.emplace_back(g2::ConvolvedGaussians::Data());
    data_vec.push_back(data);

    auto gcs = g2::ConvolvedGaussians(data);
    size_t i = 0;
    for(g2::ConvolvedGaussians::const_iterator it = gcs.cbegin(); it != gcs.cend(); ++it) {
        CHECK(*it == data.at(i++));
    }
    i = 0;
    for(const auto & gauss : gcs) {
        CHECK(gauss == data.at(i++));
    }
}
