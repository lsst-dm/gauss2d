#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"

#include <memory>

#include "lsst/gauss2d/image.h"
#include "lsst/gauss2d/python/image.h"

#include <pybind11/embed.h>

namespace g2d = lsst::gauss2d;

typedef g2d::python::Image<double> Image;
typedef g2d::ImageArray<double, Image> ImageArray;
typedef g2d::python::Image<bool> Mask;

// Import numpy, otherwise all Image calls will segfault
py::scoped_interpreter guard{};
py::module_ numpy = py::module_::import("numpy");

TEST_CASE("Image") {
    size_t n_rows = 3, n_cols = 2;
    std::shared_ptr<Image> image = std::make_shared<Image>(n_rows, n_cols);
    Image zeros{n_rows, n_cols};
    for (size_t row = 0; row < n_rows; ++row) {
        for (size_t col = 0; col < n_cols; ++col) {
            image->set_value(row, col, 0);
            zeros.set_value(row, col, 0);
        }
    }
    CHECK(*image == zeros);
    CHECK(image->get_coordsys() == g2d::COORDS_DEFAULT);
    CHECK(&(image->get_coordsys()) == &(g2d::COORDS_DEFAULT));
    CHECK(image->get_n_cols() == n_cols);
    CHECK(image->get_n_rows() == n_rows);
    CHECK(image->size() == n_cols * n_rows);

    image->set_value_unchecked(0, 0, 0);
    CHECK(image->get_value_unchecked(0, 0) == 0);
    image->add_value_unchecked(0, 0, 1);
    CHECK(image->get_value_unchecked(0, 0) == 1);
    CHECK(*image != zeros);
    CHECK_THROWS_AS(image->get_value(0, n_cols), std::out_of_range);
    CHECK_THROWS_AS(image->set_value(n_rows, 0, 1), std::out_of_range);
    double& value = image->_get_value_unchecked(1, 1);
    value = -1;
    CHECK(image->get_value_unchecked(1, 1) == -1);
}

TEST_CASE("ImageArray") {
    size_t n_rows = 3, n_cols = 2;
    std::shared_ptr<Image> image = std::make_shared<Image>(n_rows, n_cols);
    ImageArray::Data data = {image, image};
    ImageArray arr = ImageArray(&data);

    for (ImageArray::const_iterator it = arr.cbegin(); it != arr.cend(); ++it) {
        CHECK(*it == image);
    }
    for (const auto& img : arr) {
        CHECK(img == image);
    }
}

TEST_CASE("Mask") {
    auto image = Image(2, 2);
    auto mask = Mask(2, 2);
    mask.fill(false);
    CHECK_EQ(mask.get_value_unchecked(0, 0), false);
    mask._get_value_unchecked(1, 1) = true;
    CHECK_EQ(mask.get_value_unchecked(1, 1), true);
    CHECK_EQ(g2d::images_compatible<double, Image, bool, Mask>(image, mask), true);
}

TEST_CASE("Evaluate") {
    auto gauss1 = std::make_shared<g2d::Gaussian>(nullptr, std::make_shared<g2d::Ellipse>(3, 6, 0));
    auto gauss2 = std::make_shared<g2d::Gaussian>(nullptr, std::make_shared<g2d::Ellipse>(4, 8, 0));
    auto gauss_conv = std::make_shared<g2d::ConvolvedGaussian>(gauss1, gauss2);
    g2d::ConvolvedGaussians::Data data{gauss_conv, std::make_shared<g2d::ConvolvedGaussian>(gauss2, gauss1)};
    auto gaussians_conv = std::make_shared<g2d::ConvolvedGaussians>(data);
    auto output = g2d::make_gaussians_pixel<double, Image, g2d::python::Image<size_t>>(gaussians_conv,
                                                                                       nullptr, 5, 7);
    CHECK_EQ(output->get_n_rows(), 5);
}
