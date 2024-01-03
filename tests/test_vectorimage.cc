#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <memory>

#include "image.h"
#include "vectorimage.h"

namespace g2 = gauss2d;

typedef g2::VectorImage<double> Image;
typedef g2::ImageArray<double, Image> ImageArray;
typedef g2::VectorImage<bool> Mask;

TEST_CASE("VectorImage") {
    const double value_init = -42.1;
    Image nonzero{1, 1, &value_init};
    CHECK(nonzero.get_value(0, 0) == value_init);

    size_t n_rows = 3, n_cols = 2;
    std::shared_ptr<Image> image = std::make_shared<Image>(n_rows, n_cols);
    Image zeros{n_rows, n_cols};
    auto x = Image(n_rows, n_cols) == zeros;

    CHECK(x);

    CHECK(image->get_coordsys() == g2::COORDS_DEFAULT);
    CHECK(&(image->get_coordsys()) == &(g2::COORDS_DEFAULT));

    CHECK(image->get_n_cols() == n_cols);
    CHECK(image->get_n_rows() == n_rows);
    CHECK(image->size() == n_cols * n_rows);

    image->add_value_unchecked(0, 0, 1);
    CHECK(*image != zeros);
    CHECK(image->get_value_unchecked(0, 0) == 1);
    CHECK_THROWS_AS(image->get_value(n_rows, n_cols), std::out_of_range);
    CHECK_THROWS_AS(image->set_value(n_rows, n_cols, 1), std::out_of_range);
    double& value = image->_get_value_unchecked(1, 1);
    value = -1;
    CHECK(image->get_value_unchecked(1, 1) == -1);
}

TEST_CASE("VectorImageArray") {
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

TEST_CASE("VectorMask") {
    size_t n_cr = 2;
    auto image = Image(n_cr, n_cr);
    auto mask = Mask(n_cr, n_cr);
    CHECK_THROWS_AS(mask.get_value(n_cr, n_cr-1), std::out_of_range);
    CHECK_THROWS_AS(mask.get_value(n_cr-1, n_cr), std::out_of_range);
    mask.set_value_unchecked(0, 0, false);
    CHECK(mask.get_value_unchecked(0, 0) == false);
    // Can't do this because vector<bool> can't return references
    // mask._get_value_unchecked(1, 1) = true;
    mask.set_value(1, 1, true);
    CHECK(mask.get_value_unchecked(1, 1) == true);
    mask.set_value(1, 1, false);
    CHECK_THROWS_AS(mask.set_value(n_cr, n_cr-1, false), std::out_of_range);
    CHECK_THROWS_AS(mask.set_value(n_cr-1, n_cr, true), std::out_of_range);
    CHECK(mask.get_value(1, 1) == false);
    CHECK(g2::images_compatible<double, Image, bool, Mask>(image, mask));
}
