// -*- LSST-C++ -*-
/*
 * This file is part of gauss2d.
 *
 * Developed for the LSST Data Management System.
 * This product includes software developed by the LSST Project
 * (https://www.lsst.org).
 * See the COPYRIGHT file at the top-level directory of this distribution
 * for details of code ownership.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"

#include <experimental/iterator>
#include <memory>
#include <sstream>

#include "lsst/gauss2d/centroid.h"
#include "lsst/gauss2d/ellipse.h"
#include "lsst/gauss2d/evaluate.h"
#include "lsst/gauss2d/vectorimage.h"

namespace g2d = lsst::gauss2d;

typedef g2d::VectorImage<double> Image;
typedef g2d::ImageArray<double, Image> ImageArray;
typedef g2d::VectorImage<size_t> Indices;
typedef g2d::GaussianEvaluator<double, Image, Indices> Evaluator;

typedef std::shared_ptr<double> Value;

TEST_CASE("Evaluator") {
    const size_t n_rows = 5, n_cols = 4, n_comp = 2;
    auto image = std::make_shared<Image>(n_rows, n_cols);

    g2d::ConvolvedGaussians::Data data{};

    std::vector<Value> values;
    values.reserve(n_comp * g2d::N_PARAMS_GAUSS2D);

    for (size_t i = 0; i < n_comp; ++i) {
        auto x = std::make_shared<double>(n_cols / 2);
        auto y = std::make_shared<double>(n_rows / 2);
        values.push_back(x);
        values.push_back(y);

        auto integral = std::make_shared<double>(i);
        values.push_back(integral);

        auto sigma_x = std::make_shared<double>(n_cols / (i + 2));
        auto sigma_y = std::make_shared<double>(n_rows / (i + 2));
        auto rho = std::make_shared<double>(0);
        values.push_back(sigma_x);
        values.push_back(sigma_y);
        values.push_back(rho);

        data.emplace_back(std::make_shared<g2d::ConvolvedGaussian>(
                std::make_shared<const g2d::Gaussian>(
                        std::make_shared<g2d::Centroid>(
                                std::make_shared<g2d::CentroidValues>(std::move(x), std::move(y))),
                        std::make_shared<g2d::Ellipse>(std::make_shared<g2d::EllipseValues>(
                                std::move(sigma_x), std::move(sigma_y), std::move(rho))),
                        std::make_shared<g2d::GaussianIntegralValue>(std::move(integral))),
                std::make_shared<const g2d::Gaussian>(std::make_shared<g2d::Centroid>(0, 0),
                                                      std::make_shared<g2d::Ellipse>(0, 0, 0),
                                                      std::make_shared<g2d::GaussianIntegralValue>())));
    }
    const size_t n_params = values.size();
    CHECK_EQ(n_params, n_comp * g2d::N_PARAMS_GAUSS2D);

    auto gaussians = std::make_shared<const g2d::ConvolvedGaussians>(data);
    auto sigma_inv = std::make_shared<Image>(n_rows, n_cols);
    const double err = 1e-5;
    sigma_inv->fill(1 / err);

    auto eval_img = std::make_shared<Evaluator>(gaussians, nullptr, nullptr, image);

    CHECK_NE(image, nullptr);
    CHECK_EQ(eval_img->get_n_cols(), n_cols);
    CHECK_EQ(eval_img->get_n_rows(), n_rows);
    CHECK_EQ(eval_img->get_size(), n_cols * n_rows);

    double loglike_img = eval_img->loglike_pixel();
    CHECK_EQ(loglike_img, 0);

    double x_min = 2.0;
    double y_min = 1.0;
    auto coordsys2 = std::make_shared<g2d::CoordinateSystem>(1., 1., x_min, y_min);
    auto image2 = std::make_shared<Image>(n_rows, n_cols, Image::_value_default_ptr(), coordsys2);
    CHECK_EQ(*coordsys2, image2->get_coordsys());
    auto eval_offset = std::make_shared<Evaluator>(gaussians, nullptr, nullptr, image2);
    CHECK_EQ(eval_offset->get_coordsys(), *coordsys2);

    eval_offset->loglike_pixel();
    CHECK_EQ(image2->get_value(0, 0), image->get_value(y_min, x_min));

    // Add a small offset so the chi is not zero everywhere
    *image += err;

    auto eval_like = std::make_shared<Evaluator>(gaussians, image, sigma_inv);

    CHECK_NE(sigma_inv, nullptr);

    double loglike_like = eval_like->loglike_pixel();
    CHECK_NE(loglike_like, 0);

    auto img_loglike_grads = std::make_shared<Image>(1, n_params);
    ImageArray::Data data_loglike_grads = {img_loglike_grads};

    Evaluator eval_loglike_grad(gaussians, image, sigma_inv, nullptr,
                                nullptr,  // residual,
                                std::make_shared<ImageArray>(&data_loglike_grads));
    eval_loglike_grad.loglike_pixel();

    double dx = 1e-8;

    for (size_t idx_param = 0; idx_param < n_params; idx_param++) {
        double& value = *(values[idx_param]);
        double value_old = value;
        value += dx;
        double dloglike_findif = eval_like->loglike_pixel();
        value = value_old - dx;
        dloglike_findif = (dloglike_findif - eval_like->loglike_pixel()) / (2 * dx);
        double dloglike_eval = img_loglike_grads->get_value(0, idx_param);
        double delta_dloglike_max = 1e-3 * std::abs((dloglike_eval + dloglike_findif) / 2.) + 1e-4;
        CHECK_LT(std::abs(dloglike_eval - dloglike_findif), delta_dloglike_max);
        value = value_old;
    }

    ImageArray::Data data_jacs;
    for (size_t idx_param = 0; idx_param < n_params; idx_param++) {
        data_jacs.push_back(std::make_shared<Image>(n_rows, n_cols));
    }

    auto jacs = std::make_shared<ImageArray>(&data_jacs);

    Evaluator eval_jacob(gaussians, image, sigma_inv, nullptr,
                         nullptr,  // residual,
                         jacs
                         /*      map_grad_in,
                                 factors_grad_in,
                                 map_extra_in,
                                 factors_extra_in
                                 //, background
                         */
    );

    CHECK_NE(gaussians, nullptr);

    double loglike_jacob = eval_jacob.loglike_pixel();
    CHECK_EQ(loglike_jacob, loglike_like);

    // Subtract that small offset back from the image
    *image += -err;

    auto image_param = std::make_shared<Image>(n_rows, n_cols);
    auto eval_img_new = std::make_shared<Evaluator>(gaussians, nullptr, nullptr, image_param);

    const double eps = 1e-6;
    const double atol = 1e-4;
    const double rtol = 1e-4;

    std::vector<std::string> errors;

    for (size_t idx_param = 0; idx_param < 10; idx_param++) {
        const auto& jacs_param = jacs->at(idx_param);
        const auto& param = values[idx_param];
        *param += eps;

        eval_img_new->loglike_pixel();

        size_t n_failed = 0;
        for (unsigned int i = 0; i < n_cols; ++i) {
            for (unsigned int j = 0; j < n_rows; ++j) {
                double delta = sigma_inv->get_value(j, i) / eps
                               * (image_param->get_value(j, i) - image->get_value(j, i));
                double jac = jacs_param.get_value(j, i);
                if (jac != delta) {
                    double diff_grad = delta - jac;
                    double diff_grad_abs = std::abs(diff_grad);
                    // Same as numpy.is_close
                    const double margin = atol + rtol * std::abs(jac);
                    n_failed += !(diff_grad_abs <= margin);
                }
            }
        }

        if ((n_failed > 0)) {
            errors.push_back("Param[" + std::to_string(idx_param) + "] n_failed=" + std::to_string(n_failed));
        }

        *param -= eps;
    }

    std::string errormsg = "";
    if (!errors.empty()) {
        std::stringstream ss;
        ss << "";
        std::copy(std::begin(errors), std::end(errors), std::experimental::make_ostream_joiner(ss, "\n"));
        errormsg = ss.str();
    }
    CHECK_EQ(errormsg, "");
}
