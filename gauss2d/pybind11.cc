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

#ifndef __GAUSS2D_EVALUATE_H_
#include "../include/evaluate.h"
#endif

#ifndef __GAUSS2D_GAUSSIAN_H_
#include "../include/gaussian.h"
#endif

/*#ifndef GAUSS2D_GAUSSIAN_INTEGRATOR_H
#include "../src/gaussian_integrator.h"
#endif*/

#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(_gauss2d, m)
{
    m.doc() = "Gauss2D Pybind11 functions"; // optional module docstring
    m.attr("M_HWHM_SIGMA") = py::float_(gauss2d::M_HWHM_SIGMA);
    m.attr("M_SIGMA_HWHM") = py::float_(gauss2d::M_SIGMA_HWHM);

    py::class_<gauss2d::Centroid, std::shared_ptr<gauss2d::Centroid>>(m, "Centroid")
        .def(py::init<double, double>());
    py::class_<gauss2d::Covariance, std::shared_ptr<gauss2d::Covariance> >(m, "Covariance")
        .def(py::init<double, double, double>(), "sigma_x"_a=0, "sigma_y"_a=0, "rho"_a=0)
        .def(py::init<gauss2d::Ellipse&>())
        .def("make_ellipse_major", &gauss2d::Covariance::make_ellipse_major, "degrees"_a=false)
        .def_property("sigma_x_sq", &gauss2d::Covariance::get_sigma_x_sq, &gauss2d::Covariance::set_sigma_x_sq)
        .def_property("sigma_y_sq", &gauss2d::Covariance::get_sigma_y_sq, &gauss2d::Covariance::set_sigma_y_sq)
        .def_property("cov_xy", &gauss2d::Covariance::get_cov_xy, &gauss2d::Covariance::set_cov_xy)
        ;
    py::class_<gauss2d::Ellipse, std::shared_ptr<gauss2d::Ellipse> >(m, "Ellipse")
        .def(py::init<gauss2d::EllipseTerms&>())
        .def(py::init<double, double, double>(), "sigma_x"_a=0, "sigma_y"_a=0, "rho"_a=0)
        .def("_set", &gauss2d::Ellipse::_set)
        .def("_set_rho", &gauss2d::Ellipse::_set_rho)
        .def("_set_sigma_x", &gauss2d::Ellipse::_set_sigma_x)
        .def("_set_sigma_y", &gauss2d::Ellipse::_set_sigma_y)
        .def_static("check", &gauss2d::Ellipse::check)
        .def("convolve", &gauss2d::Ellipse::convolve)
        .def("get", &gauss2d::Ellipse::get)
        .def("get_cov_xy", &gauss2d::Ellipse::get_cov_xy)
        .def("get_radius_trace", &gauss2d::Ellipse::get_radius_trace)
        .def("make_convolution", &gauss2d::Ellipse::make_convolution)
        .def("make_ellipse_major", &gauss2d::Ellipse::make_ellipse_major, "degrees"_a=false)
        .def_property("rho", &gauss2d::Ellipse::get_rho, &gauss2d::Ellipse::set_rho)
        .def_property("sigma_x", &gauss2d::Ellipse::get_sigma_x, &gauss2d::Ellipse::set_sigma_x)
        .def_property("sigma_y", &gauss2d::Ellipse::get_sigma_y, &gauss2d::Ellipse::set_sigma_y)
    ;
    py::class_<gauss2d::EllipseMajor, std::shared_ptr<gauss2d::EllipseMajor> >(m, "EllipseMajor")
        .def(py::init<double, double, double, bool>(), "r_major"_a=0, "axrat"_a=1, "angle"_a=0, "degrees"_a=false)
        .def_static("check", &gauss2d::EllipseMajor::check)
        .def_property("r_major", &gauss2d::EllipseMajor::get_r_major, &gauss2d::EllipseMajor::set_r_major)
        .def_property("axrat", &gauss2d::EllipseMajor::get_axrat, &gauss2d::EllipseMajor::set_axrat)
        .def_property("angle", &gauss2d::EllipseMajor::get_angle, &gauss2d::EllipseMajor::set_angle)
        .def_property("degrees", &gauss2d::EllipseMajor::is_degrees, &gauss2d::EllipseMajor::set_degrees)
    ;
    py::class_<gauss2d::EllipseTerms, std::shared_ptr<gauss2d::EllipseTerms>>(m, "EllipseTerms")
        .def(py::init<double, double, double>(), "sigma_x"_a=0, "sigma_y"_a=0, "rho"_a=0)
        ;
    py::class_<gauss2d::Gaussian, std::shared_ptr<gauss2d::Gaussian>>(m, "Gaussian")
        .def(py::init<std::shared_ptr<gauss2d::Centroid>, std::shared_ptr<gauss2d::Ellipse>>())
        ;

/*    m.def(
        "make_gaussian", &gauss2d::make_gaussian,
        "Integrate a 2D Gaussian over a rectangular grid.",
        "cen_x"_a, "cen_x"_a, "mag"_a, "r_eff"_a, "axrat"_a, "ang"_a,
        "x_min"_a, "x_max"_a, "y_min"_a, "y_max"_a, "dim_x"_a, "dim_y"_a, "acc"_a
    );*/

    m.def(
        "make_gaussian_pixel", &gauss2d::make_gaussian_pixel,
        "Evaluate a 2D Gaussian at the centers of pixels on a rectangular grid using the standard bivariate"
        "Gaussian PDF.",
        "cen_x"_a, "cen_x"_a, "L"_a, "r_eff"_a, "axrat"_a, "ang"_a,
        "x_min"_a, "x_max"_a, "y_min"_a, "y_max"_a, "dim_x"_a, "dim_y"_a
    );

    m.def(
        "make_gaussian_pixel_sersic", &gauss2d::make_gaussian_pixel_sersic,
        "Evaluate a 2D Gaussian at the centers of pixels on a rectangular grid using the 2D Sersic PDF.",
        "cen_x"_a, "cen_y"_a, "L"_a, "r_eff"_a,  "axrat"_a, "ang"_a,
        "x_min"_a, "x_max"_a, "y_min"_a, "y_max"_a, "dim_x"_a, "dim_y"_a
    );

    m.def(
        "make_gaussian_pixel_covar", &gauss2d::make_gaussian_pixel_covar,
        "Evaluate a 2D Gaussian at the centers of pixels on a rectangular grid using the standard bivariate"
        "Gaussian PDF and given a covariance matrix.",
        "cen_x"_a, "cen_y"_a, "L"_a, "sig_x"_a, "sig_y"_a, "rho"_a,
        "x_min"_a, "x_max"_a, "y_min"_a, "y_max"_a, "dim_x"_a, "dim_y"_a
    );

    m.def(
        "make_gaussians_pixel", &gauss2d::make_gaussians_pixel,
        "Evaluate 2D Gaussians at the centers of pixels on a rectangular grid using the standard bivariate"
        "Gaussian PDF.",
        "gaussians"_a.noconvert(), "x_min"_a, "x_max"_a, "y_min"_a, "y_max"_a, "dim_x"_a, "dim_y"_a
    );

    m.def(
        "add_gaussians_pixel", &gauss2d::add_gaussians_pixel,
        "Evaluate 2D Gaussians at the centers of pixels on a rectangular grid using the standard bivariate"
        "Gaussian PDF, adding to an existing matrix.",
        "gaussians"_a.noconvert(), "x_min"_a, "x_max"_a, "y_min"_a, "y_max"_a, "output"_a.noconvert()
    );

    m.def(
        "loglike_gaussians_pixel", &gauss2d::loglike_gaussians_pixel,
        "Evaluate the log likelihood of a 2D Gaussian mixture model at the centers of pixels on a rectangular"
        "grid using the standard bivariate Gaussian PDF.",
        "data"_a.noconvert(), "sigma_inv"_a.noconvert(), "gaussians"_a.noconvert(),
        "x_min"_a, "x_max"_a, "y_min"_a, "y_max"_a, "to_add"_a, "output"_a.noconvert(),
        "residual"_a.noconvert(), "grad"_a.noconvert(),
        "grad_param_map"_a.noconvert(), "grad_param_factor"_a.noconvert(),
        "sersic_param_map"_a.noconvert(), "sersic_param_factor"_a.noconvert(),
        "background"_a.noconvert()
    );
}
