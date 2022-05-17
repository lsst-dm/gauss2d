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

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <memory>
#include <stdexcept>

#include "pybind11.h"

#include "gauss2d/ellipse.h"

namespace py = pybind11;
using namespace pybind11::literals;

void bind_ellipse(py::module &m)
{
    py::class_<gauss2d::Covariance, std::shared_ptr<gauss2d::Covariance>>(m, "Covariance")
        .def(py::init<double, double, double>(), "sigma_x_sq"_a=0, "sigma_y_sq"_a=0, "cov_xy"_a=0)
        .def(py::init<gauss2d::Ellipse&>())
        .def("convolve", &gauss2d::Covariance::convolve)
        .def("make_convolution", &gauss2d::Covariance::make_convolution)
        .def("set", static_cast<void (gauss2d::Covariance::*)(double, double, double)>(
            &gauss2d::Covariance::set),
            "sigma_x_sq"_a=0, "sigma_y_sq"_a=0, "cov_xy"_a=0)
        .def("set", static_cast<void (gauss2d::Covariance::*)(const gauss2d::Ellipse &)>(
            &gauss2d::Covariance::set))
        .def_property("sigma_x_sq", &gauss2d::Covariance::get_sigma_x_sq,
            &gauss2d::Covariance::set_sigma_x_sq)
        .def_property("sigma_y_sq", &gauss2d::Covariance::get_sigma_y_sq,
            &gauss2d::Covariance::set_sigma_y_sq)
        .def_property("cov_xy", &gauss2d::Covariance::get_cov_xy, &gauss2d::Covariance::set_cov_xy)
        .def_property("xyc", &gauss2d::Covariance::get_xyc, &gauss2d::Covariance::set_xyc)
        .def("__repr__", &gauss2d::Covariance::str)
    ;
    auto _e = py::class_<gauss2d::EllipseData, std::shared_ptr<gauss2d::EllipseData>>(m, "EllipseData");
    py::class_<gauss2d::EllipseValues,
        std::shared_ptr<gauss2d::EllipseValues>,
        gauss2d::EllipseData
    >(m, "EllipseValues")
        .def(py::init<double, double, double>(), "sigma_x"_a=0, "sigma_y"_a=0, "rho"_a=0)
        .def("set", &gauss2d::EllipseValues::set)
        .def_property("rho", &gauss2d::EllipseValues::get_rho, &gauss2d::EllipseValues::set_rho)
        .def_property("sigma_x", &gauss2d::EllipseValues::get_sigma_x, &gauss2d::EllipseValues::set_sigma_x)
        .def_property("sigma_y", &gauss2d::EllipseValues::get_sigma_y, &gauss2d::EllipseValues::set_sigma_y)
        .def_property("xyr", &gauss2d::EllipseValues::get_xyr, &gauss2d::EllipseValues::set_xyr)
        .def("__repr__", &gauss2d::EllipseValues::str)
    ;
    py::class_<gauss2d::Ellipse, std::shared_ptr<gauss2d::Ellipse>>(m, "Ellipse")
        .def(py::init<std::shared_ptr<gauss2d::EllipseData>>(), "data"_a)
        .def(py::init<gauss2d::Covariance&>())
        .def(py::init<gauss2d::EllipseMajor&>())
        .def(py::init<double, double, double>(), "sigma_x"_a=0, "sigma_y"_a=0, "rho"_a=0)
        .def_static("check", &gauss2d::Ellipse::check)
        .def("convolve", &gauss2d::Ellipse::convolve)
        .def("get_cov_xy", &gauss2d::Ellipse::get_cov_xy)
        .def("get_radius_trace", &gauss2d::Ellipse::get_radius_trace)
        .def("make_convolution", &gauss2d::Ellipse::make_convolution)
        .def("set", static_cast<void (gauss2d::Ellipse::*)(double, double, double)>(&gauss2d::Ellipse::set))
        .def("set", static_cast<void (gauss2d::Ellipse::*)(const gauss2d::Covariance&)>(&gauss2d::Ellipse::set))
        .def("set", static_cast<void (gauss2d::Ellipse::*)(const gauss2d::EllipseMajor&)>(&gauss2d::Ellipse::set))
        .def_property("rho", &gauss2d::Ellipse::get_rho, &gauss2d::Ellipse::set_rho)
        .def_property("sigma_x", &gauss2d::Ellipse::get_sigma_x, &gauss2d::Ellipse::set_sigma_x)
        .def_property("sigma_y", &gauss2d::Ellipse::get_sigma_y, &gauss2d::Ellipse::set_sigma_y)
        .def_property("xyr", &gauss2d::Ellipse::get_xyr, &gauss2d::Ellipse::set_xyr)
        .def("__repr__", &gauss2d::Ellipse::str)
    ;
    py::class_<gauss2d::EllipseMajor, std::shared_ptr<gauss2d::EllipseMajor> >(m, "EllipseMajor")
        .def(py::init<gauss2d::Covariance&, bool>(), "covariance"_a, "degrees"_a=false)
        .def(py::init<gauss2d::Ellipse&, bool>(), "ellipse"_a, "degrees"_a=false)
        .def(py::init<double, double, double, bool>(), "r_major"_a=0, "axrat"_a=1, "angle"_a=0,
            "degrees"_a=false)
        .def_static("check", &gauss2d::EllipseMajor::check)
        .def("get_angle_degrees", &gauss2d::EllipseMajor::get_angle_degrees)
        .def("get_angle_radians", &gauss2d::EllipseMajor::get_angle_radians)
        .def_property("r_major", &gauss2d::EllipseMajor::get_r_major, &gauss2d::EllipseMajor::set_r_major)
        .def_property("axrat", &gauss2d::EllipseMajor::get_axrat, &gauss2d::EllipseMajor::set_axrat)
        .def_property("angle", &gauss2d::EllipseMajor::get_angle, &gauss2d::EllipseMajor::set_angle)
        .def_property("degrees", &gauss2d::EllipseMajor::is_degrees, &gauss2d::EllipseMajor::set_degrees)
        .def_property("rqa", &gauss2d::EllipseMajor::get_rqa, &gauss2d::EllipseMajor::set_rqa)
        .def("__repr__", &gauss2d::EllipseMajor::str)
    ;
}

