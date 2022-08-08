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

#include "gauss2d/gaussian.h"

namespace py = pybind11;
using namespace pybind11::literals;

void bind_gaussian(py::module &m)
{
    m.doc() = "Gauss2D Gaussian Python bindings";

    m.attr("M_HWHM_SIGMA") = py::float_(gauss2d::M_HWHM_SIGMA);
    m.attr("M_SIGMA_HWHM") = py::float_(gauss2d::M_SIGMA_HWHM);

    auto _g = py::class_<gauss2d::GaussianIntegral, std::shared_ptr<gauss2d::GaussianIntegral>>(
        m, "GaussianIntegral");
    py::class_<gauss2d::GaussianIntegralValue, gauss2d::GaussianIntegral,
            std::shared_ptr<gauss2d::GaussianIntegralValue>>(m, "GaussianIntegralValue")
        .def(py::init<double>(), "value"_a=1)
        .def_property("value", &gauss2d::GaussianIntegralValue::get_value, &gauss2d::GaussianIntegralValue::set_value)
        .def("__repr__", &gauss2d::GaussianIntegralValue::str)
    ;

    py::class_<gauss2d::Gaussian, std::shared_ptr<gauss2d::Gaussian>>(m, "Gaussian")
        .def(py::init<std::shared_ptr<gauss2d::Centroid>, std::shared_ptr<gauss2d::Ellipse>,
            std::shared_ptr<gauss2d::GaussianIntegral>>(),
            "centroid"_a=nullptr, "ellipse"_a=nullptr, "integral"_a=nullptr)
        .def_property_readonly("centroid", &gauss2d::Gaussian::get_centroid)
        .def_property_readonly("ellipse", &gauss2d::Gaussian::get_ellipse)
        .def_property_readonly("integral", &gauss2d::Gaussian::get_integral)
        .def_property("const_normal", &gauss2d::Gaussian::get_const_normal,
            &gauss2d::Gaussian::set_const_normal)
        .def_property("integral_value", &gauss2d::Gaussian::get_integral_value,
            &gauss2d::Gaussian::set_integral_value)
        .def("__repr__", &gauss2d::Gaussian::str)
    ;
    py::class_<gauss2d::Gaussians, std::shared_ptr<gauss2d::Gaussians>>(m, "Gaussians")
        .def(py::init<std::optional<const gauss2d::Gaussians::Data>>(), "gaussians"_a)
        .def("at", &gauss2d::Gaussians::at)
        .def_property_readonly("size", &gauss2d::Gaussians::size)
        .def("__len__", &gauss2d::Gaussians::size)
        .def("__repr__", &gauss2d::Gaussians::str)
    ;
    py::class_<gauss2d::ConvolvedGaussian, std::shared_ptr<gauss2d::ConvolvedGaussian>>(
        m, "ConvolvedGaussian")
        .def(py::init<std::shared_ptr<gauss2d::Gaussian>, std::shared_ptr<gauss2d::Gaussian>>(),
            "source"_a = nullptr, "kernel"_a = nullptr)
        .def_property_readonly("kernel", &gauss2d::ConvolvedGaussian::get_kernel)
        .def_property_readonly("source", &gauss2d::ConvolvedGaussian::get_source)
        .def("__repr__", &gauss2d::ConvolvedGaussian::str)
    ;
    py::class_<gauss2d::ConvolvedGaussians, std::shared_ptr<gauss2d::ConvolvedGaussians>>(
        m, "ConvolvedGaussians")
        .def(py::init<std::optional<const gauss2d::ConvolvedGaussians::Data>>(), "convolvedbgaussians"_a)
        .def("at", &gauss2d::ConvolvedGaussians::at)
        .def_property_readonly("size", &gauss2d::ConvolvedGaussians::size)
        .def("__len__", &gauss2d::ConvolvedGaussians::size)
        .def("__repr__", &gauss2d::ConvolvedGaussians::str)
    ;
}

