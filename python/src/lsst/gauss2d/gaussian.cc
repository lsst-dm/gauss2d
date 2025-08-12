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
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include <memory>
#include <stdexcept>

#include "pybind11.h"

#include "lsst/gauss2d/gaussian.h"

#define PYBIND11_DETAILED_ERROR_MESSAGES

namespace py = pybind11;
using namespace pybind11::literals;

namespace gauss2d = lsst::gauss2d;

void bind_gaussian(py::module &m) {
    m.doc() = "Gauss2D Gaussian Python bindings";

    m.attr("M_HWHM_SIGMA") = py::float_(gauss2d::M_HWHM_SIGMA);
    m.attr("M_SIGMA_HWHM") = py::float_(gauss2d::M_SIGMA_HWHM);

    auto _g = py::classh<gauss2d::GaussianIntegral>(m, "GaussianIntegral");
    py::classh<gauss2d::GaussianIntegralValue, gauss2d::GaussianIntegral>(m, "GaussianIntegralValue")
            .def(py::init<double>(), "value"_a = 1)
            .def_property("value", &gauss2d::GaussianIntegralValue::get_value,
                          &gauss2d::GaussianIntegralValue::set_value)
            .def(py::self == py::self)
            .def(py::self != py::self)
            .def("__repr__",
                 [](const gauss2d::GaussianIntegralValue &self) {
                     return self.repr(true, self.PY_NAMESPACE_SEPARATOR);
                 })
            .def("__str__", &gauss2d::GaussianIntegralValue::str);

    py::classh<gauss2d::Gaussian>(m, "Gaussian")
            .def(py::init<std::shared_ptr<gauss2d::Centroid>, std::shared_ptr<gauss2d::Ellipse>,
                          std::shared_ptr<gauss2d::GaussianIntegral>>(),
                 "centroid"_a = nullptr, "ellipse"_a = nullptr, "integral"_a = nullptr)
            .def_property_readonly("centroid", &gauss2d::Gaussian::get_centroid)
            .def_property_readonly("ellipse", &gauss2d::Gaussian::get_ellipse)
            .def_property_readonly("integral", &gauss2d::Gaussian::get_integral)
            .def_property("const_normal", &gauss2d::Gaussian::get_const_normal,
                          &gauss2d::Gaussian::set_const_normal)
            .def_property("integral_value", &gauss2d::Gaussian::get_integral_value,
                          &gauss2d::Gaussian::set_integral_value)
            .def(py::self == py::self)
            .def(py::self != py::self)
            .def("__repr__",
                 [](const gauss2d::Gaussian &self) { return self.repr(true, self.PY_NAMESPACE_SEPARATOR); })
            .def("__str__", &gauss2d::Gaussian::str);
    py::classh<gauss2d::Gaussians>(m, "Gaussians")
            .def(py::init<std::optional<const gauss2d::Gaussians::Data>>(), "gaussians"_a)
            .def("at", &gauss2d::Gaussians::at, py::return_value_policy::copy)
            .def_property_readonly("size", &gauss2d::Gaussians::size)
            .def("__getitem__", &gauss2d::Gaussians::at_ptr, py::return_value_policy::copy)
            .def("__len__", &gauss2d::Gaussians::size)
            .def("__repr__",
                 [](const gauss2d::Gaussians &self) { return self.repr(true, self.PY_NAMESPACE_SEPARATOR); })
            .def("__str__", &gauss2d::Gaussians::str);
    py::classh<gauss2d::ConvolvedGaussian>(m, "ConvolvedGaussian")
            .def(py::init<std::shared_ptr<gauss2d::Gaussian>, std::shared_ptr<gauss2d::Gaussian>>(),
                 "source"_a = nullptr, "kernel"_a = nullptr)
            .def_property_readonly("kernel", &gauss2d::ConvolvedGaussian::get_kernel)
            .def_property_readonly("source", &gauss2d::ConvolvedGaussian::get_source)
            .def(py::self == py::self)
            .def(py::self != py::self)
            .def("__repr__",
                 [](const gauss2d::ConvolvedGaussian &self) {
                     return self.repr(true, self.PY_NAMESPACE_SEPARATOR);
                 })
            .def("__str__", &gauss2d::ConvolvedGaussian::str);
    py::classh<gauss2d::ConvolvedGaussians>(m, "ConvolvedGaussians")
            .def(py::init<std::optional<const gauss2d::ConvolvedGaussians::Data>>(), "convolvedbgaussians"_a)
            .def("at", &gauss2d::ConvolvedGaussians::at_ptr, py::return_value_policy::copy)
            .def_property_readonly("size", &gauss2d::ConvolvedGaussians::size)
            .def("__getitem__", &gauss2d::ConvolvedGaussians::at_ptr, py::return_value_policy::copy)
            .def("__len__", &gauss2d::ConvolvedGaussians::size)
            .def("__repr__",
                 [](const gauss2d::ConvolvedGaussians &self) {
                     return self.repr(true, self.PY_NAMESPACE_SEPARATOR);
                 })
            .def("__str__", &gauss2d::ConvolvedGaussians::str);
}
