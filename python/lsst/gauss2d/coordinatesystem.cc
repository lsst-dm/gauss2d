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

#include "pybind11.h"

#include "lsst/gauss2d/coordinatesystem.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace gauss2d = lsst::gauss2d;

void bind_coordinatesystem(py::module &m) {
    py::class_<gauss2d::CoordinateSystem, std::shared_ptr<gauss2d::CoordinateSystem>>(m, "CoordinateSystem")
            .def(py::init<double, double, double, double>(), "dx1"_a = 1., "dy2"_a = 1., "x_min"_a = 0.,
                 "y_min"_a = 0.)
            .def_property_readonly("dx1", &gauss2d::CoordinateSystem::get_dx1)
            .def_property_readonly("dy2", &gauss2d::CoordinateSystem::get_dy2)
            .def_property_readonly("x_min", &gauss2d::CoordinateSystem::get_x_min)
            .def_property_readonly("y_min", &gauss2d::CoordinateSystem::get_y_min)
            .def(py::self == py::self)
            .def(py::self != py::self)
            .def("__repr__", [](const gauss2d::CoordinateSystem &self) { return self.repr(true); })
            .def("__str__", &gauss2d::CoordinateSystem::str);
}
