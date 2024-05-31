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

#include "lsst/gauss2d/centroid.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace gauss2d = lsst::gauss2d;

void bind_centroid(py::module &m) {
    auto _c = py::class_<gauss2d::CentroidData, std::shared_ptr<gauss2d::CentroidData>>(m, "CentroidData");
    py::class_<gauss2d::CentroidValues, std::shared_ptr<gauss2d::CentroidValues>, gauss2d::CentroidData>(
            m, "CentroidValues")
            // Note: These can't be shared_ptrs
            .def(py::init<double, double>(), "x"_a = 0, "y"_a = 0)
            .def_property("x", &gauss2d::CentroidValues::get_x, &gauss2d::CentroidValues::set_x)
            .def_property("y", &gauss2d::CentroidValues::get_y, &gauss2d::CentroidValues::set_y)
            .def_property("xy", &gauss2d::CentroidValues::get_xy, &gauss2d::CentroidValues::set_xy)
            .def("__repr__", [](const gauss2d::CentroidValues &self) { return self.repr(true); })
            .def("__str__", &gauss2d::CentroidValues::str);
    py::class_<gauss2d::Centroid, std::shared_ptr<gauss2d::Centroid>>(m, "Centroid")
            .def(py::init<std::shared_ptr<gauss2d::CentroidData>>(), "data"_a)
            .def(py::init<double, double>(), "x"_a = 0, "y"_a = 0)
            .def_property("x", &gauss2d::Centroid::get_x, &gauss2d::Centroid::set_x)
            .def_property("y", &gauss2d::Centroid::get_y, &gauss2d::Centroid::set_y)
            .def_property("xy", &gauss2d::Centroid::get_xy, &gauss2d::Centroid::set_xy)
            .def("__repr__", [](const gauss2d::Centroid &self) { return self.repr(true); })
            .def("__str__", &gauss2d::Centroid::str);
}
