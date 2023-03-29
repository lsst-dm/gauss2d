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
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <memory>
#include <stdexcept>

#include "pybind11.h"

#include "gauss2d/evaluate.h"
#include "gauss2d/image.h"
#include "gauss2d/python/pyimage.h"
//#include "gauss2d/integrator.h"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace gauss2d::python;

void bind_image(py::module &m)
{
    py::class_<gauss2d::CoordinateSystem, std::shared_ptr<gauss2d::CoordinateSystem>>(m, "CoordinateSystem")
        .def(py::init<double, double>(), "d1"_a=1., "d2"_a=1.)
        .def("__repr__", [](const gauss2d::CoordinateSystem & self) { return self.repr(true); })
        .def("__str__", &gauss2d::CoordinateSystem::str)
    ;
    declare_image<bool>(m, "B");
    declare_image<float>(m, "F");
    declare_image<double>(m, "D");
    declare_image<int>(m, "I");
    declare_image<unsigned int>(m, "U");
    declare_image<size_t>(m, "S");
    declare_image_array<bool>(m, "B");
    declare_image_array<float>(m, "F");
    declare_image_array<double>(m, "D");
    declare_image_array<int>(m, "I");
    declare_image_array<unsigned int>(m, "U");
    declare_image_array<size_t>(m, "S");
    declare_evaluator<float>(m, "F");
    declare_evaluator<double>(m, "D");
    declare_maker<float, PyImage<float>, PyImage<size_t>>(m, "F");
    declare_maker<double, PyImage<double>, PyImage<size_t>>(m, "D");
}
