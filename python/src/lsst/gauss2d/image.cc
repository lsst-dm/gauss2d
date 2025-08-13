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

#include "lsst/gauss2d/evaluate.h"
#include "lsst/gauss2d/image.h"
#include "lsst/gauss2d/python/image.h"

#include "pybind11.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace g2d = lsst::gauss2d;

void bind_image(py::module &m) {
    g2d::python::declare_image<bool>(m, "B");
    g2d::python::declare_image<float>(m, "F");
    g2d::python::declare_image<double>(m, "D");
    g2d::python::declare_image<int>(m, "I");
    g2d::python::declare_image<unsigned int>(m, "U");
    g2d::python::declare_image<size_t>(m, "S");
    g2d::python::declare_image_array<bool>(m, "B");
    g2d::python::declare_image_array<float>(m, "F");
    g2d::python::declare_image_array<double>(m, "D");
    g2d::python::declare_image_array<int>(m, "I");
    g2d::python::declare_image_array<unsigned int>(m, "U");
    g2d::python::declare_image_array<size_t>(m, "S");
    g2d::python::declare_evaluator<float>(m, "F");
    g2d::python::declare_evaluator<double>(m, "D");
    g2d::python::declare_maker<float, g2d::python::Image<float>, g2d::python::Image<size_t>>(m, "F");
    g2d::python::declare_maker<double, g2d::python::Image<double>, g2d::python::Image<size_t>>(m, "D");
}
