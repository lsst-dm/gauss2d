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

#ifndef LSST_GAUSS2D_PYTHON_PYIMAGE_H
#define LSST_GAUSS2D_PYTHON_PYIMAGE_H

#include <memory>
#include <stdexcept>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "lsst/gauss2d/evaluate.h"
#include "lsst/gauss2d/image.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace gauss2d = lsst::gauss2d;

namespace lsst::gauss2d::python {

/*
    Convenience functions for binding concrete image types.
*/

#pragma GCC visibility push(hidden)
template <typename t>
class PyImage : public gauss2d::Image<t, PyImage<t>> {
private:
    const size_t _n_rows;
    const size_t _n_cols;

    std::unique_ptr<py::array_t<t>> _ptr_own = nullptr;
    py::array_t<t> _data;
    py::detail::unchecked_mutable_reference<t, 2> _data_ref;

    void _validate() const {
        if ((_ptr_own == nullptr ? _data : *_ptr_own).ndim() != 2) {
            throw std::runtime_error("Input data must have 2 dimensions");
        }
    }

public:
    inline t &_get_value_unchecked(size_t row, size_t col) { return this->_data_ref(row, col); };

    py::array_t<t> &get_data() {
        py::array_t<t> &rv = _ptr_own == nullptr ? _data : *_ptr_own;
        return rv;
    }

    size_t get_n_cols() const { return _n_cols; };
    size_t get_n_rows() const { return _n_rows; };

    // void add_value(size_t row, size_t col, t value) { this->_get_value(row, col) += value;}
    // void add_value_unchecked(size_t row, size_t col, t value) {
    //     _get_value_unchecked(row, col) += value;
    // }
    const inline t get_value_unchecked(size_t row, size_t col) const { return this->_data_ref(row, col); };
    void set_value_unchecked(size_t row, size_t col, t value) { this->_data_ref(row, col) = value; }
    // void set_value_unchecked(size_t row, size_t col, t value) { _get_value_unchecked(row, col) = value;};

    PyImage(size_t n_rows, size_t n_cols,
            const std::shared_ptr<const gauss2d::CoordinateSystem> coordsys = nullptr)
            : gauss2d::Image<t, PyImage<t>>(coordsys),
              _n_rows(n_rows),
              _n_cols(n_cols),
              _ptr_own(std::make_unique<py::array_t<t>>(py::array::ShapeContainer({n_rows, n_cols}))),
              _data_ref(this->get_data().template mutable_unchecked<2>()) {
        _validate();
    }
    PyImage(py::array_t<t> data, const std::shared_ptr<const gauss2d::CoordinateSystem> coordsys = nullptr)
            : gauss2d::Image<t, PyImage<t>>(coordsys),
              _n_rows(data.shape(0)),
              _n_cols(data.shape(1)),
              _data(data),
              _data_ref(this->get_data().template mutable_unchecked<2>()) {
        _validate();
    }

    ~PyImage(){};
};
#pragma GCC visibility pop

template <typename T>
void declare_image(py::module &m, std::string typestr) {
    using Class = PyImage<T>;
    std::string pyclass_name = std::string("Image") + typestr;
    py::class_<Class, std::shared_ptr<Class>>(m, pyclass_name.c_str())
            .def(py::init<size_t, size_t, const std::shared_ptr<const gauss2d::CoordinateSystem>>(),
                 "n_rows"_a, "n_cols"_a, "coordsys"_a = gauss2d::COORDS_DEFAULT)
            .def(py::init<py::array_t<T>, const std::shared_ptr<const gauss2d::CoordinateSystem>>(), "data"_a,
                 "coordsys"_a = gauss2d::COORDS_DEFAULT)
            .def_property_readonly("coordsys", &Class::get_coordsys_ptr_const)
            .def_property_readonly("data", &Class::get_data)
            .def_property_readonly("n_rows", &Class::get_n_rows)
            .def_property_readonly("n_cols", &Class::get_n_cols)
            .def("fill", &Class::fill)
            .def("get_value", &Class::get_value)
            .def("set_value", &Class::set_value)
            .def("get_value_unchecked", &Class::get_value_unchecked)
            .def("set_value_unchecked", &Class::set_value_unchecked)
            .def_property_readonly("size", &Class::size)
            .def_property_readonly("shape", &Class::shape)
            .def(py::self == py::self)
            .def(py::self != py::self)
            .def(py::self += T())
            .def("__repr__", [](const Class &self) { return self.repr(true, self.PY_NAMESPACE_SEPARATOR); })
            .def("__str__", &Class::str);
}

template <typename T>
void declare_image_array(py::module &m, std::string typestr) {
    using Class = gauss2d::ImageArray<T, PyImage<T>>;
    std::string pyclass_name = std::string("ImageArray") + typestr;
    py::class_<Class, std::shared_ptr<Class>>(m, pyclass_name.c_str())
            .def(py::init<const typename Class::Data *>(), "data"_a)
            .def("at", &Class::at, py::return_value_policy::reference)
            .def_property_readonly("size", &Class::size)
            .def("__repr__", [](const Class &self) { return self.repr(true, self.PY_NAMESPACE_SEPARATOR); })
            .def("__str__", &Class::str);
}

template <typename t>
void declare_evaluator(py::module &m, std::string typestr) {
    using Class = gauss2d::GaussianEvaluator<t, PyImage<t>, PyImage<gauss2d::idx_type>>;
    std::string pyclass_name = std::string("GaussianEvaluator") + typestr;
    py::class_<Class, std::shared_ptr<Class>>(m, pyclass_name.c_str())
            .def(py::init<const std::shared_ptr<const gauss2d::ConvolvedGaussians>,
                          const std::shared_ptr<const gauss2d::CoordinateSystem>,
                          const std::shared_ptr<const PyImage<t>>, const std::shared_ptr<const PyImage<t>>,
                          const std::shared_ptr<PyImage<t>>, const std::shared_ptr<PyImage<t>>,
                          const std::shared_ptr<gauss2d::ImageArray<t, PyImage<t>>>,
                          const std::shared_ptr<const PyImage<gauss2d::idx_type>>,
                          const std::shared_ptr<const PyImage<t>>,
                          const std::shared_ptr<const PyImage<gauss2d::idx_type>>,
                          const std::shared_ptr<const PyImage<t>>, const std::shared_ptr<const PyImage<t>>>(),
                 "gaussians"_a, "coordsys"_a = nullptr, "data"_a = nullptr, "sigma_inv"_a = nullptr,
                 "output"_a = nullptr, "residual"_a = nullptr, "grads"_a = nullptr,
                 "grad_param_map"_a = nullptr, "grad_param_factor"_a = nullptr, "extra_param_map"_a = nullptr,
                 "extra_param_factor"_a = nullptr, "background"_a = nullptr)
            .def("loglike_pixel", &Class::loglike_pixel, "to_add"_a = false)
            .def_property_readonly("n_cols", &Class::get_n_cols)
            .def_property_readonly("n_rows", &Class::get_n_rows)
            .def_property_readonly("n_size", &Class::get_size);
}

template <typename t, class Data, class Indices>
void declare_maker(py::module &m, std::string typestr) {
    m.def(("make_gaussians_pixel_" + typestr).c_str(), gauss2d::make_gaussians_pixel<t, Data, Indices>,
          "Evaluate a 2D Gaussian at the centers of pixels on a rectangular grid using the standard bivariate"
          "Gaussian PDF.",
          "gaussians"_a, "output"_a = nullptr, "n_rows"_a = 0, "n_cols"_a = 0, "coordsys"_a = nullptr,
          "to_add"_a = false);
}

}  // namespace lsst::gauss2d::python

#endif  // LSST_GAUSS2D_PYTHON_PYIMAGE_H
