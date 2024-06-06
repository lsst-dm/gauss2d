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
#include "lsst/gauss2d/string_utils.h"
#include "lsst/gauss2d/type_name.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace gauss2d = lsst::gauss2d;

namespace lsst::gauss2d::python {

/*
    Convenience functions for binding concrete image types.
*/

template <typename T>
std::string replace_type(std::string target, std::string replacement) {
    std::string token = "<" + gauss2d::type_name_str<T>() + ">";
    return replace_all(target, token, replacement);
}

std::string_view STR_NONE = "";

#pragma GCC visibility push(hidden)
template <typename T>
class Image : public lsst::gauss2d::Image<T, Image<T>> {
private:
    const size_t _n_rows;
    const size_t _n_cols;

    std::unique_ptr<py::array_t<T>> _ptr_own = nullptr;
    py::array_t<T> _data;
    py::detail::unchecked_mutable_reference<T, 2> _data_ref;

    void _validate() const {
        if ((_ptr_own == nullptr ? _data : *_ptr_own).ndim() != 2) {
            throw std::runtime_error("Input data must have 2 dimensions");
        }
    }

public:
    inline T &_get_value_unchecked(size_t row, size_t col) { return this->_data_ref(row, col); };

    py::array_t<T> &get_data() {
        py::array_t<T> &rv = _ptr_own == nullptr ? _data : *_ptr_own;
        return rv;
    }

    size_t get_n_cols() const { return _n_cols; };
    size_t get_n_rows() const { return _n_rows; };

    // void add_value(size_t row, size_t col, t value) { this->_get_value(row, col) += value;}
    // void add_value_unchecked(size_t row, size_t col, t value) {
    //     _get_value_unchecked(row, col) += value;
    // }
    const inline T get_value_unchecked(size_t row, size_t col) const { return this->_data_ref(row, col); };
    void set_value_unchecked(size_t row, size_t col, T value) { this->_data_ref(row, col) = value; }
    // void set_value_unchecked(size_t row, size_t col, t value) { _get_value_unchecked(row, col) = value;};

    Image(size_t n_rows, size_t n_cols,
          const std::shared_ptr<const gauss2d::CoordinateSystem> coordsys = nullptr)
            : gauss2d::Image<T, Image<T>>(coordsys),
              _n_rows(n_rows),
              _n_cols(n_cols),
              _ptr_own(std::make_unique<py::array_t<T>>(py::array::ShapeContainer({n_rows, n_cols}))),
              _data_ref(this->get_data().template mutable_unchecked<2>()) {
        _validate();
    }
    Image(py::array_t<T> data, const std::shared_ptr<const gauss2d::CoordinateSystem> coordsys = nullptr)
            : gauss2d::Image<T, Image<T>>(coordsys),
              _n_rows(data.shape(0)),
              _n_cols(data.shape(1)),
              _data(data),
              _data_ref(this->get_data().template mutable_unchecked<2>()) {
        _validate();
    }

    ~Image(){};
};
#pragma GCC visibility pop

template <typename T>
void declare_image(py::module &m, std::string str_type) {
    using Class = Image<T>;
    std::string pyclass_name = std::string("Image") + str_type;
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
            .def("get_value", &Class::get_value, "row"_a, "col"_a)
            .def("set_value", &Class::set_value, "row"_a, "col"_a, "value"_a)
            .def("get_value_unchecked", &Class::get_value_unchecked, "row"_a, "col"_a)
            .def("set_value_unchecked", &Class::set_value_unchecked, "row"_a, "col"_a, "value"_a)
            .def_property_readonly("size", &Class::size)
            .def_property_readonly("shape", &Class::shape)
            .def(py::self == py::self)
            .def(py::self != py::self)
            .def(py::self += T())
            .def(py::self *= T())
            .def("__repr__",
                 [str_type](const Class &self) {
                     std::string repr = self.repr(true, self.PY_NAMESPACE_SEPARATOR);
                     return replace_type<T>(repr, str_type);
                 })
            .def("__str__", [str_type](const Class &self) {
                std::string str = self.str();
                return replace_type<T>(str, str_type);
            });
}

template <typename T>
std::string replace_image_types(std::string target, std::string str_type, std::string_view separator) {
    std::string token1 = std::string("<") + type_name_str<T>() + ", ";
    target = replace_all(target, token1, STR_NONE);
    std::string token2 = type_name_str<Image<T>>(false, separator) + std::string(" >");
    target = replace_all(target, token2, str_type);
    // str appears to return e.g. Image<double> for some reason
    target = replace_type<T>(target, str_type);
    return target;
}

template <typename T>
void declare_image_array(py::module &m, std::string str_type) {
    using Class = gauss2d::ImageArray<T, Image<T>>;
    std::string pyclass_name = std::string("ImageArray") + str_type;
    py::class_<Class, std::shared_ptr<Class>>(m, pyclass_name.c_str())
            .def(py::init<const typename Class::Data *>(), "data"_a)
            .def("at", &Class::at, py::return_value_policy::reference)
            .def_property_readonly("size", &Class::size)
            .def("__getitem__", &Class::at, py::return_value_policy::reference)
            .def("__len__", &Class::size)
            .def("__repr__",
                 [str_type](const Class &self) {
                     std::string repr = self.repr(true, self.PY_NAMESPACE_SEPARATOR);
                     repr = replace_image_types<T>(repr, str_type, self.PY_NAMESPACE_SEPARATOR);
                     return repr;
                 })
            .def("__str__", [str_type](const Class &self) {
                std::string str = self.str();
                str = replace_image_types<T>(str, str_type, self.CC_NAMESPACE_SEPARATOR);
                return str;
            });
}

template <typename T>
void declare_evaluator(py::module &m, std::string str_type) {
    using Class = gauss2d::GaussianEvaluator<T, Image<T>, Image<gauss2d::idx_type>>;
    std::string pyclass_name = std::string("GaussianEvaluator") + str_type;
    py::class_<Class, std::shared_ptr<Class>>(m, pyclass_name.c_str())
            .def(py::init<const std::shared_ptr<const gauss2d::ConvolvedGaussians>,
                          const std::shared_ptr<const Image<T>>, const std::shared_ptr<const Image<T>>,
                          const std::shared_ptr<Image<T>>, const std::shared_ptr<Image<T>>,
                          const std::shared_ptr<gauss2d::ImageArray<T, Image<T>>>,
                          const std::shared_ptr<const Image<gauss2d::idx_type>>,
                          const std::shared_ptr<const Image<T>>,
                          const std::shared_ptr<const Image<gauss2d::idx_type>>,
                          const std::shared_ptr<const Image<T>>, const std::shared_ptr<const Image<T>>>(),
                 "gaussians"_a, "data"_a = nullptr, "sigma_inv"_a = nullptr, "output"_a = nullptr,
                 "residual"_a = nullptr, "grads"_a = nullptr, "grad_param_map"_a = nullptr,
                 "grad_param_factor"_a = nullptr, "extra_param_map"_a = nullptr,
                 "extra_param_factor"_a = nullptr, "background"_a = nullptr)
            .def("loglike_pixel", &Class::loglike_pixel, "to_add"_a = false)
            .def_property_readonly("n_cols", &Class::get_n_cols)
            .def_property_readonly("n_rows", &Class::get_n_rows)
            .def_property_readonly("size", &Class::get_size)
            .def("__len__", &Class::get_size)
            .def("__repr__",
                 [str_type](const Class &self) {
                     std::string repr = self.repr(true, self.PY_NAMESPACE_SEPARATOR);
                     return repr;
                 })
            .def("__str__", [str_type](const Class &self) {
                std::string str = self.str();
                return str;
            });
}

template <typename T, class Data, class Indices>
void declare_maker(py::module &m, std::string str_type) {
    m.def(("make_gaussians_pixel_" + str_type).c_str(), gauss2d::make_gaussians_pixel<T, Data, Indices>,
          "Evaluate a 2D Gaussian at the centers of pixels on a rectangular grid using the standard bivariate"
          "Gaussian PDF.",
          "gaussians"_a, "output"_a = nullptr, "n_rows"_a = 0, "n_cols"_a = 0, "coordsys"_a = nullptr,
          "to_add"_a = false);
}

}  // namespace lsst::gauss2d::python

#endif  // LSST_GAUSS2D_PYTHON_PYIMAGE_H
