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

#include "gauss2d/gaussian.h"
#ifndef __GAUSS2D_EVALUATE_H_
#include <gauss2d/evaluate.h>
#endif

#ifndef __GAUSS2D_GAUSSIAN_H_
#include <gauss2d/gaussian.h>
#endif

#ifndef __GAUSS2D_IMAGE_H_
#include "gauss2d/image.h"
#endif

/*#ifndef GAUSS2D_GAUSSIAN_INTEGRATOR_H
#include "../src/gaussian_integrator.h"
#endif*/

#include <memory>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <stdexcept>

namespace py = pybind11;
using namespace pybind11::literals;

#pragma GCC visibility push(hidden)
template <typename t>
class ImagePy : public gauss2d::Image<t, ImagePy<t>>
{
private:
    const size_t _n_rows;
    const size_t _n_cols;

    std::unique_ptr<py::array_t<t>> _ptr_own = nullptr;
    py::array_t<t> _data;
    py::detail::unchecked_mutable_reference<t, 2> _data_ref;

    void _validate() const {
        if((_ptr_own == nullptr ? _data : *_ptr_own).ndim() != 2) throw std::runtime_error("Input data must have 2 dimensions");
    }

public:
    inline t & _get_value_unchecked(size_t row, size_t col) { return this->_data_ref(row, col); };

    py::array_t<t> & get_data() {
        py::array_t<t> & rv = _ptr_own == nullptr ? _data : *_ptr_own;
        return rv;
    }

    size_t get_n_cols() const { return _n_cols;};
    size_t get_n_rows() const { return _n_rows;};

    //void add_value(size_t row, size_t col, t value) { this->_get_value(row, col) += value;}
    //void add_value_unchecked(size_t row, size_t col, t value) {
    //    _get_value_unchecked(row, col) += value;
    //}
    inline t get_value_unchecked(size_t row, size_t col) const  { return this->_data_ref(row, col);};
    void set_value(size_t row, size_t col, t value) { this->_data_ref(row, col) = value;}
    //void set_value_unchecked(size_t row, size_t col, t value) { _get_value_unchecked(row, col) = value;};

    ImagePy(size_t n_rows, size_t n_cols, const std::shared_ptr<const gauss2d::CoordinateSystem> coordsys = nullptr) :
        gauss2d::Image<t, ImagePy<t>>(coordsys),
        _n_rows(n_rows), _n_cols(n_cols),
        _ptr_own(std::make_unique<py::array_t<t>>(py::array::ShapeContainer({n_rows, n_cols}))),
        _data_ref(this->get_data().template mutable_unchecked<2>())
    {
        _validate();
    }
    ImagePy(py::array_t<t> data, const std::shared_ptr<const gauss2d::CoordinateSystem> coordsys = nullptr) : 
        gauss2d::Image<t, ImagePy<t>>(coordsys),
        _n_rows(data.shape(0)), _n_cols(data.shape(1)),
        _data(data),
        _data_ref(this->get_data().template mutable_unchecked<2>())
    {
        _validate();
    }

    ~ImagePy() {};
};
#pragma GCC visibility pop

template<typename T>
void declare_image(py::module &m, std::string typestr) {
    using Class = ImagePy<T>;
    std::string pyclass_name = std::string("ImagePy") + typestr;
    py::class_<Class, std::shared_ptr<Class>>(m, pyclass_name.c_str())
    .def(
        py::init<size_t, size_t, const std::shared_ptr<const gauss2d::CoordinateSystem>>(),
        "dim_y"_a, "dim_x"_a, "coordsys"_a = gauss2d::COORDS_DEFAULT
    )
    .def(
        py::init<py::array_t<T>, const std::shared_ptr<const gauss2d::CoordinateSystem>>(),
        "data"_a, "coordsys"_a = gauss2d::COORDS_DEFAULT
    )
    .def_property_readonly("n_rows", &Class::get_n_rows)
    .def_property_readonly("n_cols", &Class::get_n_cols)
    .def_property_readonly("data", &Class::get_data)
    .def("get_value", &Class::get_value)
    .def("set_value", &Class::set_value)
    .def("get_value_unchecked", &Class::get_value_unchecked)
    .def("set_value_unchecked", &Class::set_value_unchecked)
    .def_property_readonly("size", &Class::size);
}

template<typename T>
void declare_image_array(py::module &m, std::string typestr) {
    using Class = gauss2d::ImageArray<T, ImagePy<T>>;
    std::string pyclass_name = std::string("ImageArrayPy") + typestr;
    py::class_<Class, std::shared_ptr<Class>>(m, pyclass_name.c_str())
    .def(py::init<const typename Class::Data*>(), "data"_a)
    .def("at", &Class::at, py::return_value_policy::reference)
    .def_property_readonly("size", &Class::size);
}

template<typename t>
void declare_evaluator(py::module &m, std::string typestr) {
    using Class = gauss2d::GaussianEvaluator<t, ImagePy<t>, ImagePy<gauss2d::idx_type>>;
    std::string pyclass_name = std::string("GaussianEvaluatorPy") + typestr;
    py::class_<Class, std::shared_ptr<Class>>(m, pyclass_name.c_str())
    .def(
        py::init<
            const std::shared_ptr<const gauss2d::Gaussians>,
            const std::shared_ptr<const gauss2d::CoordinateSystem>,
            const std::shared_ptr<const ImagePy<t>>,
            const std::shared_ptr<const ImagePy<t>>,
            const std::shared_ptr<ImagePy<t>>,
            const std::shared_ptr<ImagePy<t>>,
            const std::shared_ptr<gauss2d::ImageArray<t, ImagePy<t>>>,
            const std::shared_ptr<const ImagePy<gauss2d::idx_type>>,
            const std::shared_ptr<const ImagePy<t>>,
            const std::shared_ptr<const ImagePy<gauss2d::idx_type>>,
            const std::shared_ptr<const ImagePy<t>>,
            const std::shared_ptr<const ImagePy<t>>
        >(),
        "gaussians"_a,
        "coordsys"_a = nullptr,
        "data"_a = nullptr,
        "sigma_inv"_a = nullptr,
        "output"_a = nullptr,
        "residual"_a = nullptr,
        "grads"_a = nullptr,
        "grad_param_map"_a = nullptr, 
        "grad_param_factor"_a = nullptr,
        "extra_param_map"_a = nullptr,
        "extra_param_factor"_a = nullptr,
        "background"_a = nullptr
    )
    .def("loglike_pixel", &Class::loglike_pixel, "to_add"_a = false)
    ;
}

template<typename t, class Data, class Indices>
void declare_maker(py::module &m, std::string typestr) {
    m.def(
        ("make_gaussian_pixel_py_" + typestr).c_str(), gauss2d::make_gaussians_pixel<t, Data, Indices>,
        "Evaluate a 2D Gaussian at the centers of pixels on a rectangular grid using the standard bivariate"
        "Gaussian PDF.",
        "gaussians"_a, "output"_a=nullptr, "n_rows"_a=0, "n_cols"_a=0, "coordsys"_a=nullptr
    );
}

PYBIND11_MODULE(_gauss2d, m)
{
    m.doc() = "Gauss2D Python bindings";
    m.attr("M_HWHM_SIGMA") = py::float_(gauss2d::M_HWHM_SIGMA);
    m.attr("M_SIGMA_HWHM") = py::float_(gauss2d::M_SIGMA_HWHM);

    py::class_<gauss2d::Centroid, std::shared_ptr<gauss2d::Centroid>>(m, "Centroid")
        .def(py::init<double, double>(), "x"_a=0, "y"_a=0);
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
        .def("__repr__", &gauss2d::Covariance::str)
    ;
    py::class_<gauss2d::Ellipse, std::shared_ptr<gauss2d::Ellipse> >(m, "Ellipse")
        .def(py::init<gauss2d::Covariance&>())
        .def(py::init<gauss2d::EllipseMajor&>())
        .def(py::init<double, double, double>(), "sigma_x"_a=0, "sigma_y"_a=0, "rho"_a=0)
        .def("_set", &gauss2d::Ellipse::_set)
        .def("_set_rho", &gauss2d::Ellipse::_set_rho)
        .def("_set_sigma_x", &gauss2d::Ellipse::_set_sigma_x)
        .def("_set_sigma_y", &gauss2d::Ellipse::_set_sigma_y)
        .def_static("check", &gauss2d::Ellipse::check)
        .def("convolve", &gauss2d::Ellipse::convolve)
        .def("get_cov_xy", &gauss2d::Ellipse::get_cov_xy)
        .def("get_radius_trace", &gauss2d::Ellipse::get_radius_trace)
        .def("make_convolution", &gauss2d::Ellipse::make_convolution)

        .def_property("rho", &gauss2d::Ellipse::get_rho, &gauss2d::Ellipse::set_rho)
        .def_property("sigma_x", &gauss2d::Ellipse::get_sigma_x, &gauss2d::Ellipse::set_sigma_x)
        .def_property("sigma_y", &gauss2d::Ellipse::get_sigma_y, &gauss2d::Ellipse::set_sigma_y)
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
        .def("__repr__", &gauss2d::EllipseMajor::str)
    ;

    py::class_<gauss2d::Gaussian, std::shared_ptr<gauss2d::Gaussian>>(m, "Gaussian")
        .def(py::init<std::shared_ptr<gauss2d::Centroid>, std::shared_ptr<gauss2d::Ellipse>, double>(),
            "centroid"_a, "ellipse"_a, "integral"_a=1.)
        .def_property("const_normal", &gauss2d::Gaussian::get_const_normal, &gauss2d::Gaussian::set_const_normal)
        .def_property("integral", &gauss2d::Gaussian::get_integral, &gauss2d::Gaussian::set_integral)
        .def("__repr__", &gauss2d::Gaussian::str)
    ;
    py::class_<gauss2d::ConvolvedGaussian, std::shared_ptr<gauss2d::ConvolvedGaussian>>(m, "ConvolvedGaussian")
        .def(py::init<std::shared_ptr<gauss2d::Gaussian>, std::shared_ptr<gauss2d::Gaussian>>(),
            "source"_a = nullptr, "kernel"_a = nullptr)
        .def("__repr__", &gauss2d::ConvolvedGaussian::str)
    ;
    py::class_<gauss2d::Gaussians, std::shared_ptr<gauss2d::Gaussians>>(m, "Gaussians")
        .def(py::init<const gauss2d::Gaussians::Data *>(), "gaussians"_a)
        .def("__repr__", &gauss2d::Gaussians::str)
    ;
    py::class_<gauss2d::CoordinateSystem, std::shared_ptr<gauss2d::CoordinateSystem>>(m, "CoordinateSystem")
        .def(py::init<double, double>(), "d1"_a=1., "d2"_a=1.)
    ;
    declare_image<float>(m, "F");
    declare_image<double>(m, "D");
    declare_image<int>(m, "I");
    declare_image<unsigned int>(m, "U");
    declare_image<size_t>(m, "S");
    declare_image_array<float>(m, "F");
    declare_image_array<double>(m, "D");
    declare_image_array<int>(m, "I");
    declare_image_array<unsigned int>(m, "U");
    declare_image_array<size_t>(m, "S");
    declare_evaluator<float>(m, "F");
    declare_evaluator<double>(m, "D");
    declare_maker<float, ImagePy<float>, ImagePy<size_t>>(m, "F");
    declare_maker<double, ImagePy<double>, ImagePy<size_t>>(m, "D");
}

