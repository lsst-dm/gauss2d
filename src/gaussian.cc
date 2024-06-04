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

#ifndef LSST_GAUSS2D_GAUSSIAN_H

#include <optional>
#include <stdexcept>

#include "lsst/gauss2d/gaussian.h"
#include "lsst/gauss2d/type_name.h"

namespace lsst::gauss2d {
std::string GaussianIntegralValue::repr(bool name_keywords, std::string_view namespace_separator) const {
    return type_name_str<GaussianIntegralValue>(false, namespace_separator) + "("
           + (name_keywords ? "value=" : "") + std::to_string(*_value) + ")";
}
std::string GaussianIntegralValue::str() const {
    return type_name_str<GaussianIntegralValue>(true) + "(value=" + std::to_string(*_value) + ")";
}

GaussianIntegralValue::GaussianIntegralValue(double value) : _value(std::make_shared<double>(value)){};
GaussianIntegralValue::GaussianIntegralValue(std::shared_ptr<double> value)
        : _value(value == nullptr ? std::make_shared<double>(1) : std::move(value)){};

double Gaussian::get_const_normal() const { return _integral->get_value() / (2 * _ellipse->get_area()); }
double Gaussian::get_integral_value() const { return _integral->get_value(); };

Centroid& Gaussian::get_centroid() { return *_centroid; }
Ellipse& Gaussian::get_ellipse() { return *_ellipse; }
GaussianIntegral& Gaussian::get_integral() { return *_integral; }

std::shared_ptr<Centroid> Gaussian::get_centroid_ptr() { return _centroid; }
std::shared_ptr<Ellipse> Gaussian::get_ellipse_ptr() { return _ellipse; }
std::shared_ptr<GaussianIntegral> Gaussian::get_integral_ptr() { return _integral; }

const Centroid& Gaussian::get_centroid_const() const { return *_centroid; }
const Ellipse& Gaussian::get_ellipse_const() const { return *_ellipse; }
const GaussianIntegral& Gaussian::get_integral_const() const { return *_integral; }

void Gaussian::set_const_normal(double const_normal) {
    _integral->set_value(get_const_normal() * 2 * _ellipse->get_area());
}
void Gaussian::set_integral_value(double integral) { _integral->set_value(integral); }

void Gaussian::set_centroid_ptr(std::shared_ptr<Centroid> centroid) {
    _centroid = this->_check_not_nullptr<Centroid>(centroid, "centroid");
}
void Gaussian::set_ellipse_ptr(std::shared_ptr<Ellipse> ellipse) {
    _ellipse = this->_check_not_nullptr<Ellipse>(ellipse, "ellipse");
}
void Gaussian::set_integral_ptr(std::shared_ptr<GaussianIntegral> integral) {
    _integral = this->_check_not_nullptr<GaussianIntegral>(integral, "integral");
}

std::string Gaussian::repr(bool name_keywords, std::string_view namespace_separator) const {
    return type_name_str<Gaussian>(false, namespace_separator) + "(" + (name_keywords ? "centroid=" : "")
           + _centroid->str() + ", " + (name_keywords ? "ellipse=" : "") + _ellipse->str() + ", "
           + (name_keywords ? "integral=" : "") + _integral->str() + ")";
}

std::string Gaussian::str() const {
    return type_name_str<Gaussian>(true) + "(centroid=" + _centroid->str() + ", ellipse=" + _ellipse->str()
           + ", integral=" + _integral->str() + ")";
}

bool Gaussian::operator==(const Gaussian& other) const {
    return (*_centroid == other.get_centroid_const()) && (*_ellipse == other.get_ellipse_const())
           && (*_integral == other.get_integral_const());
}

bool Gaussian::operator!=(const Gaussian& other) const { return !(*this == other); }

std::ostream& operator<<(std::ostream& out, const Gaussian& g) {
    out << g.str();
    return out;
}

Gaussian::Gaussian(std::shared_ptr<Centroid> centroid, std::shared_ptr<Ellipse> ellipse,
                   std::shared_ptr<GaussianIntegral> integral)
        : _centroid(centroid != nullptr ? std::move(centroid) : std::make_shared<Centroid>()),
          _ellipse(ellipse != nullptr ? std::move(ellipse) : std::make_shared<Ellipse>()),
          _integral(integral != nullptr ? std::move(integral) : std::make_shared<GaussianIntegralValue>()) {}
Gaussian::~Gaussian(){};

Gaussian& Gaussians::operator[](size_t i) { return *(_data[i]); }
const Gaussian& Gaussians::operator[](size_t i) const { return *(_data[i]); }

size_t Gaussians::assign(const Data& data, size_t i) {
    const size_t n_data = this->_data.size();
    const size_t i_max = i + data.size();
    if (!(i_max <= n_data)) {
        throw std::out_of_range("data i_max=" + std::to_string(i_max) + ">= this._data.size()="
                                + std::to_string(n_data) + " for this=" + this->str());
    }
    size_t i_begin = i;
    for (const auto& gauss : data) {
        if (gauss == nullptr)
            throw std::runtime_error("Gaussians data[" + std::to_string(i - i_begin) + "] can't be null");
        _data[i++] = std::move(gauss);
    }
    return i;
}

Gaussian& Gaussians::at(size_t i) const { return *(_data.at(i)); }
const Gaussian& Gaussians::at_const(size_t i) const { return *(_data.at(i)); }

typename Gaussians::Data::iterator Gaussians::begin() noexcept { return _data.begin(); }
typename Gaussians::Data::const_iterator Gaussians::begin() const noexcept { return _data.cbegin(); }
typename Gaussians::Data::const_iterator Gaussians::cbegin() const noexcept { return _data.cbegin(); }

typename Gaussians::Data::iterator Gaussians::end() noexcept { return _data.end(); }
typename Gaussians::Data::const_iterator Gaussians::end() const noexcept { return _data.cend(); }
typename Gaussians::Data::const_iterator Gaussians::cend() const noexcept { return _data.cend(); }

Gaussians::Data Gaussians::get_data() const { return _data; }

size_t Gaussians::size() const { return _data.size(); }

std::string Gaussians::repr(bool name_keywords, std::string_view namespace_separator) const {
    std::string str
            = type_name_str<Gaussians>(false, namespace_separator) + "(" + (name_keywords ? "data=" : "");
    str += repr_iter_ptr(_data, name_keywords, namespace_separator);
    return str + ")";
}

std::string Gaussians::str() const {
    std::string str = type_name_str<Gaussians>(true) + "(data=";
    str += str_iter_ptr(_data);
    return str + ")";
}

Gaussians::Gaussians(std::optional<const Data> data) {
    if (data) {
        size_t n_data = data->size();
        if (n_data > 0) {
            _data.resize(n_data, nullptr);
            this->assign(*data);
        }
    }
}
Gaussians::Gaussians(std::vector<std::optional<const Data>> data) {
    size_t n_data = 0;
    for (const auto& datum : data) {
        if (datum) n_data += datum->size();
    }
    if (n_data > 0) {
        size_t i = 0;
        _data.resize(n_data);
        for (const auto& datum : data) i = this->assign(*datum, i);
    }
}

const Gaussian& ConvolvedGaussian::get_source() const { return *_source; }
const Gaussian& ConvolvedGaussian::get_kernel() const { return *_kernel; }

std::unique_ptr<Gaussian> ConvolvedGaussian::make_convolution() const {
    return std::make_unique<Gaussian>(
            _source->get_centroid_const().make_convolution(_kernel->get_centroid_const()),
            _source->get_ellipse_const().make_convolution(_kernel->get_ellipse_const()),
            std::make_shared<GaussianIntegralValue>(_source->get_integral_value()
                                                    + _kernel->get_integral_value()));
}

std::string ConvolvedGaussian::repr(bool name_keywords, std::string_view namespace_separator) const {
    return type_name_str<ConvolvedGaussian>(false, namespace_separator) + "("
           + (name_keywords ? "source=" : "") + _source->repr(name_keywords, namespace_separator) + ", "
           + (name_keywords ? "kernel=" : "") + _kernel->repr(name_keywords, namespace_separator) + ")";
}

std::string ConvolvedGaussian::str() const {
    return type_name_str<ConvolvedGaussian>(true) + "(source=" + _source->str() + ", kernel=" + _kernel->str()
           + ")";
}

static const std::shared_ptr<const Gaussian> GAUSS_ZERO = std::make_shared<const Gaussian>();

ConvolvedGaussian::ConvolvedGaussian(std::shared_ptr<const Gaussian> source,
                                     std::shared_ptr<const Gaussian> kernel)
        : _source(source != nullptr ? source : GAUSS_ZERO),
          _kernel(kernel != nullptr ? kernel : GAUSS_ZERO) {}

size_t ConvolvedGaussians::assign(const Data& data, size_t i) {
    const size_t n_data = this->_data.size();
    const size_t i_max = i + data.size();
    if (!(i_max <= n_data)) {
        throw std::out_of_range("data i_max=" + std::to_string(i_max) + ">= this._data.size()="
                                + std::to_string(n_data) + " for this=" + this->str());
    }
    size_t i_begin = i;
    for (const auto& gauss : data) {
        if (gauss == nullptr)
            throw std::runtime_error("ConvolvedGaussian data[" + std::to_string(i - i_begin)
                                     + "] can't be null");
        _data[i++] = std::move(gauss);
    }
    return i;
}

ConvolvedGaussian& ConvolvedGaussians::at(size_t i) const { return *(_data.at(i)); }
const ConvolvedGaussian& ConvolvedGaussians::at_const(size_t i) const { return *(_data.at(i)); }

typename ConvolvedGaussians::Data::iterator ConvolvedGaussians::begin() noexcept { return _data.begin(); }
typename ConvolvedGaussians::Data::iterator ConvolvedGaussians::end() noexcept { return _data.end(); }

typename ConvolvedGaussians::Data::const_iterator ConvolvedGaussians::begin() const noexcept {
    return _data.begin();
}
typename ConvolvedGaussians::Data::const_iterator ConvolvedGaussians::end() const noexcept {
    return _data.cend();
}

typename ConvolvedGaussians::Data::const_iterator ConvolvedGaussians::cbegin() const noexcept {
    return _data.begin();
}
typename ConvolvedGaussians::Data::const_iterator ConvolvedGaussians::cend() const noexcept {
    return _data.cend();
}

ConvolvedGaussians::Data ConvolvedGaussians::get_data() const { return _data; }

size_t ConvolvedGaussians::size() const { return _data.size(); }

std::string ConvolvedGaussians::repr(bool name_keywords, std::string_view namespace_separator) const {
    std::string str = type_name_str<ConvolvedGaussians>(false, namespace_separator) + "("
                      + (name_keywords ? "data=" : "");
    str += repr_iter_ptr(_data, name_keywords, namespace_separator);
    return str + ")";
}

std::string ConvolvedGaussians::str() const {
    std::string str = type_name_str<ConvolvedGaussians>(true) + "(data=";
    str += str_iter_ptr(_data);
    return str + ")";
}

ConvolvedGaussian& ConvolvedGaussians::operator[](size_t i) { return *(_data[i]); }
const ConvolvedGaussian& ConvolvedGaussians::operator[](size_t i) const { return *(_data[i]); }

ConvolvedGaussians::ConvolvedGaussians(std::optional<const Data> data) {
    if (data) {
        size_t n_data = data->size();
        if (n_data > 0) {
            _data.resize(n_data);
            this->assign(*data);
        }
    }
}
ConvolvedGaussians::ConvolvedGaussians(std::vector<std::optional<const Data>> data) {
    size_t n_data = 0;
    for (const auto& datum : data) {
        if (datum) n_data += datum->size();
    }
    if (n_data > 0) {
        size_t i = 0;
        _data.resize(n_data);
        for (const auto& datum : data) i = this->assign(*datum, i);
    }
}
}  // namespace lsst::gauss2d

#endif
