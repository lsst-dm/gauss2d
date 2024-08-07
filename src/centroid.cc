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

#ifndef LSST_GAUSS2D_CENTROID_H

#include <memory>
#include <string>

#include "lsst/gauss2d/centroid.h"
#include "lsst/gauss2d/to_string.h"
#include "lsst/gauss2d/type_name.h"

namespace lsst::gauss2d {

double CentroidValues::get_x() const { return *_x; }
std::array<double, 2> CentroidValues::get_xy() const { return {*_x, *_y}; }
double CentroidValues::get_y() const { return *_y; }

void CentroidValues::set_x(double x) { *_x = x; }
void CentroidValues::set_xy(const std::array<double, 2>& xy) {
    *_x = xy[0];
    *_y = xy[1];
};
void CentroidValues::set_y(double y) { *_y = y; }

std::string CentroidValues::repr(bool name_keywords, std::string_view namespace_separator) const {
    return type_name_str<CentroidValues>(false, namespace_separator) + "(" + (name_keywords ? "x=" : "")
           + to_string_float(*_x) + ", " + (name_keywords ? "y=" : "") + to_string_float(*_y) + ")";
}

std::string CentroidValues::str() const {
    return type_name_str<CentroidValues>(true) + "(x=" + to_string_float(*_x) + ", y=" + to_string_float(*_y)
           + ")";
}

bool CentroidValues::operator==(const CentroidValues& other) const {
    return (this->get_x() == other.get_x()) && (this->get_y() == other.get_y());
};

bool CentroidValues::operator!=(const CentroidValues& other) const { return !(*this == other); }

CentroidValues::CentroidValues(std::shared_ptr<double> x, std::shared_ptr<double> y)
        : _x(x == nullptr ? std::make_shared<double>(0) : std::move(x)),
          _y(y == nullptr ? std::make_shared<double>(0) : std::move(y)) {}

CentroidValues::CentroidValues(double x, double y)
        : _x(std::make_shared<double>(x)), _y(std::make_shared<double>(y)) {}

void Centroid::convolve(const Centroid& cen) {
    this->set_x(this->get_x() + cen.get_x());
    this->set_y(this->get_y() + cen.get_y());
}

const CentroidData& Centroid::get_data() const { return *_data; }

double Centroid::get_x() const { return _data->get_x(); }

std::array<double, 2> Centroid::get_xy() const { return _data->get_xy(); }

double Centroid::get_y() const { return _data->get_y(); }

std::shared_ptr<Centroid> Centroid::make_convolution(const Centroid& cen) const {
    return this->make_convolution_uniq(cen);
}
std::unique_ptr<Centroid> Centroid::make_convolution_uniq(const Centroid& cen) const {
    // TODO: Replace with cloning derived data
    std::unique_ptr<Centroid> cen_ret = std::make_unique<Centroid>(this->get_x(), this->get_y());
    cen_ret->convolve(cen);
    return cen_ret;
}

void Centroid::set_x(double x) { _data->set_x(x); }
void Centroid::set_xy(const std::array<double, 2>& xy) { _data->set_xy(xy); }
void Centroid::set_y(double y) { _data->set_y(y); }

std::string Centroid::repr(bool name_keywords, std::string_view namespace_separator) const {
    return type_name_str<Centroid>(false, namespace_separator) + "(" + (name_keywords ? "data=" : "")
           + _data->repr(name_keywords, namespace_separator) + ")";
}

std::string Centroid::str() const { return type_name_str<Centroid>(true) + "(data=" + _data->str() + ")"; }

bool Centroid::operator==(const Centroid& other) const {
    return (this->get_x() == other.get_x()) && (this->get_y() == other.get_y());
};

bool Centroid::operator!=(const Centroid& other) const { return !(*this == other); }

Centroid::Centroid(std::shared_ptr<CentroidData> data)
        : _data(data == nullptr ? std::make_shared<CentroidValues>() : std::move(data)) {}

Centroid::Centroid(double x, double y) : _data(std::make_shared<CentroidValues>()) {
    set_x(x);
    set_y(y);
}

}  // namespace lsst::gauss2d

#endif
