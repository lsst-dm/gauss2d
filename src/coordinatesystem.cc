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

#ifndef LSST_GAUSS2D_COORDINATESYSTEM_H

#include <cmath>
#include <stdexcept>
#include <string>

#include "lsst/gauss2d/coordinatesystem.h"
#include "lsst/gauss2d/type_name.h"
#include "lsst/gauss2d/to_string.h"

namespace lsst::gauss2d {

double CoordinateSystem::get_dx1() const { return _dx1; }
double CoordinateSystem::get_dy2() const { return _dy2; }

double CoordinateSystem::get_x_min() const { return _x_min; }
double CoordinateSystem::get_y_min() const { return _y_min; }

bool CoordinateSystem::is_xy_aligned() const { return true; }

bool CoordinateSystem::operator==(const CoordinateSystem& other) const {
    return (_dx1 == other.get_dx1()) && (_dy2 == other.get_dy2()) && (_x_min == other.get_x_min())
           && (_y_min == other.get_y_min());
}

bool CoordinateSystem::operator!=(const CoordinateSystem& other) const { return !(*this == other); }

std::string CoordinateSystem::repr(bool name_keywords, std::string_view namespace_separator) const {
    return type_name_str<CoordinateSystem>(false, namespace_separator) + "(" + (name_keywords ? "dx1=" : "")
           + to_string_float(_dx1) + ", " + (name_keywords ? "dy2=" : "") + to_string_float(_dy2) + ", "
           + (name_keywords ? "x_min=" : "") + to_string_float(_x_min) + ", "
           + (name_keywords ? "y_min=" : "") + to_string_float(_y_min) + ")";
}

std::string CoordinateSystem::str() const {
    return type_name_str<CoordinateSystem>(true) + "(dx1=" + to_string_float(_dx1)
           + ", dy2=" + to_string_float(_dy2) + ", x_min=" + to_string_float(_x_min)
           + ", y_min=" + to_string_float(_y_min) + ")";
}

CoordinateSystem::CoordinateSystem(double dx1, double dy2, double x_min, double y_min)
        : _dx1(dx1), _dy2(dy2), _x_min(x_min), _y_min(y_min) {
    std::string errmsg = "";
    if (!((dx1 > 0) && std::isfinite(dx1) && (dy2 > 0) && std::isfinite(dy2))) {
        errmsg += "dx1, dy2 = " + to_string_float(_dx1) + ", " + to_string_float(_dy2) + " !>0 or !finite; ";
    }
    if (!(std::isfinite(_x_min) && std::isfinite(_y_min))) {
        errmsg += "x_min, y_min = " + to_string_float(_x_min) + ", " + to_string_float(_y_min) + " !finite; ";
    }
    if (errmsg != "") {
        throw std::invalid_argument(errmsg);
    }
}

CoordinateSystem::~CoordinateSystem() {}

}  // namespace lsst::gauss2d

#endif
