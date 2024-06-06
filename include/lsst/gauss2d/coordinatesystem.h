// -*- LSST-C++ -*-
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
#define LSST_GAUSS2D_COORDINATESYSTEM_H

#include <string>

#include "object.h"
#include "type_name.h"

namespace lsst::gauss2d {

/**
 * A coordinate system specifying image scale and orientation.
 *
 * This is intended to mimic some of the functionality of e.g. basic
 * FITS headers, in order for evaluators to draw images at different
 * scales and/or with rotation/translation. Rotation is not yet supported.
 *
 **/
class CoordinateSystem : public Object {
private:
    double _dx1;
    double _dy2;

    double _x_min;
    double _y_min;

public:
    // x,y coords of the bottom left corner of the image if not rotated
    // TBD how to implement this if there is rotation...

    double get_dx1() const;
    double get_dy2() const;

    double get_x_min() const;
    double get_y_min() const;

    bool is_xy_aligned() const;

    bool operator==(const CoordinateSystem& other) const;
    bool operator!=(const CoordinateSystem& other) const;

    std::string repr(bool name_keywords = false,
                     std::string_view namespace_separator = Object::CC_NAMESPACE_SEPARATOR) const override;
    std::string str() const override;

    CoordinateSystem(double dx1 = 1., double dy2 = 1, double x_min = 0, double y_min = 0);
    ~CoordinateSystem();
};

static const CoordinateSystem COORDS_DEFAULT{};

}  // namespace lsst::gauss2d
#endif  // LSST_GAUSS2D_COORDINATESYSTEM_H
