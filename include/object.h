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

#ifndef GAUSS2D_OBJECT_H
#define GAUSS2D_OBJECT_H

#include <string>

namespace gauss2d {

class Object {
public:
    virtual std::string repr(bool name_keywords = false) const = 0;
    virtual std::string str() const = 0;

    friend std::ostream &operator<<(std::ostream &out, const Object &obj) {
        out << obj.str();
        return out;
    }

    virtual ~Object() = default;
};

}  // namespace gauss2d
#endif
