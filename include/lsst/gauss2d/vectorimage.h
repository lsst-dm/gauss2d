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

#ifndef LSST_GAUSS2D_VECTORIMAGE_H
#define LSST_GAUSS2D_VECTORIMAGE_H

#include <deque>
#include <vector>

#include "image.h"

namespace lsst::gauss2d {

#pragma GCC visibility push(hidden)
/*
    A basic Image class using a vector of deques.

    This very basic implementation is mainly for testing purposes.
    It almost certainly does not perform well and should not be used
    in production.
*/
template <typename t>
class VectorImage : public gauss2d::Image<t, VectorImage<t>> {
private:
    const size_t _n_rows;
    const size_t _n_cols;

    // This is a workaround for the C++98 specialization of vector<bool>
    std::vector<std::deque<t>> _data;

public:
    inline t& _get_value_unchecked(size_t row, size_t col) { return this->_data[row][col]; };
    const inline t get_value_unchecked(size_t row, size_t col) const { return this->_data[row][col]; };

    size_t get_n_cols() const { return _n_cols; };
    size_t get_n_rows() const { return _n_rows; };

    VectorImage(size_t n_rows, size_t n_cols,
                const std::shared_ptr<const gauss2d::CoordinateSystem> coordsys = nullptr)
            : gauss2d::Image<t, VectorImage<t>>(coordsys), _n_rows(n_rows), _n_cols(n_cols) {
        _data.resize(n_rows);
        for (size_t row = 0; row < n_rows; row++) {
            _data[row].resize(n_cols);
        }
    }
    ~VectorImage() = default;
};
#pragma GCC visibility pop

}  // namespace lsst::gauss2d
#endif