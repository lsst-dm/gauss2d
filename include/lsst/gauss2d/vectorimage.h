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

/*
    A basic Image class using a vector of deques.

    This very basic implementation is mainly for testing purposes.
    It almost certainly does not perform well and should not be used
    in production.
*/
template <typename T>
class VectorImage : public gauss2d::Image<T, VectorImage<T>> {
public:
    explicit VectorImage(size_t n_rows, size_t n_cols,
                         const T* value_init = Image<T, VectorImage<T>>::_value_default_ptr(),
                         std::shared_ptr<const CoordinateSystem> coordsys = nullptr)
            : Image<T, VectorImage<T>>(coordsys), _n_rows(n_rows), _n_cols(n_cols) {
        _data.resize(n_rows);
        // No real option but to resize with vectors
        // Just reserving would cause _unchecked calls to break
        T value_init_override
                = value_init != nullptr ? *value_init : (Image<T, VectorImage<T>>::_value_default);
        for (size_t row = 0; row < _n_rows; row++) {
            _data[row].resize(n_cols, value_init_override);
        }
    }
    ~VectorImage() = default;

    T& _get_value_impl(size_t row, size_t col) {
        // This doesn't work on vector<bool> because it's bit-packed
        // One could specialize: if constexpr (std::is_same_v<bool, T>)
        // ... if there were a workable alternative, but there isn't.
        if constexpr (std::is_same_v<bool, T>) {
            throw std::logic_error("VectorImage<bool> cannot use _get_value");
        } else {
            // TODO: Check at performance vs default implementation
            return this->_data.at(row).at(col);
        }
    }
    inline T& _get_value_unchecked_impl(size_t row, size_t col) {
        // See note in _get_value
        if constexpr (std::is_same_v<bool, T>) {
            throw std::logic_error("VectorImage<bool> cannot use _get_value_unchecked");
        } else {
            return this->_data[row][col];
        }
    }
    void add_value_unchecked_impl(size_t row, size_t col, T value) { this->_data[row][col] += value; }
    inline T get_value_unchecked_impl(size_t row, size_t col) const { return this->_data[row][col]; }
    inline void set_value_impl(size_t row, size_t col, T value) { this->_data.at(row).at(col) = value; }
    inline void set_value_unchecked_impl(size_t row, size_t col, T value) { this->_data[row][col] = value; }

    size_t get_n_cols_impl() const { return _n_cols; };
    size_t get_n_rows_impl() const { return _n_rows; };

private:
    const size_t _n_rows;
    const size_t _n_cols;

    // This is a workaround for the C++98 specialization of vector<bool>
    std::vector<std::deque<T>> _data;
};

}  // namespace lsst::gauss2d
#endif