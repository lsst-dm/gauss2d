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

#ifndef GAUSS2D_CARRAYIMAGE_H
#define GAUSS2D_CARRAYIMAGE_H

#include "image.h"

namespace gauss2d {

/*
    An Image stored as C-style raw arrays.

    The bool implementation of this class will not be memory-efficient,
    but the others ought to be performant.
*/
template <typename t>
class CArrayImage : public gauss2d::Image<t, CArrayImage<t>> {
private:
    const size_t _n_rows;
    const size_t _n_cols;

    t ** _data;

public:
    inline t& _get_value_unchecked_impl(size_t row, size_t col) {
        return this->_data[row][col];
    };
    inline t get_value_unchecked_impl(size_t row, size_t col) const {
        return this->_data[row][col];
    };

    size_t get_n_cols_impl() const { return _n_cols; };
    size_t get_n_rows_impl() const { return _n_rows; };

    CArrayImage(size_t n_rows, size_t n_cols,
                const t* value_init = Image<t, CArrayImage<t>>::_value_default_ptr(),
                const std::shared_ptr<const CoordinateSystem> coordsys = nullptr)
            : Image<t, CArrayImage<t>>(coordsys), _n_rows(n_rows), _n_cols(n_cols) {
        _data = new t*[_n_rows];
        bool do_init = value_init != nullptr;
        t value_init_override = do_init ? *value_init : 0;
        for (size_t row = 0; row < _n_rows; row++) {
            _data[row] = new t[_n_cols];
            if(do_init) {
                for (size_t col = 0; col < _n_cols; col++) {
                    _data[row][col] = value_init_override;
                }
            }
        }
    }
    ~CArrayImage() {
        for (size_t row = 0; row < _n_rows; row++) {
            delete [] _data[row];
        }
        delete [] _data;
    }
};

}  // namespace gauss2d
#endif