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

#ifndef __GAUSS2D_IMAGE_H_
#define __GAUSS2D_IMAGE_H_

#include <cstddef>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace gauss2d {

class CoordinateSystem
{
private:
    double _dx1;
    double _dy2;

    double _x_min = 0.;
    double _y_min = 0.;

public:
    // x,y coords of the bottom left corner of the image if not rotated
    // TBD how to implement this if there is rotation...

    double get_dx1() const { return _dx1;}
    double get_dy2() const { return _dy2;}

    double get_x_min() const { return _x_min;}
    double get_y_min() const { return _y_min;}

    bool is_xy_aligned() const { return true;};

    bool operator==(const CoordinateSystem & other) const {
        return (_dx1 == other.get_dx1()) && (_dy2 == other.get_dy2()) && (_x_min == other.get_x_min()) && (_y_min == other.get_y_min());
    }

    CoordinateSystem(double d1=1., double d2=1) : _dx1(d1), _dy2(d2) {};
};

static const CoordinateSystem COORDS_DEFAULT {}; 

template<typename t, class Data, class Indices>
class GaussianEvaluator;

template <typename t, class C>
class Image
{
private:
    friend GaussianEvaluator<t, class Data, class Indices>;
    const std::shared_ptr<const CoordinateSystem> _coordsys_ptr;
    const gauss2d::CoordinateSystem & _coordsys;

public:
    t & _get_value(size_t row, size_t col) {
        _check_row_col(row, col);
        return this->_get_value_unchecked(row, col);
    }
    inline t & _get_value_unchecked(size_t row, size_t col) { return static_cast<C&>(*this)._get_value_unchecked(row,  col); }

    void _check_row_col(size_t row, size_t col) const {
        if(!(row < this->get_n_rows()) && !(row < this->get_n_rows())) {
            throw std::out_of_range("row,col = " + std::to_string(row) + ","  + std::to_string(col)
                + "n_rows,n_cols = " + std::to_string(this->get_n_rows()) + ","
                + std::to_string(this->get_n_rows())
            );
        }
    }

    virtual const CoordinateSystem & get_coordsys() const { return _coordsys;};

    size_t get_n_cols() const { return static_cast<const C&>(*this).get_n_cols(); }
    size_t get_n_rows() const { return static_cast<const C&>(*this).get_n_rows(); }

    void add_value(size_t row, size_t col, t value) { this->_get_value(row, col) += value;}
    void add_value_unchecked(size_t row, size_t col, t value) {
        static_cast<C&>(*this)._get_value_unchecked(row, col) += value;
    }
    inline t get_value(size_t row, size_t col) const {
        _check_row_col(row, col);
        return static_cast<const C&>(*this).get_value_unchecked(row,  col);;
    };
    inline t get_value_unchecked(size_t row, size_t col) const { return static_cast<const C&>(*this).get_value_unchecked(row,  col);};
    inline void set_value(size_t row, size_t col, t value) { _get_value(row, col) = value;}
    inline void set_value_unchecked(size_t row, size_t col, t value) { static_cast<C&>(*this)._get_value_unchecked(row, col) = value;};

    size_t size() const { return this->get_n_rows() * this->get_n_cols();};

    Image(const std::shared_ptr<const CoordinateSystem> coordsys = nullptr) : 
        _coordsys_ptr(coordsys == nullptr ? nullptr : std::move(coordsys)),
        _coordsys(coordsys == nullptr ? COORDS_DEFAULT : *_coordsys_ptr)
    {}
    virtual ~Image() = default;
};

template<typename t1, class C1, typename t2, class C2>
bool images_compatible(const Image<t1, C1> & img1, const Image<t2, C2> & img2)
{
    return (img1.get_coordsys() == img2.get_coordsys()) && (img1.get_n_cols() == img2.get_n_cols()) &&
        (img1.get_n_rows() == img2.get_n_rows());
}

template <typename t, class C>
class ImageArray
{
public:
    typedef Image<t, C> ImageT;
    typedef std::vector<std::shared_ptr<C>> Data;

private:
    Data _images {};

public:
    ImageT& operator[](size_t i) {return *(_images[i]);}
    const ImageT& operator[](size_t i) const {return *(_images[i]);}

    typename Data::iterator begin() noexcept {return _images.begin();}
    typename Data::const_iterator cbegin() const noexcept {return _images.begin();}

    typename Data::iterator end() noexcept {return _images.end();}
    typename Data::const_iterator cend() const noexcept {return _images.cend();}

    ImageT & at(size_t i=0) const {return (*_images.at(i));}
    size_t size() const {return _images.size();}

    ImageArray(const Data * data_in)
    {
        if(data_in != nullptr)
        {
            const Data & data = *data_in;
            size_t n_data = data.size();
            if(n_data > 0)
            {
                _images.resize(n_data);
                for(size_t i = 0; i < n_data; ++i)
                {
                    if(data[i] == nullptr) throw std::runtime_error("ImageArray data[" + std::to_string(i) + "] can't be null");
                    if(!images_compatible<t, C, t, C>(*data[i], *data[0])) throw std::runtime_error("ImageArray data[" + 
                        std::to_string(i) + "] must be compatible with data[0] (and all others)");
                    _images[i] = data[i];
                }
            }
        }
    };
};

}
#endif