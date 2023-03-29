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

#ifndef GAUSS2D_IMAGE_H
#define GAUSS2D_IMAGE_H

#include <array>
#include <cstddef>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "object.h"
#include "type_name.h"

namespace gauss2d {

/**
 * A coordinate system specifying image scale and orientation.
 * 
 * This is intended to mimic some of the functionality of e.g. basic
 * FITS header, in order for evaluators to draw images at different
 * scales and/or with rotation/translation.
 * 
**/
class CoordinateSystem : public Object
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
        return (
            _dx1 == other.get_dx1()) && (_dy2 == other.get_dy2())
            && (_x_min == other.get_x_min()) && (_y_min == other.get_y_min()
        );
    }

    std::string repr(bool name_keywords) const {
        return std::string("CoordinateSystem(")
            + (name_keywords ? "d1=" : "") + std::to_string(_dx1) + ", "
            + (name_keywords ? "d2=" : "") + std::to_string(_dy2) + ")";
    }

    std::string str() const override {
        return "CoordinateSystem(d1=" + std::to_string(_dx1)
            + ", d2=" + std::to_string(_dy2) + ")";
    }

    CoordinateSystem(double d1=1., double d2=1) : _dx1(d1), _dy2(d2) {};
};

static const CoordinateSystem COORDS_DEFAULT {}; 

template<typename t, class Data, class Indices>
class GaussianEvaluator;

template <typename t, class C>
class Image;

// TODO: Check if this can/should be a static method in Image
// Consider that t1/t2 might be reversed in two interchangeable functions.

/**
 * @brief Return if two images are compatible.
 * 
 * Compatible means that they have the same dimensions in both axes and
 * equivalent coordinate systems.
 * 
 * @tparam t1 
 * @tparam C1 
 * @tparam t2 
 * @tparam C2 
 * @param img1 
 * @param img2 
 * @param msg 
 * @return true If images are compatible. 
**/
template<typename t1, class C1, typename t2, class C2>
bool images_compatible(const Image<t1, C1> & img1, const Image<t2, C2> & img2,
    std::string * msg = nullptr)
{
    bool coordsys_equal = img1.get_coordsys() == img2.get_coordsys();
    bool return_msg = msg != nullptr;
    if(!return_msg && !coordsys_equal) return false;
    bool cols_equal = img1.get_n_cols() == img2.get_n_cols();
    bool rows_equal = img1.get_n_rows() == img2.get_n_rows();
    bool passed = coordsys_equal && cols_equal && rows_equal;
    if(!passed) {
        if(return_msg) {
            if(!coordsys_equal) *msg += img1.get_coordsys().str() + "!=" + img2.get_coordsys().str() + ",";
            if(!cols_equal) *msg += std::to_string(img1.get_n_cols()) + "!=" + 
                std::to_string(img2.get_n_cols()) + ",";
            if(!rows_equal) *msg += std::to_string(img1.get_n_rows()) + "!=" + 
                std::to_string(img2.get_n_rows());
        }
        return false;
    }
    return true;
}

/**
 * A 2D image with scalar numeric values, using CRTP.
 * 
 * Basic implementations of most functions are provided. Derived classes
 * should override any and all if the default implementations are not
 * effecient enough.
 * 
 * @tparam t The numeric type.
 * @tparam C The specialized class.
 * 
**/
template <typename t, class C>
class Image : public Object
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
    virtual inline t & _get_value_unchecked(size_t row, size_t col) = 0;

    void _check_row_col(size_t row, size_t col) const {
        if(!(row < this->get_n_rows()) && !(row < this->get_n_rows())) {
            throw std::out_of_range("row,col = " + std::to_string(row) + ","  + std::to_string(col)
                + " n_rows,n_cols = " + std::to_string(this->get_n_rows()) + ","
                + std::to_string(this->get_n_cols())
            );
        }
    }

    virtual const CoordinateSystem & get_coordsys() const { return _coordsys; };
    virtual std::shared_ptr<const CoordinateSystem> get_coordsys_ptr_const() const { return _coordsys_ptr; };

    virtual size_t get_n_cols() const = 0;
    virtual size_t get_n_rows() const = 0;

    void add_value(size_t row, size_t col, t value) { this->_get_value(row, col) += value; }
    void add_value_unchecked(size_t row, size_t col, t value) {
        static_cast<C&>(*this)._get_value_unchecked(row, col) += value;
    }
    virtual void fill(t value) {
        const size_t n_rows = get_n_rows();
        const size_t n_cols = get_n_cols();
        for(size_t row = 0; row < n_rows; ++row) {
            for(size_t col = 0; col < n_cols; ++col) {
                this->set_value_unchecked(row, col, value);
            }
        }
    }
    inline t get_value(size_t row, size_t col) const {
        _check_row_col(row, col);
        return static_cast<const C&>(*this).get_value_unchecked(row, col);
    }
    virtual const inline t get_value_unchecked(size_t row, size_t col) const = 0;
    inline void set_value(size_t row, size_t col, t value) { _get_value(row, col) = value; }
    inline void set_value_unchecked(size_t row, size_t col, t value) {
        static_cast<C&>(*this)._get_value_unchecked(row, col) = value;
    }

    std::array<size_t, 2> shape() const { return {this->get_n_rows(), this->get_n_cols()}; }

    size_t size() const { return this->get_n_rows() * this->get_n_cols();};

    std::string repr(bool name_keywords) const {
        return std::string(type_name<C>()) + "("
            + (name_keywords ? "coordsys=" : "") +  _coordsys.repr(name_keywords) + ", "
            + (name_keywords ? "n_rows=" : "") + std::to_string(this->get_n_rows()) + ", "
            + (name_keywords ? "n_cols=" : "") + std::to_string(this->get_n_cols()) + ")";
    }

    std::string str() const override {
        return std::string(type_name<C>()) + "(coordsys=" + _coordsys.str()
            + ", n_rows=" + std::to_string(this->get_n_rows())
            + ", n_cols=" + std::to_string(this->get_n_cols())
            + ")";
    }

    virtual void operator += (t value) {
        const size_t n_rows = get_n_rows();
        const size_t n_cols = get_n_cols();
        for(size_t row = 0; row < n_rows; ++row) {
            for(size_t col = 0; col < n_cols; ++col) {
                this->add_value_unchecked(row, col, value);
            }
        }
    }

    bool operator == (const Image &other) const {
        if(images_compatible<t, C, t, C>(*this, other)) {
            const size_t n_rows = get_n_rows();
            const size_t n_cols = get_n_cols();
            for(size_t row = 0; row < n_rows; ++row) {
                for(size_t col = 0; col < n_cols; ++col) {
                    if(this->get_value_unchecked(row, col) != other.get_value_unchecked(row, col)) {
                        return false;
                    }
                }
            }
            return true;
        }
        return false;
    }

    const bool operator != (const Image &other) const {
        return !(*this == other);
    }

    // TODO: Figure out if there's any point to this in CRTP (or otherwise)
    Image(
        size_t n_rows,
        size_t n_cols,
        const std::shared_ptr<const CoordinateSystem> coordsys = nullptr
    ) = delete;

    Image(const std::shared_ptr<const CoordinateSystem> coordsys = nullptr) : 
        _coordsys_ptr(coordsys == nullptr ? nullptr : std::move(coordsys)),
        _coordsys(coordsys == nullptr ? COORDS_DEFAULT : *_coordsys_ptr)
    {}
    virtual ~Image() = default;
};

/**
 * @brief An array of compatible Images.
 *
 * See images_compatible() for the definition of compatibility.
 *
 * @tparam t The numeric type.
 * @tparam C The specialized class.
 */
template <typename t, class C>
class ImageArray : public Object
{
public:
    typedef Image<t, C> ImageT;
    typedef std::vector<std::shared_ptr<C>> Data;

private:
    Data _images {};

public:
    ImageT& operator [] (size_t i) {return *(_images[i]);}
    const ImageT& operator [] (size_t i) const {return *(_images[i]);}

    ImageT & at(size_t i=0) const {return (*_images.at(i));}

    using iterator = typename Data::iterator;
    using const_iterator = typename Data::const_iterator;

    typename Data::iterator begin() noexcept { return _images.begin(); }
    typename Data::iterator end() noexcept  { return _images.end(); };

    typename Data::const_iterator cbegin() const noexcept  { return _images.cbegin(); };
    typename Data::const_iterator cend() const noexcept  { return _images.begin(); };


    size_t size() const {return _images.size();}

    std::string repr(bool name_keywords) const {
        std::string str = std::string(type_name<C>()) + "(" + (name_keywords ? "data=[" : "[");
        for(auto img = this->cbegin(); img != this->cend(); ++img) str += (*img)->repr(name_keywords) + ",";
        return str + "])";
    }

    std::string str() const override {
        std::string str = std::string(type_name<C>()) + "(data=[";
        for(auto img = this->cbegin(); img != this->cend(); ++img) str += (*img)->str() + ",";
        return str + "])";
    }

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
#endif // GAUSS2D_IMAGE_H
