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

#ifndef LSST_GAUSS2D_IMAGE_H
#define LSST_GAUSS2D_IMAGE_H

#include <array>
#include <cstddef>
#include <iostream>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "coordinatesystem.h"
#include "object.h"
#include "type_name.h"

namespace lsst::gauss2d {

typedef size_t idx_type;

template <typename T, class Data, class Indices>
class GaussianEvaluator;

template <typename T, class C>
class Image;

// TODO: Check if this can/should be a static method in Image
// Consider that T1/T2 might be reversed in two interchangeable functions.

/**
 * @brief Return if two images are compatible.
 *
 * Compatible means that they have the same dimensions in both axes and
 * equivalent coordinate systems.
 *
 * @tparam T1 The data type of C1.
 * @tparam C1 The class of the first Image.
 * @tparam T2 The data type of C2.
 * @tparam C2 The class of the second Image.
 * @param img1 The first image.
 * @param img2 The second image.
 * @param msg A string to append error messages to, if not null.
 * @return true If images are compatible.
 **/
template <typename T1, class C1, typename T2, class C2>
bool images_compatible(const Image<T1, C1>& img1, const Image<T2, C2>& img2, bool compare_coordsys = true,
                       std::string* msg = nullptr) {
    bool coordsys_equal = !compare_coordsys || (img1.get_coordsys() == img2.get_coordsys());
    bool return_msg = msg != nullptr;
    if (!return_msg && !coordsys_equal) return false;
    bool cols_equal = img1.get_n_cols() == img2.get_n_cols();
    bool rows_equal = img1.get_n_rows() == img2.get_n_rows();
    bool passed = coordsys_equal && cols_equal && rows_equal;
    if (!passed) {
        if (return_msg) {
            if (!coordsys_equal) {
                *msg += img1.get_coordsys().str() + "!=" + img2.get_coordsys().str() + ",";
            }
            if (!cols_equal) {
                *msg += std::to_string(img1.get_n_cols()) + "!=" + std::to_string(img2.get_n_cols()) + ",";
            }
            if (!rows_equal) {
                *msg += std::to_string(img1.get_n_rows()) + "!=" + std::to_string(img2.get_n_rows());
            }
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
 * efficient enough.
 *
 * @tparam t The numeric type.
 * @tparam C The specialized class.
 *
 **/
template <typename T, class C>
class Image : public Object {
public:
    // TODO: Figure out if there's any point to this in CRTP (or otherwise)
    explicit Image(size_t n_rows, size_t n_cols, const T* value_init = _value_default_ptr(),
                   std::shared_ptr<const CoordinateSystem> coordsys = nullptr)
            = delete;

    // Convenience initializer for a coordsys (could be private method?)
    explicit Image(std::shared_ptr<const CoordinateSystem> coordsys = nullptr)
            : _coordsys_ptr(coordsys == nullptr ? nullptr : std::move(coordsys)),
              _coordsys(_coordsys_ptr == nullptr ? COORDS_DEFAULT : *_coordsys_ptr) {}
    ~Image() = default;

    static constexpr T _value_default = 0;
    static const T* _value_default_ptr() { return &_value_default; };

    T& _get_value(size_t row, size_t col) { return static_cast<C&>(*this)._get_value_impl(row, col); }
    T& _get_value_impl(size_t row, size_t col) {
        _check_row_col(row, col);
        return this->_get_value_unchecked(row, col);
    }
    inline T& _get_value_unchecked(size_t row, size_t col) {
        return self()._get_value_unchecked_impl(row, col);
    }
    inline T& _get_value_unchecked_impl(size_t row, size_t col) = delete;

    void _check_row_col(size_t row, size_t col) const { return self_const()._check_row_col_impl(row, col); }
    void _check_row_col_impl(size_t row, size_t col) const {
        if (!((row < this->get_n_rows()) && (col < this->get_n_cols()))) {
            throw std::out_of_range("row,col = " + std::to_string(row) + "," + std::to_string(col)
                                    + " n_rows,n_cols = " + std::to_string(this->get_n_rows()) + ","
                                    + std::to_string(this->get_n_cols()));
        }
    }

    const CoordinateSystem& get_coordsys() const { return _coordsys; };
    std::shared_ptr<const CoordinateSystem> get_coordsys_ptr_const() const { return _coordsys_ptr; };

    size_t get_n_cols() const { return static_cast<const C&>(*this).get_n_cols_impl(); }
    size_t get_n_cols_impl() const = delete;
    size_t get_n_rows() const { return static_cast<const C&>(*this).get_n_rows_impl(); }
    size_t get_n_rows_impl() = delete;

    void add_value(size_t row, size_t col, T value) { self().add_value_impl(row, col, value); }
    void add_value_impl(size_t row, size_t col, T value) { this->_get_value(row, col) += value; }
    void add_value_unchecked(size_t row, size_t col, T value) {
        self().add_value_unchecked_impl(row, col, value);
    }
    void add_value_unchecked_impl(size_t row, size_t col, T value) {
        self()._get_value_unchecked(row, col) += value;
    }
    void fill(T value) { self().fill_impl(value); }
    void fill_impl(T value) {
        const size_t n_rows = get_n_rows();
        const size_t n_cols = get_n_cols();
        for (size_t row = 0; row < n_rows; ++row) {
            for (size_t col = 0; col < n_cols; ++col) {
                this->set_value_unchecked(row, col, value);
            }
        }
    }
    inline T get_value(size_t row, size_t col) const { return self_const().get_value_impl(row, col); }
    inline T get_value_impl(size_t row, size_t col) const {
        _check_row_col(row, col);
        return self_const().get_value_unchecked(row, col);
    }
    inline T get_value_unchecked(size_t row, size_t col) const {
        return self_const().get_value_unchecked_impl(row, col);
    }
    inline T get_value_unchecked_impl(size_t row, size_t col) const = delete;
    inline void set_value(size_t row, size_t col, T value) { return self().set_value_impl(row, col, value); }
    inline void set_value_impl(size_t row, size_t col, T value) { self()._get_value(row, col) = value; }
    inline void set_value_unchecked(size_t row, size_t col, T value) {
        self().set_value_unchecked_impl(row, col, value);
    }
    inline void set_value_unchecked_impl(size_t row, size_t col, T value) {
        self()._get_value_unchecked(row, col) = value;
    }

    std::array<size_t, 2> shape() const { return {this->get_n_rows(), this->get_n_cols()}; }

    size_t size() const { return this->get_n_rows() * this->get_n_cols(); };

    std::string repr(bool name_keywords, std::string_view namespace_separator) const override {
        return type_name_str<C>(false, namespace_separator) + "(" + (name_keywords ? "coordsys=" : "")
               + _coordsys.repr(name_keywords, namespace_separator) + ", " + (name_keywords ? "n_rows=" : "")
               + std::to_string(this->get_n_rows()) + ", " + (name_keywords ? "n_cols=" : "")
               + std::to_string(this->get_n_cols()) + ")";
    }

    std::string str() const override {
        return type_name_str<C>(true) + "(coordsys=" + _coordsys.str() + ", n_rows="
               + std::to_string(this->get_n_rows()) + ", n_cols=" + std::to_string(this->get_n_cols()) + ")";
    }

    Image<T, C>& operator+=(T value) {
        const size_t n_rows = get_n_rows();
        const size_t n_cols = get_n_cols();
        for (size_t row = 0; row < n_rows; ++row) {
            for (size_t col = 0; col < n_cols; ++col) {
                this->add_value_unchecked(row, col, value);
            }
        }
        return *this;
    }

    Image<T, C>& operator*=(T value) {
        const size_t n_rows = get_n_rows();
        const size_t n_cols = get_n_cols();
        for (size_t row = 0; row < n_rows; ++row) {
            for (size_t col = 0; col < n_cols; ++col) {
                // Avoid annoying warning for bool case
                if constexpr (std::is_same_v<T, bool>) {
                    this->set_value_unchecked(row, col, value && this->get_value_unchecked(row, col));
                } else {
                    this->set_value_unchecked(row, col, value * this->get_value_unchecked(row, col));
                }
            }
        }
        return *this;
    }

    // TODO: Implement if deemed worthwhile
    // void operator+=(T value);

    bool operator==(const Image& other) const {
        if (images_compatible<T, C, T, C>(*this, other)) {
            const size_t n_rows = get_n_rows();
            const size_t n_cols = get_n_cols();
            for (size_t row = 0; row < n_rows; ++row) {
                for (size_t col = 0; col < n_cols; ++col) {
                    if (this->get_value_unchecked(row, col) != other.get_value_unchecked(row, col)) {
                        return false;
                    }
                }
            }
            return true;
        }
        return false;
    }

    const bool operator!=(const Image& other) const { return !(*this == other); }

private:
    friend GaussianEvaluator<T, class Data, class Indices>;
    const std::shared_ptr<const CoordinateSystem> _coordsys_ptr;
    const gauss2d::CoordinateSystem& _coordsys;

    inline C& self() { return static_cast<C&>(*this); };
    inline const C& self_const() const { return static_cast<const C&>(*this); };
};

/**
 * @brief An array of compatible Images.
 *
 * See images_compatible() for the definition of compatibility.
 *
 * @tparam T The numeric type.
 * @tparam C The specialized class.
 */
template <typename T, class C>
class ImageArray : public Object {
public:
    typedef Image<T, C> ImageT;
    typedef std::vector<std::shared_ptr<C>> Data;

    ImageT& operator[](size_t i) { return *(_images[i]); }
    const ImageT& operator[](size_t i) const { return *(_images[i]); }

    ImageT& at(size_t i = 0) const { return (*_images.at(i)); }

    using iterator = typename Data::iterator;
    using const_iterator = typename Data::const_iterator;

    typename Data::iterator begin() noexcept { return _images.begin(); }
    typename Data::iterator end() noexcept { return _images.end(); };

    typename Data::const_iterator cbegin() const noexcept { return _images.cbegin(); };
    typename Data::const_iterator cend() const noexcept { return _images.begin(); };

    size_t size() const { return _images.size(); }

    std::string repr(bool name_keywords, std::string_view namespace_separator) const override {
        std::string str = type_name_str<ImageArray<T, C>>(false, namespace_separator) + "("
                          + (name_keywords ? "data=" : "")
                          + repr_iter_ptr(_images, name_keywords, namespace_separator) + ")";
        return str;
    }

    std::string str() const override {
        std::string str = type_name_str<ImageArray<T, C>>(true) + "(data=" + str_iter_ptr(_images) + ")";
        return str;
    }

    explicit ImageArray(const Data* data_in) {
        if (data_in != nullptr) {
            const Data& data = *data_in;
            size_t n_data = data.size();
            if (n_data > 0) {
                _images.resize(n_data);
                for (size_t i = 0; i < n_data; ++i) {
                    if (data[i] == nullptr) {
                        throw std::invalid_argument("ImageArray data[" + std::to_string(i)
                                                    + "] can't be null");
                    }
                    if (!images_compatible<T, C, T, C>(*data[i], *data[0])) {
                        throw std::invalid_argument("ImageArray data[" + std::to_string(i)
                                                    + "] must be compatible with data[0] (and all others)");
                    }
                    _images[i] = data[i];
                }
            }
        }
    };

private:
    Data _images{};
};

}  // namespace lsst::gauss2d
#endif  // LSST_GAUSS2D_IMAGE_H
