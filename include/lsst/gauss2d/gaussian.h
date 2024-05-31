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

#ifndef LSST_GAUSS2D_GAUSSIAN_H
#define LSST_GAUSS2D_GAUSSIAN_H

#include <algorithm>
#include <iostream>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "centroid.h"
#include "ellipse.h"
#include "object.h"

namespace lsst::gauss2d {

/**
 * @brief Interface for the normalization (total integrated value) of a 2D Gaussian.
 *
 */
class GaussianIntegral : public Object {
public:
    virtual double get_value() const = 0;
    virtual void set_value(double value) = 0;

    virtual std::string repr(bool name_keywords = false) const override = 0;
    virtual std::string str() const override = 0;

    virtual bool operator==(const GaussianIntegral& other) const {
        return this->get_value() == other.get_value();
    }

    virtual ~GaussianIntegral() = default;
};

/**
 * @brief A GaussianIntegral storing a float value.
 *
 * @note At the moments, limits are not enforced as there is no way to prevent
 * external methods setting the value negative and/or non-finite.
 *
 **/
class GaussianIntegralValue : public GaussianIntegral {
private:
    // TODO: Add some value safety to this
    // Probably must be a thin wrapper with a getter/setter enforcing >= 0
    // either that or delete the shared_ptr constructor and add a copy constructor
    std::shared_ptr<double> _value;

public:
    double get_value() const override { return *_value; }
    void set_value(double value) override { *_value = value; }

    std::string repr(bool name_keywords = false) const override;
    std::string str() const override;

    GaussianIntegralValue(double value = 1.);
    GaussianIntegralValue(std::shared_ptr<double> value);

    ~GaussianIntegralValue(){};
};

/**
 * @brief A 2D Gaussian with a Centroid, Ellipse, and integral.
 *
 * Gaussian offers some convenience functions but is otherwise a container
 * for its three component subclasses.
 **/
class Gaussian : public Object {
private:
    std::shared_ptr<Centroid> _centroid;
    std::shared_ptr<Ellipse> _ellipse;
    std::shared_ptr<GaussianIntegral> _integral;

    /**
     * @brief Check if a pointer is null and throw if so
     *
     * @tparam T The type of the object to check
     * @param ptr The pointer to check
     * @param name The name of the pointer to include in the error message (if null)
     * @return std::shared_ptr<T> The pointer that passed the check
     */
    template <typename T>
    std::shared_ptr<T> _check_not_nullptr(std::shared_ptr<T> ptr, std::string name) {
        if (ptr == nullptr) throw std::invalid_argument(this->str() + "Can't set " + name + " to nullptr");
        return ptr;
    }

public:
    /// Get the multiplicative factor for Gaussian function evaluations: integral/(2*area)
    double get_const_normal() const;
    /// Get the integral value
    double get_integral_value() const;

    /// Get the centroid object
    Centroid& get_centroid();
    /// Get the ellipse object
    Ellipse& get_ellipse();
    /// Get the integral object
    GaussianIntegral& get_integral();

    std::shared_ptr<Centroid> get_centroid_ptr();
    std::shared_ptr<Ellipse> get_ellipse_ptr();
    std::shared_ptr<GaussianIntegral> get_integral_ptr();

    const Centroid& get_centroid_const() const;
    const Ellipse& get_ellipse_const() const;
    const GaussianIntegral& get_integral_const() const;

    void set_const_normal(double const_normal);
    void set_integral_value(double integral);

    void set_centroid_ptr(std::shared_ptr<Centroid> centroid);
    void set_ellipse_ptr(std::shared_ptr<Ellipse> ellipse);
    void set_integral_ptr(std::shared_ptr<GaussianIntegral> integral);

    std::string repr(bool name_keywords = false) const override;
    std::string str() const override;

    bool operator==(const Gaussian& other) const;
    bool operator!=(const Gaussian& other) const;

    /**
     * @brief Construct a new Gaussian object
     *
     * @param centroid The centroid. Defaults to a new, default Centroid.
     * @param ellipse The ellipse. Defaults to a new, default Ellipse.
     * @param integral The integral. Defaults to a new, default GaussianIntegralValue.
     */
    Gaussian(std::shared_ptr<Centroid> centroid = nullptr, std::shared_ptr<Ellipse> ellipse = nullptr,
             std::shared_ptr<GaussianIntegral> integral = nullptr);
    ~Gaussian();
};

/**
 * @brief An array of Gaussian objects.
 *
 * This class exists partly to be an immutable container of Gaussians with
 * convenient constructors, but also so that it can be neatly wrapped with
 * pybind11.
 *
 */
class Gaussians : public Object {
public:
    typedef std::vector<std::shared_ptr<Gaussian>> Data;

private:
    Data _data = {};

    size_t assign(const Data& data, size_t i = 0);

public:
    Gaussian& operator[](size_t i);
    const Gaussian& operator[](size_t i) const;

    Gaussian& at(size_t i) const;
    const Gaussian& at_const(size_t i) const;

    using iterator = typename Data::iterator;
    using const_iterator = typename Data::const_iterator;

    typename Data::iterator begin() noexcept;
    typename Data::const_iterator cbegin() const noexcept;

    typename Data::const_iterator begin() const noexcept;
    typename Data::const_iterator end() const noexcept;

    typename Data::iterator end() noexcept;
    typename Data::const_iterator cend() const noexcept;

    Data get_data() const;

    size_t size() const;

    std::string repr(bool name_keywords = false) const override;
    std::string str() const override;

    // These constructors explicitly copy inputs rather than moving
    Gaussians(std::optional<const Data> data);
    Gaussians(std::vector<std::optional<const Data>> data);
};

/**
 * A convolution of a Gaussian source and kernel.
 *
 **/
class ConvolvedGaussian : public Object {
private:
    std::shared_ptr<const Gaussian> _source;
    std::shared_ptr<const Gaussian> _kernel;

public:
    const Gaussian& get_source() const;
    const Gaussian& get_kernel() const;

    std::unique_ptr<Gaussian> make_convolution() const;

    std::string repr(bool name_keywords = false) const override;
    std::string str() const override;

    ConvolvedGaussian(std::shared_ptr<const Gaussian> source = nullptr,
                      std::shared_ptr<const Gaussian> kernel = nullptr);
};

/**
 * A collection of ConvolvedGaussian objects.
 *
 * This class exists largely for the same reason as Gaussians, and to be passed
 * to evaluators.
 *
 */
class ConvolvedGaussians : public Object {
public:
    typedef std::vector<std::shared_ptr<ConvolvedGaussian>> Data;

private:
    Data _data = {};

    size_t assign(const Data& data, size_t i = 0);

public:
    ConvolvedGaussian& at(size_t i) const;
    const ConvolvedGaussian& at_const(size_t i) const;

    using iterator = typename Data::iterator;
    using const_iterator = typename Data::const_iterator;

    typename Data::iterator begin() noexcept;
    typename Data::iterator end() noexcept;

    typename Data::const_iterator begin() const noexcept;
    typename Data::const_iterator end() const noexcept;

    typename Data::const_iterator cbegin() const noexcept;
    typename Data::const_iterator cend() const noexcept;

    Data get_data() const;

    size_t size() const;

    std::string repr(bool name_keywords = false) const override;
    std::string str() const override;

    ConvolvedGaussian& operator[](size_t i);
    const ConvolvedGaussian& operator[](size_t i) const;

    // These constructors explicitly copy inputs rather than moving
    ConvolvedGaussians(std::optional<const Data> data);
    ConvolvedGaussians(std::vector<std::optional<const Data>> data);
};

}  // namespace lsst::gauss2d
#endif
