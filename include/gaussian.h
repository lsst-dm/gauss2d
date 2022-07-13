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

#ifndef GAUSS2D_GAUSSIAN_H
#define GAUSS2D_GAUSSIAN_H

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

namespace gauss2d {

/**
 * Interface for the normalization (total integrated value) of a 2D Gaussian.
 * 
**/
class GaussianIntegral : public Object
{
public:
    virtual double get_value() const = 0;
    virtual void set_value(double value) = 0;

    virtual std::string str() const override = 0;

    virtual bool operator==(const GaussianIntegral& other) const {
        return this->get_value() == other.get_value();
    }

    virtual ~GaussianIntegral() = default;
};

/**
 * A GaussianIntegral storing a float value.
 * 
**/
class GaussianIntegralValue : public GaussianIntegral
{
private:
    double _value;

public:
    double get_value() const override { return _value; }
    void set_value(double value) override { _value = value; }

    std::string str() const override { return "GaussianIntegralValue(" + std::to_string(_value) + ")"; }

    GaussianIntegralValue(double value=1.): _value(value) {};
    ~GaussianIntegralValue() {};
};

/**
 * A 2D Gaussian with a Centroid, Ellipse, and integral.
 * 
 * Gaussian offers some convenience functions but is otherwise a container
 * for its three component subclasses.
**/
class Gaussian : public Object
{
private:
    std::shared_ptr<Centroid> _centroid;
    std::shared_ptr<Ellipse> _ellipse;
    std::shared_ptr<GaussianIntegral> _integral;

    template <typename T>
    std::shared_ptr<T> _check_not_nullptr(std::shared_ptr<T> ptr, std::string name) {
        if(ptr == nullptr) throw std::invalid_argument(this->str() + "Can't set " + name + " to nullptr");
        return ptr;
    }

public:
    double get_const_normal() const;
    double get_integral_value() const;

    Centroid & get_centroid();
    Ellipse & get_ellipse();
    GaussianIntegral & get_integral();

    std::shared_ptr<Centroid> get_centroid_ptr();
    std::shared_ptr<Ellipse> get_ellipse_ptr();
    std::shared_ptr<GaussianIntegral> get_integral_ptr();

    const Centroid & get_centroid_const() const;
    const Ellipse & get_ellipse_const() const;
    const GaussianIntegral & get_integral_const() const;

    void set_const_normal(double const_normal);
    void set_integral_value(double integral);

    void set_centroid_ptr(std::shared_ptr<Centroid> centroid);
    void set_ellipse_ptr(std::shared_ptr<Ellipse> ellipse);
    void set_integral_ptr(std::shared_ptr<GaussianIntegral> integral);

    std::string str() const override;

    bool operator == (const Gaussian& other) const;
    bool operator != (const Gaussian& other) const;

    Gaussian(std::shared_ptr<Centroid> centroid = nullptr, std::shared_ptr<Ellipse> ellipse = nullptr,
             std::shared_ptr<GaussianIntegral> integral = nullptr);
    ~Gaussian();
};

/**
 * A collection of Gaussian objects.
 * 
 * This class exists partly to be an immutable container of Gaussians with 
 * convenient constructors, but also so that it can be neatly wrapped with
 * pybind11.
 * 
 */
class Gaussians : public Object
{
public:
    typedef std::vector<std::shared_ptr<Gaussian>> Data;

private:
    Data _data = {};

    size_t assign(const Data & data, size_t i = 0);

public:
    Gaussian & operator [] (size_t i);
    const Gaussian & operator [] (size_t i) const;

    Gaussian & at(size_t i) const;
    const Gaussian & at_const(size_t i) const;

    using iterator = typename Data::iterator;
    using const_iterator = typename Data::const_iterator;

    typename Data::iterator begin() noexcept;
    typename Data::const_iterator cbegin() const noexcept;

    typename Data::iterator end() noexcept;
    typename Data::const_iterator cend() const noexcept;

    Data get_data() const;

    size_t size() const;

    std::string str() const override;

    // These constructors explicitly copy inputs rather than moving
    Gaussians(std::optional<const Data> data);
    Gaussians(std::vector<std::optional<const Data>> data);
};

/**
 * A convolution of a Gaussian source and kernel.
 * 
**/
class ConvolvedGaussian : public Object
{
private:
    std::shared_ptr<Gaussian> _source;
    std::shared_ptr<Gaussian> _kernel;

public:
    Gaussian & get_source();
    Gaussian & get_kernel();

    const Gaussian & get_source_const() const;
    const Gaussian & get_kernel_const() const;

    std::unique_ptr<Gaussian> make_convolution() const;

    std::string str() const override;

    ConvolvedGaussian(
        std::shared_ptr<Gaussian> source = nullptr,
        std::shared_ptr<Gaussian> kernel = nullptr
    );
};

/**
 * A collection of ConvolvedGaussian objects.
 * 
 * This class exists largely for the same reason as Gaussians, and to be passed
 * to evaluators.
 * 
 */
class ConvolvedGaussians : public Object
{
public:
    typedef std::vector<std::shared_ptr<ConvolvedGaussian>> Data;

private:
    Data _data = {};

    size_t assign(const Data & data, size_t i = 0);

public:
    ConvolvedGaussian & at(size_t i) const;
    const ConvolvedGaussian & at_const(size_t i) const;
    
    using iterator = typename Data::iterator;
    using const_iterator = typename Data::const_iterator;

    typename Data::iterator begin() noexcept;
    typename Data::iterator end() noexcept;

    typename Data::const_iterator cbegin() const noexcept;
    typename Data::const_iterator cend() const noexcept;

    Data get_data() const;

    size_t size() const;

    std::string str() const override;

    ConvolvedGaussian & operator [] (size_t i);
    const ConvolvedGaussian & operator [] (size_t i) const;

    // These constructors explicitly copy inputs rather than moving
    ConvolvedGaussians(std::optional<const Data> data);
    ConvolvedGaussians(std::vector<std::optional<const Data>> data);
};

} // namespace gauss2d
#endif
