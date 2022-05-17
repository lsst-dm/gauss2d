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
#ifndef GAUSS2D_CENTROID_H
#include "centroid.h"
#endif

#ifndef GAUSS2D_ELLIPSE_H
#include "ellipse.h"
#endif

#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace gauss2d {

class GaussianIntegral
{
public:
    virtual double get_value() const = 0;
    virtual void set_value(double value) = 0;

    virtual std::string str() const = 0;
    virtual ~GaussianIntegral() = default;
};

class GaussianIntegralValue : public GaussianIntegral
{
private:
    double _value;

public:
    double get_value() const { return _value; }
    void set_value(double value) { _value = value; }

    std::string str() const { return "GaussianIntegralValue(" + std::to_string(_value) + ")"; }

    GaussianIntegralValue(double value=1.): _value(value) {};
    ~GaussianIntegralValue() {};
};

class Gaussian
{
private:
    std::shared_ptr<Centroid> _centroid;
    std::shared_ptr<Ellipse> _ellipse;
    std::shared_ptr<GaussianIntegral> _integral;

    template <typename T>
    void _set_if_not_nullptr(std::shared_ptr<T> x, std::shared_ptr<T> to_set, std::string name) {
        if(to_set == nullptr) throw std::invalid_argument(this->str() + "Can't set " + name + " to nullptr");
        x = std::move(to_set);
    }

public:
    double get_const_normal() const { return _integral->get_value()/(2*_ellipse->get_area()); }
    double get_integral_value() const {return _integral->get_value();};

    Centroid & get_centroid() { return *_centroid; }
    Ellipse & get_ellipse() { return *_ellipse; }
    GaussianIntegral & get_integral() { return *_integral; }

    std::shared_ptr<Centroid> get_centroid_ptr() { return _centroid; }
    std::shared_ptr<Ellipse> get_ellipse_ptr() { return _ellipse; }
    std::shared_ptr<GaussianIntegral> get_integral_ptr() { return _integral; }

    const Centroid & get_centroid_const() const { return *_centroid;}
    const Ellipse & get_ellipse_const() const { return *_ellipse;}

    void set_const_normal(double const_normal) {
        _integral->set_value(get_const_normal()*2*_ellipse->get_area());
    }
    void set_integral_value(double integral) {
        _integral->set_value(integral);
    }

    void set_centroid_ptr(std::shared_ptr<Centroid> centroid) {
        this->_set_if_not_nullptr<Centroid>(_centroid, centroid, "centroid");
    }
    void set_ellipse_ptr(std::shared_ptr<Ellipse> ellipse) {
        this->_set_if_not_nullptr<Ellipse>(_ellipse, ellipse, "ellipse");
    }
    void set_integral_ptr(std::shared_ptr<GaussianIntegral> integral) {
        this->_set_if_not_nullptr<GaussianIntegral>(_integral, integral, "integral");
    }

    std::string str() const {
        return "Gaussian(centroid=" + _centroid->str() + ", ellipse=" + _ellipse->str() + ", integral=" + _integral->str() + ")";
    }

    Gaussian(std::shared_ptr<Centroid> centroid = nullptr, std::shared_ptr<Ellipse> ellipse = nullptr,
             std::shared_ptr<GaussianIntegral> integral = nullptr) :
        _centroid(centroid != nullptr ? std::move(centroid): std::make_shared<Centroid>()),
        _ellipse(ellipse != nullptr ? std::move(ellipse): std::make_shared<Ellipse>()),
        _integral(integral != nullptr ? std::move(integral): std::make_shared<GaussianIntegralValue>()
    ) {}
    ~Gaussian() {};
};

class Gaussians
{
public:
    typedef std::vector<std::shared_ptr<Gaussian>> Data;

private:
    Data _data = {};

    size_t assign(const Data * data, size_t i = 0)
    {
        size_t i_begin = i;
        for(const auto & gauss : *data)
        {
            if(gauss == nullptr) throw std::runtime_error("ConvolvedGaussian data["
                + std::to_string(i - i_begin) + "] can't be null");
            _data[i++] = gauss;
        }
        return i;
    }

public:
    Gaussian& operator[](size_t i) { return *(_data[i]); }
    const Gaussian& operator[](size_t i) const { return *(_data[i]); }

    typename Data::iterator begin() noexcept {return _data.begin();}
    typename Data::const_iterator cbegin() const noexcept {return _data.begin();}

    typename Data::iterator end() noexcept {return _data.end();}
    typename Data::const_iterator cend() const noexcept {return _data.cend();}

    Data get_data() const { return _data; }

    inline std::shared_ptr<Gaussian> & at(size_t i) {return _data.at(i);}
    size_t size() const {return _data.size();}

    std::string str() const {
        std::string s = "Gaussians([";
        for(const auto & g : _data) s += g->str() + ",";
        return s + "])";
    }

    Gaussians(const Data * data)
    {
        if(data != nullptr)
        {
            size_t n_data = data->size();
            if(n_data > 0)
            {
                _data.resize(n_data);
                this->assign(data);
            }
        }
    }
    Gaussians(std::vector<const Data *> data)
    {
        size_t n_data = 0;
        for(const auto & datum : data)
        {
            if(datum != nullptr) n_data += datum->size();
        }
        if(n_data > 0)
        {
            size_t i = 0;
            _data.resize(n_data);
            for(const auto & datum : data) i = this->assign(datum, i);
        }
    }
};

class ConvolvedGaussian
{
private:
    std::shared_ptr<Gaussian> _source;
    std::shared_ptr<Gaussian> _kernel;

public:
    Gaussian & get_source() { return *_source; }
    Gaussian & get_kernel() { return *_kernel; }

    const Gaussian & get_source_const() const { return *_source; }
    const Gaussian & get_kernel_const() const { return *_kernel; }

    std::string str() const {
        return "ConvolvedGaussian(source=" + _source->str() + ", kernel=" + _kernel->str() + ")";
    }

    ConvolvedGaussian(std::shared_ptr<Gaussian> source = nullptr, std::shared_ptr<Gaussian> kernel = nullptr) :
        _source(source != nullptr ? source : std::make_shared<Gaussian>()),
        _kernel(kernel != nullptr ? kernel : std::make_shared<Gaussian>())
    {}
};

class ConvolvedGaussians
{
public:
    typedef std::vector<std::shared_ptr<ConvolvedGaussian>> Data;

private:
    Data _data = {};

    size_t assign(const Data * data, size_t i = 0)
    {
        size_t i_begin = i;
        for(const auto & gauss : *data)
        {
            if(gauss == nullptr) throw std::runtime_error("ConvolvedGaussian data["
                + std::to_string(i - i_begin) + "] can't be null");
            _data[i++] = gauss;
        }
        return i;
    }

public:
    ConvolvedGaussian& operator[](size_t i) { return *(_data[i]); }
    const ConvolvedGaussian& operator[](size_t i) const { return *(_data[i]); }

    typename Data::iterator begin() noexcept {return _data.begin();}
    typename Data::const_iterator cbegin() const noexcept {return _data.begin();}

    typename Data::iterator end() noexcept {return _data.end();}
    typename Data::const_iterator cend() const noexcept {return _data.cend();}

    Data get_data() const { return _data; }

    inline std::shared_ptr<ConvolvedGaussian> & at(size_t i) {return _data.at(i);}
    size_t size() const {return _data.size();}

    std::string str() const {
        std::string s = "ConvolvedGaussians([";
        for(const auto & g : _data) s += g->str() + ",";
        return s + "])";
    }

    ConvolvedGaussians(const Data * data)
    {
        if(data != nullptr)
        {
            size_t n_data = data->size();
            if(n_data > 0)
            {
                _data.resize(n_data);
                this->assign(data);
            }
        }
    }
    ConvolvedGaussians(std::vector<const Data *> data)
    {
        size_t n_data = 0;
        for(const auto & datum : data)
        {
            if(datum != nullptr) n_data += datum->size();
        }
        if(n_data > 0)
        {
            size_t i = 0;
            _data.resize(n_data);
            for(const auto & datum : data) i = this->assign(datum, i);
        }
    }
};

} // namespace gauss2d
#endif
