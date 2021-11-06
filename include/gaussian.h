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
    std::shared_ptr<Centroid> _cen;
    std::shared_ptr<Ellipse> _ell;
    std::shared_ptr<GaussianIntegral> _integral;

public:
    double get_const_normal() const { return _integral->get_value()/(2*_ell->get_area()); }
    double get_integral() const {return _integral->get_value();};

    Centroid & get_centroid() { return *_cen;}
    Ellipse & get_ellipse() { return *_ell;}

    const Centroid & get_centroid_const() const { return *_cen;}
    const Ellipse & get_ellipse_const() const { return *_ell;}

    void set_const_normal(double const_normal) { _integral->set_value(get_const_normal()*2*_ell->get_area()); }
    void set_integral(double integral) {
        _integral->set_value(integral);
    }

    std::string str() const {
        return "Gaussian(cen=" + _cen->str() + ", ell=" + _ell->str() + ", integral=" + _integral->str() + ")";
    }

    Gaussian(std::shared_ptr<Centroid> cen = nullptr, std::shared_ptr<Ellipse> ell = nullptr,
             std::shared_ptr<GaussianIntegral> integral = nullptr) :
        _cen(cen != nullptr ? std::move(cen): std::make_shared<Centroid>()),
        _ell(ell != nullptr ? std::move(ell): std::make_shared<Ellipse>()),
        _integral(integral != nullptr ? std::move(integral): std::make_shared<GaussianIntegralValue>()
    ) {}
    ~Gaussian() {};
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

class Gaussians
{
public:
    typedef std::vector<std::shared_ptr<ConvolvedGaussian>> Data;

private:
    Data _data = {};

public:
    ConvolvedGaussian& operator[](size_t i) {return *(_data[i]);}
    const ConvolvedGaussian& operator[](size_t i) const {return *(_data[i]);}

    typename Data::iterator begin() noexcept {return _data.begin();}
    typename Data::const_iterator cbegin() const noexcept {return _data.begin();}

    typename Data::iterator end() noexcept {return _data.end();}
    typename Data::const_iterator cend() const noexcept {return _data.cend();}

    inline std::shared_ptr<ConvolvedGaussian> & at(size_t i) {return _data.at(i);}
    size_t size() const {return _data.size();}

    std::string str() const {
        std::string s = "Gaussians([";
        for(const auto & g : _data) s += g->str() + ",";
        return s + "])";
    }

    Gaussians(const Data * data_in)
    {
        if(data_in != nullptr)
        {
            const Data & data = *data_in;
            size_t n_data = data.size();
            if(n_data > 0)
            {
                _data.resize(n_data);
                for(size_t i = 0; i < n_data; ++i)
                {
                    if(data[i] == nullptr) throw std::runtime_error("ConvolvedGaussian data[" + std::to_string(i) + "] can't be null");
                    _data[i] = data[i];
                }
            }
        }
    }
};

} // namespace gauss2d
#endif
