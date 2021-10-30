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

#ifndef __GAUSS2D_GAUSSIAN_H_
#define __GAUSS2D_GAUSSIAN_H_

#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace gauss2d {

const double M_HWHM_SIGMA = 1.1774100225154746910115693264599;
const double M_SIGMA_HWHM = 0.84932180028801904272150283410295;
const double M_PI_180 = M_PI/180.;
const double M_180_PI = 180./M_PI;

class Ellipse;
class EllipseMajor;

class Centroid
{
public:
    double x = 0;
    double y = 0;

    std::string str() const {
        return "Centroid(x=" + std::to_string(x) + ", y=" + std::to_string(y) + ")";
    }

    Centroid(double x=0, double y=0) : x(x), y(y) {}
};

class Covariance
{
private:
    double _sigma_x_sq = 0;
    double _sigma_y_sq = 0;
    double _cov_xy = 0;

public:
    static void check(double sigma_x_sq, double sigma_y_sq, double cov_xy)
    {
        double offdiag_max = sqrt(sigma_x_sq) * sqrt(sigma_y_sq);
        // Define implied rho as -1 for negative cov and +1 for positive cov if zero size
        // This enforces cov_xy == 0 if sigma_x_sq == sigma_y_sq == 0
        double rho = offdiag_max > 0 ? cov_xy / offdiag_max : (cov_xy > 0) - (cov_xy < 0);
        if(!(sigma_x_sq >= 0) || !(sigma_y_sq >= 0) || !(rho > -1 && rho < 1))
        {
            throw std::invalid_argument(
                "Invalid sigma_x_sq, sigma_y_sq, cov_xy=" + std::to_string(sigma_x_sq) + ","
                + std::to_string(sigma_y_sq) + "," + std::to_string(cov_xy) + " with implied rho="
                + std::to_string(rho) + "; sigma_x,y_sq >= 0 and -1 < rho < 1 required."
            );
        }
    }
    void convolve(const Covariance& cov);
    double get_sigma_x_sq() const { return _sigma_x_sq;};
    double get_sigma_y_sq() const { return _sigma_y_sq;};
    double get_cov_xy() const { return _cov_xy;};

    std::shared_ptr<Covariance> make_convolution(const Covariance& cov) const;

    void set(const Ellipse & ellipse);
    void set(double sigma_x_sq=0, double sigma_y_sq=0, double cov_xy=0);
    void set_sigma_x_sq(double sigma_x_sq);
    void set_sigma_y_sq(double sigma_y_sq);
    void set_cov_xy(double cov_xy);
    std::string str() const {
        return "Covariance(sigma_x_sq=" + std::to_string(_sigma_x_sq) + ", sigma_y_sq="
            + std::to_string(_sigma_y_sq) + ", cov_xy=" + std::to_string(_cov_xy) + ")";
    }

    Covariance(double sigma_x_sq, double sigma_y_sq, double cov_xy);
    Covariance(const Ellipse & ell);
};

class EllipseMajor;
class Ellipse
{
private:
    struct Data
    {
        double _rho = 0.;
        double _sigma_x = 0.;
        double _sigma_y = 0.;
    };
    std::unique_ptr<Data> _data;

public:
    // Set without any checks
    virtual void _set(double sigma_x, double sigma_y, double rho) {
        _data->_sigma_x=sigma_x;
        _data->_sigma_y=sigma_y;
        _data->_rho=rho;
    }
    virtual void _set_rho(double rho) { _data->_rho = rho;};
    virtual void _set_sigma_x(double sigma_x) { _data->_sigma_x = sigma_x;};
    virtual void _set_sigma_y(double sigma_y) { _data->_sigma_y = sigma_y;};

    static void check(double sigma_x, double sigma_y, double rho)
    {
        if(!(sigma_x >= 0) || !(sigma_y >= 0) || !(rho > -1 && rho < 1))
        {
            throw std::invalid_argument(
                "Invalid sigma_x, sigma_y, rho=" + std::to_string(sigma_x) + ","
                + std::to_string(sigma_y) + "," + std::to_string(rho)
                + "; sigma_x,y >= 0 and -1 < rho < 1 required."
            );
        }
    }
    void convolve(const Ellipse& ell);

    double get_cov_xy() const;
    double get_radius_trace() const;
    double get_sigma_x_sq() const;
    double get_sigma_y_sq() const;
    virtual double get_rho() const { return _data->_rho;}
    virtual double get_sigma_x() const { return _data->_sigma_x; }
    virtual double get_sigma_y() const { return _data->_sigma_y; }
    virtual double get_sigma_xy() const { return _data->_sigma_x*_data->_sigma_y; }

    std::shared_ptr<Ellipse> make_convolution(const Ellipse& ell) const;

    virtual void set(double sigma_x, double sigma_y, double rho);
    virtual void set(const Covariance & covar);
    virtual void set(const EllipseMajor & ellipse);
    virtual void set_rho(double rho);
    virtual void set_sigma_x(double sigma_x);
    virtual void set_sigma_y(double sigma_y);

    std::string str() const {
        return "Ellipse(sigma_x=" + std::to_string(this->get_sigma_x())
            + ", sigma_y=" + std::to_string(this->get_sigma_y())
            + ", rho=" + std::to_string(this->get_rho()) + ")";
    }

    Ellipse(double sigma_x=0, double sigma_y=0, double rho=0);
    Ellipse(const Covariance & covar);
    Ellipse(const EllipseMajor & ellipse);
    virtual ~Ellipse() {};
};

class EllipseMajor
{
private:
    double _r_major = 0.;
    double _axrat = 1.;
    double _angle = 0.;
    bool _degrees = false;

public:
    static void check(double r_major, double axrat, double angle)
    {
        if(!(r_major >= 0) || !(axrat >= 0 && axrat <= 1))
        {
            throw std::invalid_argument(
                "Invalid r_major, axrat, angle=" + std::to_string(r_major) + "," + std::to_string(axrat)
                + "," + std::to_string(angle) + "; r_major >= 0, 1 >= axrat >= 0 required.");
        }
    }
    double get_r_major() const {return _r_major;}
    double get_axrat() const {return _axrat;}
    double get_angle() const {return _angle;}
    double get_angle_degrees() const {return _degrees ? _angle : _angle*M_180_PI;}
    double get_angle_radians() const {return _degrees ? _angle*M_PI_180 : _angle;}
    bool is_degrees() const {return _degrees;}

    void set(double r_major, double axrat, double angle);
    void set_r_major(double r_major);
    void set_axrat(double axrat);
    void set_angle(double angle);
    void set_degrees(bool degrees);
    std::string str() const {
        return "EllipseMajor(r_major=" + std::to_string(_r_major) + ", axrat=" + std::to_string(_axrat)
            + ", angle=" + std::to_string(_angle) + ", degrees=" + std::to_string(_degrees) + ")";
    }

    EllipseMajor(double r_major, double axrat, double angle, bool degrees=false);
    EllipseMajor(Covariance & covar, bool degrees=false);
    EllipseMajor(Ellipse & ellipse, bool degrees=false);
};

class Gaussian : public std::enable_shared_from_this<Gaussian>
{
private:
    double _integral;
    double _const_normal;

public:
    std::shared_ptr<Centroid> cen;
    std::shared_ptr<Ellipse> ell;

    double get_const_normal() const {return _const_normal;};
    double get_integral() const {return _integral;};

    void set_const_normal(double const_normal) {
        _const_normal = const_normal;
        double rho = ell->get_rho();
        _integral = _const_normal*2.*M_PI*ell->get_sigma_xy()*sqrt(1-rho*rho);
    }
    void set_integral(double integral) {
        _integral = integral;
        double rho = ell->get_rho();
        _const_normal = integral/(2.*M_PI*ell->get_sigma_xy()*sqrt(1-rho*rho));
    }

    std::string str() const {
        return "Gaussian(cen=" + cen->str() + ", ell=" + ell->str() + ", integral=" + std::to_string(_integral) + ")";
    }

    Gaussian(std::shared_ptr<Centroid> cen = nullptr, std::shared_ptr<Ellipse> ell = nullptr, double integral=1.) :
        cen(cen != nullptr ? std::move(cen): std::make_shared<Centroid>()),
        ell(ell != nullptr? std::move(ell): std::make_shared<Ellipse>())
    {
        set_integral(integral); 
    }
};

class ConvolvedGaussian
{
private:
    std::shared_ptr<Gaussian> _source;
    std::shared_ptr<Gaussian> _kernel;

public:
    inline Gaussian & source() { return *_source; }
    inline Gaussian & kernel() { return *_kernel; }

    inline const Gaussian & source_const() const { return *_source; }
    inline const Gaussian & kernel_const() const { return *_kernel; }

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
