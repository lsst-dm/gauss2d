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

#ifndef GAUSS2D_ELLIPSE_H
#define GAUSS2D_ELLIPSE_H

#include <cmath>
#include <memory>

namespace gauss2d {

const double M_HWHM_SIGMA = 1.1774100225154746910115693264599;
const double M_SIGMA_HWHM = 0.84932180028801904272150283410295;
const double M_PI_180 = M_PI/180.;
const double M_180_PI = 180./M_PI;

class Ellipse;
class EllipseMajor;

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

class EllipseData
{
public:
    virtual double get_sigma_x() const = 0;
    virtual double get_sigma_y() const = 0;
    virtual double get_rho() const = 0;

    virtual void set_sigma_x(double sigma_x) = 0;
    virtual void set_sigma_y(double sigma_y) = 0;
    virtual void set_rho(double rho) = 0;
    virtual void set(double sigma_x, double sigma_y, double rho) = 0;

    virtual std::string str() const = 0;
    virtual ~EllipseData() = default;
};

class EllipseValues : public EllipseData
{
private:
    std::shared_ptr<double> _sigma_x;
    std::shared_ptr<double> _sigma_y;
    std::shared_ptr<double> _rho;

public:
    double get_sigma_x() const { return *_sigma_x; }
    double get_sigma_y() const { return *_sigma_y; }
    double get_rho() const { return *_rho; }

    void set_sigma_x(double sigma_x);
    void set_sigma_y(double sigma_y);
    void set_rho(double rho);
    void set(double sigma_x, double sigma_y, double rho);

    std::string str() const;

    EllipseValues(std::shared_ptr<double> sigma_x, std::shared_ptr<double> sigma_y, std::shared_ptr<double> rho=nullptr) :
        _sigma_x(sigma_x == nullptr ? std::make_shared<double>(0) : std::move(sigma_x)),
        _sigma_y(sigma_y == nullptr ? std::make_shared<double>(0) : std::move(sigma_y)),
        _rho(rho == nullptr ? std::make_shared<double>(0) : std::move(rho)) {};
    EllipseValues(double sigma_x=0, double sigma_y=0, double rho=0) :
        _sigma_x(std::make_shared<double>(sigma_x)), _sigma_y(std::make_shared<double>(sigma_y)),
        _rho(std::make_shared<double>(rho)) {};
};

class Ellipse
{
private:
    std::shared_ptr<EllipseData> _data;

public:
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

    double get_area() const;
    double get_cov_xy() const;
    double get_radius_trace() const;
    double get_sigma_x_sq() const;
    double get_sigma_y_sq() const;
    double get_rho() const { return _data->get_rho();}
    double get_sigma_x() const { return _data->get_sigma_x(); }
    double get_sigma_y() const { return _data->get_sigma_y(); }
    double get_sigma_xy() const { return _data->get_sigma_x()*_data->get_sigma_y(); }

    std::shared_ptr<Ellipse> make_convolution(const Ellipse& ell) const;

    void set(double sigma_x, double sigma_y, double rho);
    void set(const Covariance & covar);
    void set(const EllipseMajor & ellipse);
    void set_rho(double rho);
    void set_sigma_x(double sigma_x);
    void set_sigma_y(double sigma_y);

    std::string str() const {
        return "Ellipse(" + _data->str() + ")";
    }

    Ellipse(std::shared_ptr<EllipseData> data);
    Ellipse(double sigma_x=0, double sigma_y=0, double rho=0);
    Ellipse(const Covariance & covar);
    Ellipse(const EllipseMajor & ellipse);
    ~Ellipse() {};
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
    double get_area() const {return M_PI*_r_major*_r_major*_axrat;}
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

}
#endif