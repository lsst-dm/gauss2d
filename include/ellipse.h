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
#include <stdexcept>

#include "object.h"

namespace gauss2d {

const double M_HWHM_SIGMA = 1.1774100225154746910115693264599;
const double M_SIGMA_HWHM = 0.84932180028801904272150283410295;
const double M_PI_180 = M_PI/180.;
const double M_180_PI = 180./M_PI;

class Ellipse;
class EllipseMajor;

/**
 * A Covariance is a representation of a 2D Gaussian with two
 * directional standard deviations (sigma_x, sigma_y) and a
 *
 * This object is intended for intermediate usage and is
 * implemented directly with double members rather than
 * having an abstract Data class.
**/
class Covariance : public Object
{
private:
    double _sigma_x_sq = 0;
    double _sigma_y_sq = 0;
    double _cov_xy = 0;

public:
    static void check(double sigma_x_sq, double sigma_y_sq, double cov_xy);
    void convolve(const Covariance & cov);
    double get_sigma_x_sq() const { return _sigma_x_sq;};
    double get_sigma_y_sq() const { return _sigma_y_sq;};
    double get_cov_xy() const { return _cov_xy;};
    std::array<double, 3> get_xyc() const { return {_sigma_x_sq, _sigma_y_sq, _cov_xy}; }

    std::shared_ptr<Covariance> make_convolution(const Covariance & cov) const;
    std::unique_ptr<Covariance> make_convolution_uniq(const Covariance & cov) const;

    void set(const Ellipse & ellipse);
    void set(double sigma_x_sq=0, double sigma_y_sq=0, double cov_xy=0);
    void set_sigma_x_sq(double sigma_x_sq);
    void set_sigma_y_sq(double sigma_y_sq);
    void set_cov_xy(double cov_xy);
    void set_xyc(const std::array<double, 3> & xyc);

    std::string str() const override;

    bool operator==(const Covariance& other) const;

    friend std::ostream & operator << (std::ostream &out, const Covariance &obj);

    Covariance(double sigma_x_sq=0, double sigma_y_sq=0, double cov_xy=0);
    Covariance(const Ellipse & ell);
};

/**
 * Interface for classes storing data defining an Ellipse.
**/
class EllipseData : public Object
{
public:
    virtual double get_sigma_x() const = 0;
    virtual double get_sigma_y() const = 0;
    virtual double get_rho() const = 0;
    virtual std::array<double, 3> get_xyr() const = 0;

    virtual void set_sigma_x(double sigma_x) = 0;
    virtual void set_sigma_y(double sigma_y) = 0;
    virtual void set_rho(double rho) = 0;
    virtual void set(double sigma_x, double sigma_y, double rho) = 0;
    virtual void set_xyr(const std::array<double, 3> & xyr) = 0;

    virtual std::string str() const override = 0;

    bool operator==(const EllipseData& other) const {
        return get_xyr() == other.get_xyr();
    };

    friend std::ostream & operator << (std::ostream &out, const EllipseData &obj) {
        out << obj.str();
        return out;
    }

    virtual ~EllipseData() = default;
};

/**
 * An EllipseData storing ellipse values as shared_ptrs.
 * 
 * shared_ptr usage allows EllipseValues to share parameters if desired. 
 *
**/
class EllipseValues : public EllipseData
{
private:
    std::shared_ptr<double> _sigma_x;
    std::shared_ptr<double> _sigma_y;
    std::shared_ptr<double> _rho;

public:
    double get_sigma_x() const override { return *_sigma_x; }
    double get_sigma_y() const override { return *_sigma_y; }
    double get_rho() const override { return *_rho; }
    std::array<double, 3> get_xyr() const override { return {*_sigma_x, *_sigma_y, *_rho}; }

    void set_sigma_x(double sigma_x) override;
    void set_sigma_y(double sigma_y) override;
    void set_rho(double rho) override;
    void set(double sigma_x, double sigma_y, double rho) override;
    void set_xyr(const std::array<double, 3> & xyr) override;

    std::string str() const override;

    EllipseValues(
        std::shared_ptr<double> sigma_x,
        std::shared_ptr<double> sigma_y,
        std::shared_ptr<double> rho=nullptr
    ) :
        _sigma_x(sigma_x == nullptr ? std::make_shared<double>(0) : std::move(sigma_x)),
        _sigma_y(sigma_y == nullptr ? std::make_shared<double>(0) : std::move(sigma_y)),
        _rho(rho == nullptr ? std::make_shared<double>(0) : std::move(rho)) {};
    EllipseValues(double sigma_x=0, double sigma_y=0, double rho=0) :
        _sigma_x(std::make_shared<double>(sigma_x)), _sigma_y(std::make_shared<double>(sigma_y)),
        _rho(std::make_shared<double>(rho)) {};
};

/**
 * An Ellipse is the representation of the shape of an ellipse
 * or its respective 2D Gaussian with the three unique terms 
 * of its covariance matrix (or a weighted moment of inertia
 * tensor). These terms are the squares of the xx and yy 
 * moments (equivalent to the squares of the x- and y-axis
 * standard deviations), and the (signed) xy moment (also known as 
 * the covariance).
 *
 * The storage of the parameters is implemented in EllipseData;
 * this class serves as a storage-independent interface for usage in
 * Gaussian classes.
 *
**/
class Ellipse : public Object
{
private:
    std::shared_ptr<EllipseData> _data;

public:
    static void check(double sigma_x, double sigma_y, double rho)
    {
        if(!(sigma_x >= 0) || !(sigma_y >= 0) || !(rho >= -1 && rho <= 1))
        {
            throw std::invalid_argument(
                "Invalid sigma_x, sigma_y, rho=" + std::to_string(sigma_x) + ","
                + std::to_string(sigma_y) + "," + std::to_string(rho)
                + "; sigma_x,y >= 0 and 1 >= rho >= -1 required."
            );
        }
    }
    void convolve(const Ellipse& ell);

    double get_area() const;
    double get_cov_xy() const;
    const EllipseData & get_data() const { return *_data; }
    double get_radius_trace() const;
    double get_sigma_x_sq() const;
    double get_sigma_y_sq() const;
    double get_rho() const { return _data->get_rho();}
    double get_sigma_x() const { return _data->get_sigma_x(); }
    double get_sigma_y() const { return _data->get_sigma_y(); }
    double get_sigma_xy() const { return _data->get_sigma_x()*_data->get_sigma_y(); }
    std::array<double, 3> get_xyr() const { return  _data->get_xyr(); }

    std::shared_ptr<Ellipse> make_convolution(const Ellipse& ell) const;
    std::unique_ptr<Ellipse> make_convolution_uniq(const Ellipse & ell) const;

    void set(double sigma_x, double sigma_y, double rho);
    void set(const Covariance & covar);
    void set(const EllipseMajor & ellipse);
    void set_rho(double rho);
    void set_sigma_x(double sigma_x);
    void set_sigma_y(double sigma_y);
    void set_xyr(const std::array<double, 3> & xyr);

    std::string str() const override {
        return "Ellipse(data=" + _data->str() + ")";
    }

    bool operator==(const Ellipse& other) const {
        return get_data() == other.get_data();
    }

    Ellipse(std::shared_ptr<EllipseData> data);
    Ellipse(double sigma_x=0, double sigma_y=0, double rho=0);
    Ellipse(const Covariance & covar);
    Ellipse(const EllipseMajor & ellipse);
    ~Ellipse() {};
};

/**
 * An EllipseMajor is a representation of a 2D Gaussian with a
 * major axis length (in units of the standard deviation),
 * the axis ratio (minor axis divided by major), and a 
 * position angle, by convention here defined as counter-clockwise
 * from the positive-x axis.
 *
 * This object is intended for intermediate usage and is
 * implemented directly with double members rather than
 * having an abstract Data class.
**/
class EllipseMajor : public Object
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
    std::array<double, 3> get_rqa() const { return  {_r_major, _axrat, _angle}; }
    bool is_degrees() const {return _degrees;}

    void set(double r_major, double axrat, double angle);
    void set_r_major(double r_major);
    void set_axrat(double axrat);
    void set_angle(double angle);
    void set_degrees(bool degrees);
    void set_rqa(const std::array<double, 3> & rqa);

    std::string str() const override {
        return "EllipseMajor(r_major=" + std::to_string(_r_major) + ", axrat=" + std::to_string(_axrat)
            + ", angle=" + std::to_string(_angle) + ", degrees=" + std::to_string(_degrees) + ")";
    }

    bool operator==(const EllipseMajor& other) const {
        return (get_r_major() == other.get_r_major()) && (get_axrat() == other.get_axrat())
            && (get_angle_degrees() == other.get_angle_degrees());
    }

    EllipseMajor(double r_major, double axrat, double angle, bool degrees=false);
    EllipseMajor(Covariance & covar, bool degrees=false);
    EllipseMajor(Ellipse & ellipse, bool degrees=false);
};

}
#endif