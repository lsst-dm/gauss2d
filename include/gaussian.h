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

//#include <pybind11/pybind11.h>
//#include <pybind11/numpy.h>

//namespace py = pybind11;

//typedef py::array_t<double> ndarray;

#include <memory>
#include <utility>
#include <cmath>

namespace gauss2d {
    const double M_HWHM_SIGMA = 1.1774100225154746910115693264599;
    const double M_SIGMA_HWHM = 0.84932180028801904272150283410295;
    const double M_PI_180 = M_PI/180.;
    const double M_180_PI = 180./M_PI;

    class Ellipse;
    class EllipseMajor;
    class EllipseTerms;

    class Centroid
    {
    public:
        double x;
        double y;

        Centroid(double x, double y) : x(x), y(y) {}
    };

    class Covariance
    {
    private:
        double _sigma_x_sq;
        double _sigma_y_sq;
        double _cov_xy;
        void init(double sigma_x, double sigma_y, double rho);

    public:
        double get_sigma_x_sq() const { return _sigma_x_sq;};
        double get_sigma_y_sq() const { return _sigma_y_sq;};
        double get_cov_xy() const { return _cov_xy;};

        std::shared_ptr<EllipseMajor> make_ellipse_major(bool degrees=false) const;

        void set_sigma_x_sq(double sigma_x_sq);
        void set_sigma_y_sq(double sigma_y_sq);
        void set_cov_xy(double cov_xy);

        Covariance(double sigma_x, double sigma_y, double rho) { this->init(sigma_x, sigma_y, rho);};
        Covariance(const Ellipse & ell);
    };

    class Ellipse
    {
    private:
        double _rho;
        double _sigma_x;
        double _sigma_y;

    public:
        // Set without any checks
        inline void _set(double sigma_x, double sigma_y, double rho) {
            _sigma_x=sigma_x;
            _sigma_y=sigma_y;
            _rho=rho;
        }
        inline void _set_rho(double rho) { _rho = rho;};
        inline void _set_sigma_x(double sigma_x) { _sigma_x = sigma_x;};
        inline void _set_sigma_y(double sigma_y) { _sigma_y = sigma_y;};

        static void check(double sigma_x, double sigma_y, double rho)
        {
            if(!(rho > -1 && rho < 1) || !(sigma_x >= 0) || !(sigma_y >= 0))
            {
                throw std::invalid_argument(
                    "Invalid rho, sigma_x, sigma_y=" + std::to_string(rho) + "," + std::to_string(rho) +
                    "," + std::to_string(rho) + "; -1 < rho < 1, sigma_x,y >= 0 required.");
            }
        }
        void convolve(const Ellipse& ell);

        double get_cov_xy() const { return _sigma_x * _sigma_y * _rho;};
        double get_radius_trace() const { return sqrt(_sigma_x*_sigma_x + _sigma_y*_sigma_y);};
        EllipseTerms get() const;
        inline double get_rho() const { return _rho;};
        inline double get_sigma_x() const { return _sigma_x; };
        inline double get_sigma_y() const { return _sigma_y; };

        std::shared_ptr<Ellipse> make_convolution(const Ellipse& ell) const;
        std::shared_ptr<EllipseMajor> make_ellipse_major(bool degrees=false) const;

        void set(double sigma_x, double sigma_y, double rho);
        void set_rho(double rho);
        void set_sigma_x(double sigma_x);
        void set_sigma_y(double sigma_y);

        Ellipse(double sigma_x=0, double sigma_y=0, double rho=0)
        {
            this->set(sigma_x, sigma_y, rho);
        }
        Ellipse(EllipseTerms & terms);
    };

    class EllipseMajor
    {
    private:
        double _r_major;
        double _axrat;
        double _angle;
        bool _degrees;
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
        bool is_degrees() const {return _degrees;}

        void set_r_major(double r_major);
        void set_axrat(double axrat);
        void set_angle(double angle);
        void set_degrees(bool degrees);

        EllipseMajor(double r_major, double axrat, double angle, bool degrees=false);
    };

    class EllipseTerms
    {
    public:
        double sigma_x;
        double sigma_y;
        double rho;

        EllipseTerms(double sigma_x, double sigma_y, double rho): sigma_x(sigma_x), sigma_y(sigma_y), rho(rho)
        {}
    };

    class Gaussian
    {
    public:
        std::shared_ptr<Centroid> cen;
        std::shared_ptr<Ellipse> ell;

        Gaussian(std::shared_ptr<Centroid> cen, std::shared_ptr<Ellipse> ell) :
            cen(std::move(cen)), ell(std::move(ell)) {}
    };
}
#endif
