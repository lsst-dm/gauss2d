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
#include "ellipse.h"

#include <cmath>
#include <stdexcept>
#include <string>

namespace gauss2d
{
    template <class S>
    std::pair<S,S> sincos(S arg) { return { std::sin(arg), std::cos(arg) }; }

    Covariance::Covariance(double sigma_x_sq, double sigma_y_sq, double cov_xy)
    {
        set(sigma_x_sq, sigma_y_sq, cov_xy);
    }

    Covariance::Covariance(const Ellipse & ell)
    {
        set(ell);
    }

    void Covariance::check(double sigma_x_sq, double sigma_y_sq, double cov_xy)
    {
        double offdiag_max = sqrt(sigma_x_sq) * sqrt(sigma_y_sq);
        // Define implied rho as -1 for negative cov and +1 for positive cov if zero size
        // This enforces cov_xy == 0 if sigma_x_sq == sigma_y_sq == 0
        double rho = offdiag_max > 0 ? cov_xy / offdiag_max : (cov_xy > 0) - (cov_xy < 0);
        if(!(sigma_x_sq >= 0) || !(sigma_y_sq >= 0) || !(rho >= -1 && rho <= 1))
        {
            throw std::invalid_argument(
                "Invalid sigma_x_sq, sigma_y_sq, cov_xy=" + std::to_string(sigma_x_sq) + ","
                + std::to_string(sigma_y_sq) + "," + std::to_string(cov_xy) + " with implied rho="
                + std::to_string(rho) + "; sigma_x,y_sq >= 0 and -1 < rho < 1 required."
            );
        }
    }

    void Covariance::convolve(const Covariance& cov)
    {
        double sigma_x_sq = this->get_sigma_x_sq() + cov.get_sigma_x_sq();
        double sigma_y_sq = this->get_sigma_y_sq() + cov.get_sigma_y_sq();
        double cov_xy = this->get_cov_xy() + cov.get_cov_xy();
        this->set(sigma_x_sq, sigma_y_sq, cov_xy);
    }

    std::shared_ptr<Covariance> Covariance::make_convolution(const Covariance& cov) const
    {
        std::shared_ptr<Covariance> cov_ret = std::make_shared<Covariance>(
            this->get_sigma_x_sq(), this->get_sigma_y_sq(), this->get_cov_xy());
        cov_ret->convolve(cov);
        return cov_ret;
    }

    void Covariance::set(const Ellipse & ellipse)
    {
        double sigma_x = ellipse.get_sigma_x();
        double sigma_y = ellipse.get_sigma_y();
        double rho = ellipse.get_rho();
        this->_sigma_x_sq = sigma_x*sigma_x;
        this->_sigma_y_sq = sigma_y*sigma_y;
        // Must be done like so in case sigma_x/y^2 round to zero, and check after
        this->set_cov_xy(sigma_x * sigma_y * (_sigma_x_sq > 0) * (_sigma_y_sq > 0) * rho);
    }

    void Covariance::set(double sigma_x_sq, double sigma_y_sq, double cov_xy)
    {
        set_sigma_x_sq(sigma_x_sq);
        set_sigma_y_sq(sigma_y_sq);
        set_cov_xy(cov_xy);
    }

    void Covariance::set_sigma_x_sq(double sigma_x_sq) {
        if(!(sigma_x_sq >= 0))
        {
            throw std::invalid_argument(
                this->str() + " can't set invalid sigma_x_sq=" + std::to_string(sigma_x_sq)
                 + "; sigma_x_sq >= 0 required.");
        }
        _sigma_x_sq = sigma_x_sq;
    }

    void Covariance::set_sigma_y_sq(double sigma_y_sq) {
        if(!(sigma_y_sq >= 0))
        {
            throw std::invalid_argument(
                this->str() + " can't set invalid sigma_y_sq=" + std::to_string(sigma_y_sq)
                 + "; sigma_y_sq >= 0 required.");
        }
        _sigma_y_sq = sigma_y_sq;
    }

    void Covariance::set_cov_xy(double cov_xy) {
        // Take individual sqrt just to be safe and avoid potential overflow
        double offdiag_max = sqrt(_sigma_x_sq) * sqrt(_sigma_y_sq);
        // If offdiag_max is zero, we can only accept cov_xy exactly zero
        double rho = offdiag_max > 0 ? cov_xy / offdiag_max : -!(cov_xy >= 0) + !(cov_xy <= 0);
        if(!(rho > -1 && rho < 1))
        {
            throw std::invalid_argument(
                this->str() + "can't set invalid cov_xy=" + std::to_string(cov_xy)
                + " (>0=" + std::to_string(cov_xy > 0)
                + ", <0=" + std::to_string(cov_xy < 0)
                + ", offdiag_max=" + std::to_string(offdiag_max)
                + " with implied rho=" + std::to_string(rho) + "; -1 < rho < 1 required.");
        }
        _cov_xy = cov_xy;
    }

    void Covariance::set_xyc(const std::array<double, 3> & xyc) {
        this->set(xyc[0], xyc[1], xyc[2]);
    }

    std::string Covariance::str() const {
        return "Covariance(sigma_x_sq=" + std::to_string(_sigma_x_sq) + ", sigma_y_sq="
            + std::to_string(_sigma_y_sq) + ", cov_xy=" + std::to_string(_cov_xy) + ")";
    }

    bool Covariance::operator == (const Covariance& other) const {
        return get_xyc() == other.get_xyc();
    }

    std::ostream & operator << (std::ostream &out, const Covariance &obj) {
        out << obj.str();
        return out;
    }

    void EllipseValues::set_sigma_x(double sigma_x) {
        if(!(_sigma_x >= 0))
        {
            throw std::invalid_argument("Invalid sigma_x=" + std::to_string(sigma_x) +
                                        "; sigma_y >= 0 required.");
        }
         *_sigma_x = sigma_x;
    }
    
    void EllipseValues::set_sigma_y(double sigma_y) {
        if(!(sigma_y >= 0))
        {
            throw std::invalid_argument("Invalid sigma_y=" + std::to_string(sigma_y) +
                                        "; sigma_y >= 0 required.");
        }
        *_sigma_y = sigma_y;
    }

    void EllipseValues::set_rho(double rho) {
        if(!(rho > -1 && rho < 1))
        {
            throw std::invalid_argument("Invalid rho=" + std::to_string(rho) + "; -1 < rho < 1 required.");
        }
        *_rho = rho;
    }

    void EllipseValues::set(double sigma_x, double sigma_y, double rho) {
        Ellipse::check(sigma_x, sigma_y, rho);
        *_sigma_x = sigma_x;
        *_sigma_y = sigma_y;
        *_rho = rho;
    }

    void EllipseValues::set_h(double hwhm_x, double hwhm_y, double rho) {
        Ellipse::check(hwhm_x, hwhm_y, rho);
        *_sigma_x = M_HWHM_SIGMA*hwhm_x;
        *_sigma_y = M_HWHM_SIGMA*hwhm_y;
        *_rho = rho;
    }

    void EllipseValues::set_hxyr(const std::array<double, 3> & hxyr) {
        this->set_h(hxyr[0], hxyr[1], hxyr[2]);
    }

    void EllipseValues::set_xyr(const std::array<double, 3> & xyr) {
        this->set(xyr[0], xyr[1], xyr[2]);
    }

    std::string EllipseValues::str() const {
        return "EllipseValues(sigma_x=" + std::to_string(this->get_sigma_x())
            + ", sigma_y=" + std::to_string(this->get_sigma_y())
            + ", rho=" + std::to_string(this->get_rho()) + ")";
    }

    void Ellipse::convolve(const Ellipse& ell)
    {
        double sigma_x_ell = ell.get_sigma_x();
        double sigma_y_ell = ell.get_sigma_y();
        double sigma_x = this->get_sigma_x();
        sigma_x = sqrt(sigma_x*sigma_x + sigma_x_ell*sigma_x_ell);
        double sigma_y = this->get_sigma_y();
        sigma_y = sqrt(sigma_y*sigma_y + sigma_y_ell*sigma_y_ell);
        double rho = (this->get_cov_xy() + ell.get_cov_xy()) / (sigma_x * sigma_y);
        this->set(sigma_x, sigma_y, rho);
    }

    double Ellipse::get_area() const {
        double rho = this->get_rho();
        return M_PI*this->get_sigma_xy()*sqrt(1-rho*rho);
    }

    double Ellipse::get_cov_xy() const {
        return this->get_sigma_x() * this->get_sigma_y() * this->get_rho();
    }

    double Ellipse::get_radius_trace() const {
        return sqrt(this->get_sigma_x_sq() + this->get_sigma_y_sq());
    }

    double Ellipse::get_sigma_x_sq() const
    {
        double sigma = this->get_sigma_x();
        return sigma*sigma;
    }

    double Ellipse::get_sigma_y_sq() const
    {
        double sigma = this->get_sigma_y();
        return sigma*sigma;
    }

    std::shared_ptr<Ellipse> Ellipse::make_convolution(const Ellipse& ell) const
    {
        return this->make_convolution_uniq(ell);
    }
    std::unique_ptr<Ellipse> Ellipse::make_convolution_uniq(const Ellipse& ell) const
    {
        // TODO: Replace with cloning derived data
        std::unique_ptr<Ellipse> ell_ret = std::make_unique<Ellipse>(
            this->get_sigma_x(), this->get_sigma_y(), this->get_rho());
        ell_ret->convolve(ell);
        return ell_ret;
    }

    void Ellipse::set(double sigma_x, double sigma_y, double rho) {
        _data->set(sigma_x, sigma_y, rho);
    }

    void Ellipse::set(const Covariance & covar)
    {
        double sigma_x = covar.get_sigma_x_sq();
        double sigma_y = covar.get_sigma_y_sq();
        if(sigma_x == 0 && sigma_y == 0) return;
        sigma_x = sqrt(sigma_x);
        sigma_y = sqrt(sigma_y);
        double rho = (sigma_x == 0 || sigma_y == 0) ? 0 : covar.get_cov_xy()/(sigma_x*sigma_y);
        this->set(sigma_x, sigma_y, rho);
    }

    void Ellipse::set(const EllipseMajor & ellipse)
    {
        const double r_major = ellipse.get_r_major();
        if(r_major == 0) return;
        const double axrat = ellipse.get_axrat();
        if(axrat == 1)
        {
            this->set_sigma_x(r_major);
            this->set_sigma_y(r_major);
            return;
        }
        const auto [sin_th, cos_th] = sincos(ellipse.get_angle_radians());
        const double sin_th_sq = sin_th*sin_th;
        const double cos_th_sq = cos_th*cos_th;

        const double r_major_sq = r_major*r_major;
        const double r_minor_sq = r_major_sq*axrat*axrat;
        const double sigma_x = sqrt(cos_th_sq*r_major_sq + sin_th_sq*r_minor_sq);
        const double sigma_y = sqrt(sin_th_sq*r_major_sq + cos_th_sq*r_minor_sq);
        double rho = (sigma_x == 0 || sigma_y == 0) ? 0
            : sin_th*cos_th*(r_major_sq - r_minor_sq)/(sigma_x*sigma_y);
        this->set(sigma_x, sigma_y, rho);
    }

    void Ellipse::set_h(double hwhm_x, double hwhm_y, double rho) {
        _data->set_h(hwhm_x, hwhm_y, rho);
    }
    void Ellipse::set_hwhm_x(double hwhm_x) {
        _data->set_hwhm_x(hwhm_x);
    }
    void Ellipse::set_hwhm_y(double hwhm_y) {
        _data->set_hwhm_y(hwhm_y);
    }

    void Ellipse::set_rho(double rho) {
        _data->set_rho(rho);
    }

    void Ellipse::set_sigma_x(double sigma_x) {
        _data->set_sigma_x(sigma_x);
    }

    void Ellipse::set_sigma_y(double sigma_y) {
        _data->set_sigma_y(sigma_y);
    }

    void Ellipse::set_hxyr(const std::array<double, 3> & hxyr) {
        this->set_h(hxyr[0], hxyr[1], hxyr[2]);
    }

    void Ellipse::set_xyr(const std::array<double, 3> & xyr) {
        this->set(xyr[0], xyr[1], xyr[2]);
    }

    Ellipse::Ellipse(std::shared_ptr<EllipseData> data) :
        _data(data == nullptr ? std::make_shared<EllipseValues>() : std::move(data)) {};

    Ellipse::Ellipse(double sigma_x, double sigma_y, double rho) :
        _data(std::make_shared<EllipseValues>())
    {
        _data->set(sigma_x, sigma_y, rho);
    }

    Ellipse::Ellipse(const Covariance & covar) :
        _data(std::make_shared<EllipseValues>())
    {
        this->set(covar);
    }

    Ellipse::Ellipse(const EllipseMajor & ellipse) :
        _data(std::make_shared<EllipseValues>())
    {
        this->set(ellipse);
    }

    std::pair<double, double> get_x_pm(double sigma_x_sq, double sigma_y_sq, double cov_xy)
    {
        double apc = sigma_x_sq + sigma_y_sq;
        double x = apc/2;
        // TODO: Write tests for this with inputs returning close to zero
        // Probably most efficient but unstable, e.g.:
        // sigma_x, sigma_y, rho = 1.58113883008419, 1.5811388300841895, 0
        // ... yields -3.552713678800501e-15 inside the sqrt
        // double pm = sqrt(apc*apc - 4*(sigma_x_sq*sigma_y_sq - cov_xy*cov_xy))/2;
        // Two more multiplications, but seemingly more stable
        // Cancels out cross term from apc (2*sigma_x_sq*sigma_y_sq)
        double pm = sqrt(sigma_x_sq*sigma_x_sq + sigma_y_sq*sigma_y_sq
            - 2*(sigma_x_sq*sigma_y_sq - 2*cov_xy*cov_xy))/2;
        return {x, pm};
    }

    void init(EllipseMajor & ellipse, const Covariance & covar, bool degrees)
    {
        double sigma_x_sq = covar.get_sigma_x_sq();
        double sigma_y_sq = covar.get_sigma_y_sq();
        if(sigma_x_sq == 0 && sigma_y_sq == 0) return;
        double cov_xy = covar.get_cov_xy();
        auto [x, pm] = get_x_pm(sigma_x_sq, sigma_y_sq, cov_xy);
        double r_major = x + pm;
        if(r_major == 0) return;

        double axrat = sqrt((x - pm)/r_major);
        r_major = sqrt(r_major);
        double ang = atan2(2*cov_xy, sigma_x_sq - sigma_y_sq)/2;
        if(degrees) ang *= M_180_PI;
        ellipse.set(r_major, axrat, ang);
    }

    EllipseMajor::EllipseMajor(double r_major, double axrat, double angle, bool degrees) :
        _r_major(r_major), _axrat(axrat), _angle(angle), _degrees(degrees)
    {
        EllipseMajor::check(r_major, axrat, angle);
    }

    EllipseMajor::EllipseMajor(Covariance & covar, bool degrees) : _degrees(degrees)
    {
        init(*this, covar, degrees);
    }

    EllipseMajor::EllipseMajor(Ellipse & ellipse, bool degrees) : _degrees(degrees)
    {
        init(*this, Covariance(ellipse), degrees);
    }

    void EllipseMajor::set(double r_major, double axrat, double angle)
    {
        check(r_major, axrat, angle);
        _r_major = r_major;
        _axrat = axrat;
        _angle = angle;
    }

    void EllipseMajor::set_r_major(double r_major)
    {
        if(!(r_major >= 0))
        {
            throw std::invalid_argument("Invalid r_major=" + std::to_string(r_major)
                + "; r_major >= 0 required.");
        }
        _r_major = r_major;
    }

    void EllipseMajor::set_axrat(double axrat)
    {
        if(!(axrat >= 0 && axrat <= 1))
        {
            throw std::invalid_argument("Invalid axrat=" + std::to_string(axrat)
                + "; 1 >= axrat >= 0 required.");
        }
        _axrat = axrat;
    }

    void EllipseMajor::set_angle(double angle) {
        _angle = angle;
    }

    void EllipseMajor::set_degrees(bool degrees) {
        if(degrees && !_degrees) {
            _degrees = true;
            _angle *= M_180_PI;
        } else if(!degrees && _degrees) {
            _degrees = false;
            _angle *= M_PI_180;
        }
    }

    void EllipseMajor::set_rqa(const std::array<double, 3> & rqa) {
        this->set(rqa[0], rqa[1], rqa[2]);
    }

}

#endif