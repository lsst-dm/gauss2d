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
#include "gaussian.h"

#include <string>
#include <stdexcept>
#include <cmath>

namespace gauss2d
{
    void Covariance::init(double sigma_x, double sigma_y, double rho) {
        Ellipse::check(sigma_x, sigma_y, rho);
        this->_sigma_x_sq = sigma_x*sigma_x;
        this->_sigma_y_sq = sigma_y*sigma_y;
        this->_cov_xy = sigma_x * sigma_y * rho;
    }

    Covariance::Covariance(const Ellipse & ell) {
        this->init(ell.get_sigma_x(), ell.get_sigma_y(), ell.get_rho());
    }

    std::shared_ptr<EllipseMajor> Covariance::make_ellipse_major(bool degrees) const
    {
        if(_sigma_x_sq == 0 && _sigma_y_sq == 0) return std::make_shared<EllipseMajor>(0, 1, 0);
        double apc = _sigma_x_sq + _sigma_y_sq;
        double x = apc/2;
        double pm = sqrt(apc*apc - 4*(_sigma_x_sq*_sigma_y_sq - _cov_xy*_cov_xy)) / 2;

        double r_major = x + pm;
        if(r_major == 0) return std::make_shared<EllipseMajor>(0, 1, 0);
        double axrat = sqrt((x - pm)/r_major);
        r_major = sqrt(r_major);
        double ang = atan2(2 * _cov_xy, _sigma_x_sq - _sigma_y_sq) / 2;
        if(degrees) ang *= M_180_PI;

        return std::make_shared<EllipseMajor>(r_major, axrat, ang, degrees);
    }

    void Covariance::set_sigma_x_sq(double sigma_x_sq) {
        if(!(sigma_x_sq >= 0))
        {
            throw std::invalid_argument("Invalid sigma_x_sq=" + std::to_string(sigma_x_sq) +
                                        "; sigma_x_sq >= 0 required.");
        }
        _sigma_x_sq = sigma_x_sq;
    }

    void Covariance::set_sigma_y_sq(double sigma_y_sq) {
        if(!(sigma_y_sq >= 0))
        {
            throw std::invalid_argument("Invalid sigma_y_sq=" + std::to_string(sigma_y_sq) +
                                        "; sigma_y_sq >= 0 required.");
        }
        _sigma_y_sq = sigma_y_sq;
    }

    void Covariance::set_cov_xy(double cov_xy) {
        // Take individual sqrt just to be safe and avoid potential overflow
        double offdiag_max = sqrt(_sigma_x_sq) * sqrt(_sigma_y_sq);
        double rho = cov_xy / offdiag_max;
        if(!(rho > -1 && rho < 1))
        {
            throw std::invalid_argument("Invalid cov_xy=" + std::to_string(cov_xy) + " with implied rho="
                                        + std::to_string(rho) + "; -1 < rho < 1 required.");
        }
        _cov_xy = cov_xy;
    }

    void Ellipse::convolve(const Ellipse& ell)
    {
        double sigma_x_ell = ell.get_sigma_x();
        double sigma_y_ell = ell.get_sigma_y();
        double sigma_x = sqrt(_sigma_x*_sigma_x + sigma_x_ell*sigma_x_ell);
        double sigma_y = sqrt(_sigma_y*_sigma_y + sigma_y_ell*sigma_y_ell);
        double rho = (this->get_cov_xy() + ell.get_cov_xy()) / (sigma_x * sigma_y);
        this->set(sigma_x, sigma_y, rho);
    }

    EllipseTerms Ellipse::get() const
    {
        return EllipseTerms(_sigma_x, _sigma_y, _rho);
    };

    std::shared_ptr<Ellipse> Ellipse::make_convolution(const Ellipse& ell) const
    {
        std::shared_ptr<Ellipse> ell_ret = std::make_shared<Ellipse>(_sigma_x, _sigma_y, _rho);
        ell_ret->convolve(ell);
        return std::move(ell_ret);
    }

    std::shared_ptr<EllipseMajor> Ellipse::make_ellipse_major(bool degrees) const
    {
        return std::move(Covariance(*this).make_ellipse_major(degrees=degrees));
    }

    void Ellipse::set(double sigma_x, double sigma_y, double rho) {
        Ellipse::check(sigma_x, sigma_y, rho);
        _rho = rho;
        _sigma_x = sigma_x;
        _sigma_y = sigma_y;
    }

    void Ellipse::set_rho(double rho) {
        if(!(rho > -1 && rho < 1))
        {
            throw std::invalid_argument("Invalid rho=" + std::to_string(rho) + "; -1 < rho < 1 required.");
        }
        _rho = rho;
    }

    void Ellipse::set_sigma_x(double sigma_x) {
        if(!(sigma_x >= 0))
        {
            throw std::invalid_argument("Invalid sigma_x=" + std::to_string(sigma_x) +
                                        "; sigma_x >= 0 required.");
        }
        _sigma_x = sigma_x;
    }

    void Ellipse::set_sigma_y(double sigma_y) {
        if(!(sigma_y >= 0))
        {
            throw std::invalid_argument("Invalid sigma_y=" + std::to_string(sigma_y) +
                                        "; sigma_y >= 0 required.");
        }
        _sigma_y = sigma_y;
    }

    Ellipse::Ellipse(EllipseTerms & terms)
    {
        this->set(terms.sigma_x, terms.sigma_y, terms.rho);
    }

    EllipseMajor::EllipseMajor(double r_major, double axrat, double angle, bool degrees) :
        _r_major(r_major), _axrat(axrat),_angle(angle), _degrees(degrees)
    {
        EllipseMajor::check(r_major, axrat, angle);
    }

    void EllipseMajor::set_r_major(double r_major) {
        _r_major = r_major;
    }

    void EllipseMajor::set_axrat(double axrat) {
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
}

#endif