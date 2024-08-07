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

#ifndef LSST_GAUSS2D_ELLIPSE_H

#include <cmath>
#include <stdexcept>
#include <string>

#include "lsst/gauss2d/ellipse.h"
#include "lsst/gauss2d/to_string.h"
#include "lsst/gauss2d/type_name.h"

namespace lsst::gauss2d {
template <class S>
std::pair<S, S> sincos(S arg) {
    return {std::sin(arg), std::cos(arg)};
}

Covariance::Covariance(double sigma_x_sq, double sigma_y_sq, double cov_xy) {
    set(sigma_x_sq, sigma_y_sq, cov_xy);
}

Covariance::Covariance(const Ellipse& ell) { set(ell); }

void Covariance::check(double sigma_x_sq, double sigma_y_sq, double cov_xy) {
    double offdiag_max = sqrt(sigma_x_sq) * sqrt(sigma_y_sq);
    // Define implied rho as -1 for negative cov and +1 for positive cov if zero size
    // This enforces cov_xy == 0 if sigma_x_sq == sigma_y_sq == 0
    double rho = offdiag_max > 0 ? cov_xy / offdiag_max : (cov_xy > 0) - (cov_xy < 0);
    if (!(sigma_x_sq >= 0) || !(sigma_y_sq >= 0) || !(rho >= -1 && rho <= 1)) {
        throw std::invalid_argument("Invalid sigma_x_sq, sigma_y_sq, cov_xy=" + to_string_float(sigma_x_sq)
                                    + "," + to_string_float(sigma_y_sq) + "," + to_string_float(cov_xy)
                                    + " with implied rho=" + to_string_float(rho)
                                    + "; sigma_x,y_sq >= 0 and -1 < rho < 1 required.");
    }
}

void Covariance::convolve(const Covariance& cov) {
    double sigma_x_sq = this->get_sigma_x_sq() + cov.get_sigma_x_sq();
    double sigma_y_sq = this->get_sigma_y_sq() + cov.get_sigma_y_sq();
    double cov_xy = this->get_cov_xy() + cov.get_cov_xy();
    this->set(sigma_x_sq, sigma_y_sq, cov_xy);
}

std::shared_ptr<Covariance> Covariance::make_convolution(const Covariance& cov) const {
    std::shared_ptr<Covariance> cov_ret = std::make_shared<Covariance>(
            this->get_sigma_x_sq(), this->get_sigma_y_sq(), this->get_cov_xy());
    cov_ret->convolve(cov);
    return cov_ret;
}

void Covariance::set(const Ellipse& ellipse) {
    double sigma_x = ellipse.get_sigma_x();
    double sigma_y = ellipse.get_sigma_y();
    double rho = ellipse.get_rho();
    this->_sigma_x_sq = sigma_x * sigma_x;
    this->_sigma_y_sq = sigma_y * sigma_y;
    // Must be done like so in case sigma_x/y^2 round to zero, and check after
    this->set_cov_xy(sigma_x * sigma_y * (_sigma_x_sq > 0) * (_sigma_y_sq > 0) * rho);
}

void Covariance::set(double sigma_x_sq, double sigma_y_sq, double cov_xy) {
    set_sigma_x_sq(sigma_x_sq);
    set_sigma_y_sq(sigma_y_sq);
    set_cov_xy(cov_xy);
}

void Covariance::set_sigma_x_sq(double sigma_x_sq) {
    if (!(sigma_x_sq >= 0)) {
        throw std::invalid_argument(this->str() + " can't set invalid sigma_x_sq="
                                    + to_string_float(sigma_x_sq) + "; sigma_x_sq >= 0 required.");
    }
    _sigma_x_sq = sigma_x_sq;
}

void Covariance::set_sigma_y_sq(double sigma_y_sq) {
    if (!(sigma_y_sq >= 0)) {
        throw std::invalid_argument(this->str() + " can't set invalid sigma_y_sq="
                                    + to_string_float(sigma_y_sq) + "; sigma_y_sq >= 0 required.");
    }
    _sigma_y_sq = sigma_y_sq;
}

void Covariance::set_cov_xy(double cov_xy) {
    // Take individual sqrt just to be safe and avoid potential overflow
    double offdiag_max = sqrt(_sigma_x_sq) * sqrt(_sigma_y_sq);
    // If offdiag_max is zero, we can only accept cov_xy exactly zero
    double rho = offdiag_max > 0 ? cov_xy / offdiag_max : -!(cov_xy >= 0) + !(cov_xy <= 0);
    if (!(rho > -1 && rho < 1)) {
        throw std::invalid_argument(this->str() + "can't set invalid cov_xy=" + to_string_float(cov_xy)
                                    + " (>0=" + to_string_float(cov_xy > 0)
                                    + ", <0=" + to_string_float(cov_xy < 0)
                                    + ", offdiag_max=" + to_string_float(offdiag_max) + " with implied rho="
                                    + to_string_float(rho) + "; -1 < rho < 1 required.");
    }
    _cov_xy = cov_xy;
}

void Covariance::set_xyc(const std::array<double, 3>& xyc) { this->set(xyc[0], xyc[1], xyc[2]); }

std::string Covariance::repr(bool name_keywords, std::string_view namespace_separator) const {
    return type_name_str<Covariance>(false, namespace_separator) + "(" + (name_keywords ? "sigma_x_sq=" : "")
           + to_string_float(_sigma_x_sq) + ", " + (name_keywords ? "sigma_y_sq=" : "")
           + to_string_float(_sigma_y_sq) + ", " + (name_keywords ? "cov_xy=" : "") + to_string_float(_cov_xy)
           + ")";
}

std::string Covariance::str() const {
    return type_name_str<Covariance>(true) + "(sigma_x_sq=" + to_string_float(_sigma_x_sq)
           + ", sigma_y_sq=" + to_string_float(_sigma_y_sq) + ", cov_xy=" + to_string_float(_cov_xy) + ")";
}

bool Covariance::operator==(const Covariance& other) const { return get_xyc() == other.get_xyc(); }
bool Covariance::operator!=(const Covariance& other) const { return !(*this == other); }

std::ostream& operator<<(std::ostream& out, const Covariance& obj) {
    out << obj.str();
    return out;
}

double EllipseData::get_area() const {
    double rho = this->get_rho();
    return M_PI * this->get_sigma_xy() * sqrt(1 - rho * rho);
}

double EllipseData::get_sigma_xy() const { return this->get_sigma_x() * this->get_sigma_y(); }

void EllipseData::convolve(const Ellipse& ell) {
    double sigma_x_ell = ell.get_sigma_x();
    double sigma_y_ell = ell.get_sigma_y();
    double sigma_x = this->get_sigma_x();
    sigma_x = sqrt(sigma_x * sigma_x + sigma_x_ell * sigma_x_ell);
    double sigma_y = this->get_sigma_y();
    sigma_y = sqrt(sigma_y * sigma_y + sigma_y_ell * sigma_y_ell);
    double rho = (this->get_cov_xy() + ell.get_cov_xy()) / (sigma_x * sigma_y);
    this->set(sigma_x, sigma_y, rho);
}

double EllipseData::get_cov_xy() const { return this->get_sigma_x() * this->get_sigma_y() * this->get_rho(); }

double EllipseData::get_hwhm_x() const { return M_SIGMA_HWHM * this->get_sigma_x(); };

double EllipseData::get_hwhm_y() const { return M_SIGMA_HWHM * this->get_sigma_y(); };

std::array<double, 3> EllipseData::get_hxyr() const {
    return {this->get_hwhm_x(), this->get_hwhm_y(), this->get_rho()};
}
std::array<double, 3> EllipseData::get_xyr() const {
    return {this->get_sigma_x(), this->get_sigma_y(), this->get_rho()};
}

double EllipseData::get_radius_trace() const { return sqrt(this->get_sigma_x_sq() + this->get_sigma_y_sq()); }

double EllipseData::get_sigma_x_sq() const {
    double sigma = this->get_sigma_x();
    return sigma * sigma;
}

double EllipseData::get_sigma_y_sq() const {
    double sigma = this->get_sigma_y();
    return sigma * sigma;
}

double EllipseValues::get_sigma_x() const { return *_sigma_x; }
double EllipseValues::get_sigma_y() const { return *_sigma_y; }
double EllipseValues::get_rho() const { return *_rho; }
std::array<double, 3> EllipseValues::get_xyr() const { return {*_sigma_x, *_sigma_y, *_rho}; }

void EllipseValues::set_sigma_x(double sigma_x) {
    EllipseData::check_size(sigma_x, "(size=sigma_x)");
    *_sigma_x = sigma_x;
}

void EllipseValues::set_sigma_y(double sigma_y) {
    EllipseData::check_size(sigma_y, "(size=sigma_y)");
    *_sigma_y = sigma_y;
}

void EllipseValues::set_rho(double rho) {
    EllipseData::check_rho(rho);
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
    *_sigma_x = M_HWHM_SIGMA * hwhm_x;
    *_sigma_y = M_HWHM_SIGMA * hwhm_y;
    *_rho = rho;
}

std::string EllipseValues::repr(bool name_keywords, std::string_view namespace_separator) const {
    return type_name_str<EllipseValues>(false, namespace_separator) + "(" + (name_keywords ? "sigma_x=" : "")
           + to_string_float(this->get_sigma_x()) + ", " + (name_keywords ? "sigma_y=" : "")
           + to_string_float(this->get_sigma_y()) + ", " + (name_keywords ? "rho=" : "")
           + to_string_float(this->get_rho()) + ")";
}

std::string EllipseValues::str() const {
    return type_name_str<EllipseValues>(true) + "(sigma_x=" + to_string_float(this->get_sigma_x())
           + ", sigma_y=" + to_string_float(this->get_sigma_y()) + ", rho=" + to_string_float(this->get_rho())
           + ")";
}

Ellipse::Ellipse(std::shared_ptr<EllipseData> data)
        : _data(data == nullptr ? std::make_shared<EllipseValues>() : std::move(data)) {};

Ellipse::Ellipse(double sigma_x, double sigma_y, double rho) : _data(std::make_shared<EllipseValues>()) {
    _data->set(sigma_x, sigma_y, rho);
}

Ellipse::Ellipse(const Covariance& covar) : _data(std::make_shared<EllipseValues>()) { this->set(covar); }

Ellipse::Ellipse(const EllipseMajor& ellipse) : _data(std::make_shared<EllipseValues>()) {
    this->set(ellipse);
}

std::shared_ptr<Ellipse> Ellipse::make_convolution(const Ellipse& ell) const {
    return this->make_convolution_uniq(ell);
}
std::unique_ptr<Ellipse> Ellipse::make_convolution_uniq(const Ellipse& ell) const {
    // TODO: Replace with cloning derived data
    std::unique_ptr<Ellipse> ell_ret
            = std::make_unique<Ellipse>(this->get_sigma_x(), this->get_sigma_y(), this->get_rho());
    ell_ret->convolve(ell);
    return ell_ret;
}

void EllipseData::set(double sigma_x, double sigma_y, double rho) {
    this->check(sigma_x, sigma_y, rho);
    this->set_sigma_x(sigma_x);
    this->set_sigma_y(sigma_y);
    this->set_rho(rho);
}

void EllipseData::set(const Covariance& covar) {
    double sigma_x = covar.get_sigma_x_sq();
    double sigma_y = covar.get_sigma_y_sq();
    sigma_x = sqrt(sigma_x);
    sigma_y = sqrt(sigma_y);
    double rho = (sigma_x == 0 || sigma_y == 0) ? 0 : covar.get_cov_xy() / (sigma_x * sigma_y);
    this->set(sigma_x, sigma_y, rho);
}

void EllipseData::set(const EllipseMajor& ellipse) {
    const double r_major = ellipse.get_r_major();
    if (r_major == 0) {
        this->set(0, 0, 0);
        return;
    }
    const double axrat = ellipse.get_axrat();
    if (axrat == 1) {
        this->set(r_major, r_major, 0);
        return;
    }
    const auto [sin_th, cos_th] = sincos(ellipse.get_angle_radians());
    const double sin_th_sq = sin_th * sin_th;
    const double cos_th_sq = cos_th * cos_th;

    const double r_major_sq = r_major * r_major;
    const double r_minor_sq = r_major_sq * axrat * axrat;
    const double sigma_x = sqrt(cos_th_sq * r_major_sq + sin_th_sq * r_minor_sq);
    const double sigma_y = sqrt(sin_th_sq * r_major_sq + cos_th_sq * r_minor_sq);
    double rho = (sigma_x == 0 || sigma_y == 0)
                         ? 0
                         : sin_th * cos_th * (r_major_sq - r_minor_sq) / (sigma_x * sigma_y);
    this->set(sigma_x, sigma_y, rho);
}

void EllipseData::set_h(double hwhm_x, double hwhm_y, double rho) {
    EllipseData::check(hwhm_x, hwhm_y, rho);
    this->set_hwhm_x(hwhm_x);
    this->set_hwhm_y(hwhm_y);
    this->set_rho(rho);
}
void EllipseData::set_hwhm_x(double hwhm_x) { this->set_sigma_x(M_HWHM_SIGMA * hwhm_x); }
void EllipseData::set_hwhm_y(double hwhm_y) { this->set_sigma_y(M_HWHM_SIGMA * hwhm_y); }
void EllipseData::set_hxyr(const std::array<double, 3>& hxyr) { this->set_h(hxyr[0], hxyr[1], hxyr[2]); }
void EllipseData::set_xyr(const std::array<double, 3>& xyr) { this->set(xyr[0], xyr[1], xyr[2]); }

std::string Ellipse::repr(bool name_keywords, std::string_view namespace_separator) const {
    return type_name_str<Ellipse>(false, namespace_separator) + "(" + (name_keywords ? "data=" : "")
           + _data->repr(name_keywords, namespace_separator) + ")";
}

std::string Ellipse::str() const { return type_name_str<Ellipse>(true) + "(data=" + _data->str() + ")"; }

bool Ellipse::operator==(const Ellipse& other) const { return get_data() == other.get_data(); };
bool Ellipse::operator!=(const Ellipse& other) const { return !(*this == other); }

const EllipseData& Ellipse::get_data() const { return *this->_data; }
double Ellipse::get_rho() const { return this->_data->get_rho(); }
double Ellipse::get_sigma_x() const { return this->_data->get_sigma_x(); }
double Ellipse::get_sigma_y() const { return this->_data->get_sigma_y(); }

void Ellipse::set_rho(double rho) { this->_data->set_rho(rho); }
void Ellipse::set_sigma_x(double sigma_x) { this->_data->set_sigma_x(sigma_x); }
void Ellipse::set_sigma_y(double sigma_y) { this->_data->set_sigma_y(sigma_y); }

std::pair<double, double> get_x_pm(double sigma_x_sq, double sigma_y_sq, double cov_xy) {
    double apc = sigma_x_sq + sigma_y_sq;
    double x = apc / 2;
    // TODO: Write tests for this with inputs returning close to zero
    // Probably most efficient but unstable, e.g.:
    // sigma_x, sigma_y, rho = 1.58113883008419, 1.5811388300841895, 0
    // ... yields -3.552713678800501e-15 inside the sqrt
    // double pm = sqrt(apc*apc - 4*(sigma_x_sq*sigma_y_sq - cov_xy*cov_xy))/2;
    // Two more multiplications, but seemingly more stable
    // Cancels out cross term from apc (2*sigma_x_sq*sigma_y_sq)
    double pm = (sigma_x_sq * sigma_x_sq + sigma_y_sq * sigma_y_sq
                 - 2 * (sigma_x_sq * sigma_y_sq - 2 * cov_xy * cov_xy));
    // Return zero if the result is negative
    // TODO: Consider checking if < -machine_eps?
    pm = (pm > 0) ? sqrt(pm) / 2 : 0;
    return {x, pm};
}

void init(EllipseMajor& ellipse, const Covariance& covar, bool degrees) {
    double sigma_x_sq = covar.get_sigma_x_sq();
    double sigma_y_sq = covar.get_sigma_y_sq();
    if (sigma_x_sq == 0 && sigma_y_sq == 0) return;
    double cov_xy = covar.get_cov_xy();
    auto [x, pm] = get_x_pm(sigma_x_sq, sigma_y_sq, cov_xy);
    double r_major = x + pm;
    if (r_major == 0) return;

    double axrat = sqrt((x - pm) / r_major);
    r_major = sqrt(r_major);
    double ang = atan2(2 * cov_xy, sigma_x_sq - sigma_y_sq) / 2;
    if (degrees) ang *= M_180_PI;
    ellipse.set(r_major, axrat, ang);
}

EllipseMajor::EllipseMajor(double r_major, double axrat, double angle, bool degrees)
        : _r_major(r_major), _axrat(axrat), _angle(angle), _degrees(degrees) {
    EllipseMajor::check(r_major, axrat, angle);
}

EllipseMajor::EllipseMajor(const Covariance& covar, bool degrees) : _degrees(degrees) {
    init(*this, covar, degrees);
}

EllipseMajor::EllipseMajor(const Ellipse& ellipse, bool degrees) : _degrees(degrees) {
    init(*this, Covariance(ellipse), degrees);
}

void EllipseMajor::set(double r_major, double axrat, double angle) {
    check(r_major, axrat, angle);
    _r_major = r_major;
    _axrat = axrat;
    _angle = angle;
}

void EllipseMajor::set_r_major(double r_major) {
    if (!(r_major >= 0)) {
        throw std::invalid_argument("Invalid r_major=" + to_string_float(r_major)
                                    + "; r_major >= 0 required.");
    }
    _r_major = r_major;
}

void EllipseMajor::set_axrat(double axrat) {
    if (!(axrat >= 0 && axrat <= 1)) {
        throw std::invalid_argument("Invalid axrat=" + to_string_float(axrat)
                                    + "; 1 >= axrat >= 0 required.");
    }
    _axrat = axrat;
}

void EllipseMajor::set_angle(double angle) { _angle = angle; }

void EllipseMajor::set_degrees(bool degrees) {
    if (degrees && !_degrees) {
        _degrees = true;
        _angle *= M_180_PI;
    } else if (!degrees && _degrees) {
        _degrees = false;
        _angle *= M_PI_180;
    }
}

void EllipseMajor::set_rqa(const std::array<double, 3>& rqa) { this->set(rqa[0], rqa[1], rqa[2]); }

std::string EllipseMajor::repr(bool name_keywords, std::string_view namespace_separator) const {
    return type_name_str<EllipseMajor>(false, namespace_separator) + "(" + (name_keywords ? "r_major=" : "")
           + to_string_float(_r_major) + ", " + (name_keywords ? "axrat=" : "") + to_string_float(_axrat)
           + ", " + (name_keywords ? "angle=" : "") + to_string_float(_angle) + ", "
           + (name_keywords ? "degrees=" : "") + to_string_float(_degrees) + ")";
}

std::string EllipseMajor::str() const {
    return type_name_str<EllipseMajor>(true) + "(r_major=" + to_string_float(_r_major)
           + ", axrat=" + to_string_float(_axrat) + ", angle=" + to_string_float(_angle)
           + ", degrees=" + to_string_float(_degrees) + ")";
}

bool EllipseMajor::operator==(const EllipseMajor& other) const {
    return (get_r_major() == other.get_r_major()) && (get_axrat() == other.get_axrat())
           && (get_angle_degrees() == other.get_angle_degrees());
}

bool EllipseMajor::operator!=(const EllipseMajor& other) const { return !(*this == other); }

}  // namespace lsst::gauss2d

#endif
