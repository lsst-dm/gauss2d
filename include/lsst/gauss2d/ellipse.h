// -*- LSST-C++ -*-
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
#define LSST_GAUSS2D_ELLIPSE_H

#include <array>
#include <cmath>
#include <memory>
#include <stdexcept>

#include "object.h"
#include "to_string.h"

namespace lsst::gauss2d {

const double M_HWHM_SIGMA = 1.1774100225154746910115693264599;
const double M_SIGMA_HWHM = 0.84932180028801904272150283410295;
const double M_PI_180 = M_PI / 180.;
const double M_180_PI = 180. / M_PI;

class Ellipse;
class EllipseMajor;

/**
 * @brief A representation of a 2D Gaussian with x and y
 *  standard deviations and a covariance value.
 *
 * This ellipse representation is intended for intermediate
 * calculations and does not use an abstract Data class.
 *
 * @note The covariance is equal to sigma_x*sigma_y*rho,
 * where -1 < rho < 1. Covariance values implying rho values outside
 * this range are invalid.
 */
class Covariance : public Object {
public:
    /**
     * @brief Construct a new Covariance object
     *
     * @param sigma_x_sq The value of sigma_x^2
     * @param sigma_y_sq The value of sigma_y^2
     * @param cov_xy The value of the covariance
     */
    explicit Covariance(double sigma_x_sq = 0, double sigma_y_sq = 0, double cov_xy = 0);
    /// Construct a covariance using values from an ellipse instance.
    explicit Covariance(const Ellipse& ell);

    /// Check whether the supplied values are valid, throwing if not.
    static void check(double sigma_x_sq, double sigma_y_sq, double cov_xy);
    /// Convolve with another covariance, adding the values of each parameter to this.
    void convolve(const Covariance& cov);
    /// Get the square of sigma_x
    double get_sigma_x_sq() const { return _sigma_x_sq; };
    /// Get the square of sigma_y
    double get_sigma_y_sq() const { return _sigma_y_sq; };
    /// Get the covariance
    double get_cov_xy() const { return _cov_xy; };
    /// Get the array of sigma_x^2, sigma_y^2, covariance
    std::array<double, 3> get_xyc() const { return {_sigma_x_sq, _sigma_y_sq, _cov_xy}; }

    /**
     * @brief Return the convolution of this with another covariance.
     *
     * Convolution simply sums the values of each respective parameter.
     *
     * @param cov The covariance to convolve with.
     *
     * @return A new covariance with values set to the convolution.
     */
    std::shared_ptr<Covariance> make_convolution(const Covariance& cov) const;
    /// Same as make_convolution(), but returning a unique_ptr
    std::unique_ptr<Covariance> make_convolution_uniq(const Covariance& cov) const;

    /// Set values from an ellipse instance
    void set(const Ellipse& ellipse);
    /// Set all values at once
    void set(double sigma_x_sq = 0, double sigma_y_sq = 0, double cov_xy = 0);
    /// Set the square of sigma_x
    void set_sigma_x_sq(double sigma_x_sq);
    /// Set the square of sigma_y
    void set_sigma_y_sq(double sigma_y_sq);
    /// Set the off-diagonal (covariance) term
    void set_cov_xy(double cov_xy);
    /// Set sigma_x_sq, sigma_y_sq and cov_xy from an array ref
    void set_xyc(const std::array<double, 3>& xyc);

    std::string repr(bool name_keywords = false,
                     std::string_view namespace_separator = Object::CC_NAMESPACE_SEPARATOR) const override;
    std::string str() const override;

    bool operator==(const Covariance& other) const;
    bool operator!=(const Covariance& other) const;

    friend std::ostream& operator<<(std::ostream& out, const Covariance& obj);

private:
    double _sigma_x_sq = 0;
    double _sigma_y_sq = 0;
    double _cov_xy = 0;
};

/**
 * @brief Interface for an object storing Ellipse data.
 *
 * This is an abstract class designed to store data for an Ellipse to
 * retrieve. No additional restrictions are placed on implementations.
 * Some default implementations are provided.
 *
 **/
class EllipseData : public Object {
public:
    virtual ~EllipseData() = default;

    /// Check whether Ellipse parameter values are valid, throwing if not.
    static void check(double size_x, double size_y, double rho, std::string_view error_suffix = "") {
        if (!(size_x >= 0) || !(size_y >= 0) || !(rho >= -1 && rho <= 1)) {
            throw std::invalid_argument("Invalid size_x, size_y, rho=" + to_string_float(size_x) + ","
                                        + to_string_float(size_y) + "," + to_string_float(rho)
                                        + "; sigma_x,y >= 0 and 1 >= rho >= -1 required."
                                        + std::string(error_suffix));
        }
    }
    /// Check whether an x- or x-yaxis size value is valid
    static void check_size(double size, std::string_view error_suffix = "") {
        if (!(size >= 0)) {
            throw std::invalid_argument("Invalid size=" + to_string_float(size) + "; size >= 0 required."
                                        + std::string(error_suffix));
        }
    }
    /// Check whether a rho value is valid
    static void check_rho(double rho, std::string_view error_suffix = "") {
        if (!(rho >= -1 && rho <= 1)) {
            throw std::invalid_argument("Invalid rho=" + to_string_float(rho) + "; 1 >= rho >= -1 required."
                                        + std::string(error_suffix));
        }
    }
    /// Convolve this ellipse with another
    virtual void convolve(const Ellipse& ell);
    /// Return the area of this ellipse, equal to pi*sigma_major*sigma_minor
    virtual double get_area() const;
    /// Return the covariance, equal to sigma_x*sigma_y*rho
    virtual double get_cov_xy() const;
    /// Get the x-axis half-width at half-maximum
    virtual double get_hwhm_x() const;
    /// Get the y-axis half-width at half-maximum
    virtual double get_hwhm_y() const;
    /// Return the trace radius, equal to sqrt(sigma_x^2 + sigma_y^2)
    virtual double get_radius_trace() const;
    /// Get sigma_x
    virtual double get_sigma_x() const = 0;
    /// Get sigma_y
    virtual double get_sigma_y() const = 0;
    /// Get rho
    virtual double get_rho() const = 0;
    /// Get the square of sigma_x
    virtual double get_sigma_x_sq() const;
    /// Get the square of sigma_y
    virtual double get_sigma_y_sq() const;
    /// Return sigma_x*sigma_y
    virtual double get_sigma_xy() const;
    /// Get hwhm_x, hwhm_y, rho
    virtual std::array<double, 3> get_hxyr() const;
    /// Get sigma_x, sigma_y, rho
    virtual std::array<double, 3> get_xyr() const;

    /// Set sigma_x, sigma_y, rho
    virtual void set(double sigma_x, double sigma_y, double rho);
    /// Set values from a Covariance object
    virtual void set(const Covariance& covar);
    /// Set values from an EllipseMajor object
    virtual void set(const EllipseMajor& ellipse);
    /// Set hwhm_x, hwhm_y, rho (half-width at half-max)
    virtual void set_h(double hwhm_x, double hwhm_y, double rho);
    /// Set the x-axis half-width at half-max (FWHM/2)
    virtual void set_hwhm_x(double hwhm_x);
    /// Set the y-axis half-width at half-max (FWHM/2)
    virtual void set_hwhm_y(double hwhm_y);
    /// Set the x-axis dispersion (sigma)
    virtual void set_sigma_x(double sigma_x) = 0;
    /// Set the y-axis dispersion (sigma)
    virtual void set_sigma_y(double sigma_y) = 0;
    /// Set the correlation parameter (rho)
    virtual void set_rho(double rho) = 0;
    /// Set hwhm_x, hwhm_y, rho from an array
    virtual void set_hxyr(const std::array<double, 3>& hxyr);
    /// Set sigma_x, sigma_y, rho from an array
    virtual void set_xyr(const std::array<double, 3>& xyr);

    virtual std::string repr(bool name_keywords = false, std::string_view namespace_separator
                                                         = Object::CC_NAMESPACE_SEPARATOR) const override
            = 0;
    virtual std::string str() const override = 0;

    bool operator==(const EllipseData& other) const { return get_xyr() == other.get_xyr(); };
    bool operator!=(const EllipseData& other) const { return !(*this == other); };

    friend std::ostream& operator<<(std::ostream& out, const EllipseData& obj) {
        out << obj.str();
        return out;
    }
};

/**
 * @brief An EllipseData storing sigma_x, sigma_y, rho values as shared_ptrs.
 *
 * This implementation stores values in shared_ptrs, allowing sharing
 * of values with other instances.
 *
 **/
class EllipseValues : public EllipseData {
public:
    /**
     * @brief Construct a new EllipseValues object
     *
     * @param sigma_x The sigma_x shared_ptr
     * @param sigma_y The sigma_y shared_ptr
     * @param rho The rho shared_ptr, defaulting to a new zero-valued double
     */
    explicit EllipseValues(std::shared_ptr<double> sigma_x, std::shared_ptr<double> sigma_y,
                           std::shared_ptr<double> rho = nullptr)
            : _sigma_x(sigma_x == nullptr ? std::make_shared<double>(0) : std::move(sigma_x)),
              _sigma_y(sigma_y == nullptr ? std::make_shared<double>(0) : std::move(sigma_y)),
              _rho(rho == nullptr ? std::make_shared<double>(0) : std::move(rho)) {};
    /// Construct a new EllipseValues object, creating new pointers for every value
    explicit EllipseValues(double sigma_x = 0, double sigma_y = 0, double rho = 0)
            : _sigma_x(std::make_shared<double>(sigma_x)),
              _sigma_y(std::make_shared<double>(sigma_y)),
              _rho(std::make_shared<double>(rho)) {};

    double get_sigma_x() const override;
    double get_sigma_y() const override;
    double get_rho() const override;
    std::array<double, 3> get_xyr() const override;

    void set(double sigma_x, double sigma_y, double rho) override;
    void set_h(double hwhm_x, double hwhm_y, double rho) override;
    void set_sigma_x(double sigma_x) override;
    void set_sigma_y(double sigma_y) override;
    void set_rho(double rho) override;

    std::string repr(bool name_keywords = false,
                     std::string_view namespace_separator = Object::CC_NAMESPACE_SEPARATOR) const override;
    std::string str() const override;

private:
    std::shared_ptr<double> _sigma_x;
    std::shared_ptr<double> _sigma_y;
    std::shared_ptr<double> _rho;
};

/**
 * @brief An Ellipse with sigma_x, sigma_y, and rho values.
 *
 * This ellipse parameterization requires sigma_x and sigma_y >=0,
 * and the correlation -1 < rho < 1.
 *
 * @note This ellipse parameterization can be represented as a matrix,
 * with sigma_x and sigma_y on the diagonal, and the other two values
 * set to rho.
 */
class Ellipse : public EllipseData {
public:
    explicit Ellipse(std::shared_ptr<EllipseData> data);
    /**
     * @brief Construct a new Ellipse object with default float values.
     *
     * @param sigma_x The initial sigma_x value (default 0)
     * @param sigma_y The initial sigma_y value (default 0)
     * @param rho The intial rho value (default 0)
     */
    explicit Ellipse(double sigma_x = 0, double sigma_y = 0, double rho = 0);
    explicit Ellipse(const Covariance& covar);
    explicit Ellipse(const EllipseMajor& ellipse);

    /// Return a const ref to this ellipse's data
    const EllipseData& get_data() const;
    double get_rho() const override;
    double get_sigma_x() const override;
    double get_sigma_y() const override;
    /**
     * @brief Return the convolution of this with another ellipse.
     *
     * Convolution simply sums the values of the covariance matrix terms.
     *
     * @param ell The ellipse to convolve with.
     *
     * @return A new ellipse with values set to the convolution.
     */
    std::shared_ptr<Ellipse> make_convolution(const Ellipse& ell) const;
    /// Same as make_convolution(), but returning a unique_ptr
    std::unique_ptr<Ellipse> make_convolution_uniq(const Ellipse& ell) const;

    void set_rho(double rho) override;
    void set_sigma_x(double sigma_x) override;
    void set_sigma_y(double sigma_y) override;

    std::string repr(bool name_keywords = false,
                     std::string_view namespace_separator = Object::CC_NAMESPACE_SEPARATOR) const override;
    std::string str() const override;

    bool operator==(const Ellipse& other) const;
    bool operator!=(const Ellipse& other) const;

private:
    std::shared_ptr<EllipseData> _data;
};

/**
 * @brief An Ellipse with r_major, axrat and angle values.
 *
 * This ellipse representation is intended for intermediate
 * calculations and does not use an abstract Data class.
 * r_major must be >= 0, and 0 <= axrat <= 1.
 */
class EllipseMajor : public Object {
public:
    /**
     * @brief Construct a new EllipseMajor object with default float values.
     *
     * @param r_major The initial r_major value
     * @param axrat The initial axrat value
     * @param angle The intial angle value
     * @param degrees Whether the angle unit is degrees (true) or radians (false)
     */
    explicit EllipseMajor(double r_major, double axrat, double angle, bool degrees = false);
    explicit EllipseMajor(const Covariance& covar, bool degrees = false);
    explicit EllipseMajor(const Ellipse& ellipse, bool degrees = false);

    /// Check whether the supplied values are valid, throwing if not.
    static void check(double r_major, double axrat, double angle) {
        if (!(r_major >= 0) || !(axrat >= 0 && axrat <= 1)) {
            throw std::invalid_argument("Invalid r_major, axrat, angle=" + std::to_string(r_major) + ","
                                        + std::to_string(axrat) + "," + std::to_string(angle)
                                        + "; r_major >= 0, 1 >= axrat >= 0 required.");
        }
    }
    /// Return the area of this ellipse, equal to pi*sigma_major*sigma_minor
    double get_area() const { return M_PI * _r_major * _r_major * _axrat; }
    /// Get the major axis length
    double get_r_major() const { return _r_major; }
    /// Get the axis ratio
    double get_axrat() const { return _axrat; }
    /// Get the position angle in the configured units
    double get_angle() const { return _angle; }
    /// Get the position angle in degrees
    double get_angle_degrees() const { return _degrees ? _angle : _angle * M_180_PI; }
    /// Get the position angle in radians
    double get_angle_radians() const { return _degrees ? _angle * M_PI_180 : _angle; }
    /// Get the array of r_major, axrat, angle
    std::array<double, 3> get_rqa() const { return {_r_major, _axrat, _angle}; }
    /// Return if the units are degrees (true) or radians (false)
    bool is_degrees() const { return _degrees; }

    void set(double r_major, double axrat, double angle);
    void set_r_major(double r_major);
    void set_axrat(double axrat);
    void set_angle(double angle);
    void set_degrees(bool degrees);
    void set_rqa(const std::array<double, 3>& rqa);

    std::string repr(bool name_keywords = false,
                     std::string_view namespace_separator = Object::CC_NAMESPACE_SEPARATOR) const override;
    std::string str() const override;

    bool operator==(const EllipseMajor& other) const;
    bool operator!=(const EllipseMajor& other) const;

private:
    double _r_major = 0.;
    double _axrat = 1.;
    double _angle = 0.;
    bool _degrees = false;
};

}  // namespace lsst::gauss2d
#endif
