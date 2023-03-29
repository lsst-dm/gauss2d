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

#ifndef GAUSS2D_CENTROID_H
#define GAUSS2D_CENTROID_H

#include <array>
#include <memory>
#include <string>

#include "object.h"

namespace gauss2d {

/**
 * @brief Interface for an object storing Centroid data.
 * 
 * This is an abstract class designed to store data for a Centroid to
 * retrieve. No additional restrictions are placed on implementations.
 *
**/
class CentroidData : public Object
{
public:
    /// Get the x value
    virtual double get_x() const = 0;
    /// Get the x and y values
    virtual std::array<double, 2> get_xy() const = 0;
    /// Get the y value
    virtual double get_y() const = 0;

    virtual void set_x(double x) = 0;
    virtual void set_xy(const std::array<double, 2> & xy) = 0;
    virtual void set_y(double y) = 0;

    bool operator==(const CentroidData& other) const {
        return (get_x() == other.get_x()) && (get_y() == other.get_y());
    };

    virtual std::string repr(bool name_keywords=false) const override = 0;
    virtual std::string str() const override = 0;
    virtual ~CentroidData() = default;
};

/**
 * @brief A CentroidData storing centroid values as shared_ptrs.
 * 
 * This implementation stores values in shared_ptrs, allowing sharing
 * of values with other instances.
 *
**/
class CentroidValues : public virtual CentroidData
{
private:
    std::shared_ptr<double> _x;
    std::shared_ptr<double> _y;

public:
    double get_x() const override;
    std::array<double, 2> get_xy() const override;
    double get_y() const override;

    void set_x(double x) override;
    void set_xy(const std::array<double, 2> & xy) override;
    void set_y(double y) override;

    std::string repr(bool name_keywords=false) const override;
    std::string str() const override;

    /**
     * @brief Construct a new Centroid Values object
     * 
     * @param x The x-axis centroid value
     * @param y The y-axis centroid value
     */
    CentroidValues(std::shared_ptr<double> x, std::shared_ptr<double> y);
    CentroidValues(double x=0, double y=0);

    virtual ~CentroidValues() {};
};

/**
 * @brief A 2D coordinate representing the center of a plane figure.
 * 
 * This is a centroid in a 2D coordinate system, generally used for ellipses
 * in this package. Storage is implemented in CentroidData.
 *
 * @param data The centroid data
**/
class Centroid : public Object
{
private:
    std::shared_ptr<CentroidData> _data;

public:
    /**
     * @brief Convolve this with another centroid.
     *
     * Convolution simply adds the value of the other centroid to this.
     * 
     * @param cen The centroid to convolve with
     */
    void convolve(const Centroid & cen);

    /// Get this centroid's data
    const CentroidData & get_data() const;
    /// Get the x value
    double get_x() const;
    /// Get the x and y values
    std::array<double, 2> get_xy() const;
    /// Get the y value
    double get_y() const;

    /**
     * @brief Return the convolution of this with another centroid.
     *
     * Convolution simply adds the values of both centroids together
     * 
     * @param cen The centroid to convolve with
     *
     * @return A new centroid with values set to the convolution.
     */
    std::shared_ptr<Centroid> make_convolution(const Centroid& cen) const;
    /// Same as make_convolution(), but returning a unique_ptr.
    std::unique_ptr<Centroid> make_convolution_uniq(const Centroid& cen) const;

    void set_x(double x);
    void set_xy(const std::array<double, 2> & xy);
    void set_y(double y);

    std::string repr(bool name_keywords=false) const override;
    std::string str() const override;

    bool operator==(const Centroid& other) const;

    Centroid(std::shared_ptr<CentroidData> data);
    Centroid(double x=0, double y=0);
};

} // namespace gauss2d
#endif
