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
 * Interface for an object storing Centroid data.
 *
**/
class CentroidData : public Object
{
public:
    virtual double get_x() const = 0;
    virtual std::array<double, 2> get_xy() const = 0;
    virtual double get_y() const = 0;

    virtual void set_x(double x) = 0;
    virtual void set_xy(const std::array<double, 2> & xy) = 0;
    virtual void set_y(double y) = 0;

    bool operator==(const CentroidData& other) const {
        return (get_x() == other.get_x()) && (get_y() == other.get_y());
    };

    virtual std::string str() const override = 0;
    virtual ~CentroidData() = default;
};

/**
 * A CentroidData storing centroid values as shared_ptrs.
 * 
 * shared_ptr usage allows CentroidData to share parameters if desired. 
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

    std::string str() const override;

    CentroidValues(std::shared_ptr<double> x, std::shared_ptr<double> y);
    CentroidValues(double x=0, double y=0);

    virtual ~CentroidValues() {};
};

/**
 * A Centroid is a 2D coordinate representing the center of a plane figure
 * (specifically an ellipse in this package).
 *
 * The storage of the parameters is implemented in CentroidData;
 * this class serves as a storage-independent interface for usage in
 * Ellipse classes.
 *
**/
class Centroid : public Object
{
private:
    std::shared_ptr<CentroidData> _data;

public:
    void convolve(const Centroid & cen);

    const CentroidData & get_data() const;
    double get_x() const;
    std::array<double, 2> get_xy() const;
    double get_y() const;

    std::shared_ptr<Centroid> make_convolution(const Centroid& ell) const;
    std::unique_ptr<Centroid> make_convolution_uniq(const Centroid & ell) const;

    void set_x(double x);
    void set_xy(const std::array<double, 2> & xy);
    void set_y(double y);

    std::string str() const override;

    bool operator==(const Centroid& other) const;

    Centroid(std::shared_ptr<CentroidData> data);
    Centroid(double x=0, double y=0);
};

} // namespace gauss2d
#endif
