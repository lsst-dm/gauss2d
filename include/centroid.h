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

#include <memory>
#include <string>

namespace gauss2d {

class CentroidData
{
public:
    virtual double get_x() const = 0;
    virtual std::array<double, 2> get_xy() const = 0;
    virtual double get_y() const = 0;

    virtual void set_x(double x) = 0;
    virtual void set_xy(const std::array<double, 2> & xy) = 0;
    virtual void set_y(double y) = 0;

    virtual std::string str() const = 0;
    virtual ~CentroidData() = default;
};

class CentroidValues : public virtual CentroidData
{
private:
    std::shared_ptr<double> _x;
    std::shared_ptr<double> _y;

public:
    double get_x() const { return *_x; }
    std::array<double, 2> get_xy() const { return {*_x, *_y}; }
    double get_y() const { return *_y; }

    void set_x(double x) { *_x = x; }
    void set_xy(const std::array<double, 2> & xy) { *_x = xy[0]; *_y = xy[1]; };
    void set_y(double y) { *_y = y; }

    std::string str() const {
        return "Centroid(x=" + std::to_string(*_x) + ", y=" + std::to_string(*_y) + ")";
    }

    CentroidValues(std::shared_ptr<double> x, std::shared_ptr<double> y) :
        _x(x == nullptr ? std::make_shared<double>(0) : std::move(x)),
        _y(y == nullptr ? std::make_shared<double>(0) : std::move(y)) {};
    CentroidValues(double x=0, double y=0) :
        _x(std::make_shared<double>(x)), _y(std::make_shared<double>(y)) {};

    virtual ~CentroidValues() {};
};

class Centroid
{
private:
    std::shared_ptr<CentroidData> _data;

public:
    double get_x() const { return _data->get_x(); }
    std::array<double, 2> get_xy() const { return _data->get_xy(); }
    double get_y() const { return _data->get_y(); }

    void set_x(double x) { _data->set_x(x); }
    void set_xy(const std::array<double, 2> & xy) { _data->set_xy(xy); }
    void set_y(double y) { _data->set_y(y); }

    std::string str() const {
        return "Centroid(" + _data->str() + ")";
    }

    Centroid(std::shared_ptr<CentroidData> data) : _data(data == nullptr ? std::make_shared<CentroidValues>() : std::move(data)) {}
    Centroid(double x=0, double y=0) : _data(std::make_shared<CentroidValues>()) {
        set_x(x);
        set_y(y);
    }
};

} // namespace gauss2d
#endif
