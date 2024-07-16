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

#ifndef LSST_GAUSS2D_EVALUATE_H

#include "lsst/gauss2d/evaluate.h"

namespace lsst::gauss2d {
namespace detail {
Terms terms_from_covar(const double weight, const Ellipse& ell) {
    double sig_x = ell.get_sigma_x();
    double sig_y = ell.get_sigma_y();
    double rho = ell.get_rho();

    const double norm_exp = 1. / (1 - rho * rho);
    Terms rval = {.weight = weight / (2. * M_PI * sig_x * sig_y) * sqrt(norm_exp),
                  .xx = norm_exp / (2. * sig_x * sig_x),
                  .yy = norm_exp / (2. * sig_y * sig_y),
                  .xy = rho * norm_exp / (sig_x * sig_y)};
    return rval;
}
}  // namespace detail
}  // namespace lsst::gauss2d
#endif
