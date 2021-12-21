# This file is part of gauss2d.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import numpy as np


def covar_to_ellipse(sigma_x_sq, sigma_y_sq, cov_xy, degrees=False):
    """Convert covariance matrix terms to ellipse major axis, axis ratio and
    position angle representation.

    Parameters
    ----------
    sigma_x_sq, sigma_y_sq : `float` or array-like
        x- and y-axis squared standard deviations of a 2-dimensional normal
        distribution (diagonal terms of its covariance matrix).
        Must be scalar or identical length array-likes.
    cov_xy : `float` or array-like
        x-y covariance of a of a 2-dimensional normal distribution
        (off-diagonal term of its covariance matrix).
        Must be scalar or identical length array-likes.
    degrees : `bool`
        Whether to return the position angle in degrees instead of radians.

    Returns
    -------
    r_major, axrat, angle : `float` or array-like
        Converted major-axis radius, axis ratio and position angle
        (counter-clockwise from the +x axis) of the ellipse defined by
        each set of input covariance matrix terms.

    Notes
    -----
    The eigenvalues from the determinant of a covariance matrix are:
    |a-m b|
    |b c-m|
    det = (a-m)(c-m) - b^2 = ac - (a+c)m + m^2 - b^2 = m^2 - (a+c)m + (ac-b^2)
    Solving:
    m = ((a+c) +/- sqrt((a+c)^2 - 4(ac-b^2)))/2
    ...or equivalently:
    m = ((a+c) +/- sqrt((a-c)^2 + 4b^2))/2

    Unfortunately, the latter simplification is not as well-behaved
    in floating point math, leading to square roots of negative numbers when
    one of a or c is very close to zero.

    The values from this function should match those from
    `Ellipse.make_ellipse_major` to within rounding error, except in the
    special case of sigma_x == sigma_y == 0, which returns a NaN axis ratio
    here by default. This function mainly intended to be more convenient
    (and possibly faster) for array-like inputs.
    """
    apc = sigma_x_sq + sigma_y_sq
    x = apc/2
    pm = np.sqrt(apc**2 - 4*(sigma_x_sq*sigma_y_sq - cov_xy**2))/2

    r_major = x + pm
    axrat = np.sqrt((x - pm)/r_major)
    r_major = np.sqrt(r_major)
    angle = np.arctan2(2*cov_xy, sigma_x_sq - sigma_y_sq)/2
    return r_major, axrat, (np.degrees(angle) if degrees else angle)


def gauss2dint(xdivsigma):
    """Return the fraction of the total surface integral of a 2D Gaussian
    contained with a given multiple of its dispersion.

    Parameters
    ----------
    xdivsigma : `float`
        The multiple of the dispersion to integrate to.

    Returns
    -------
    frac : `float`
        The fraction of the surface integral contained within an ellipse of
        size `xdivsigma`.

    Notes
    -----
    This solution can be computed as follows:

    https://www.wolframalpha.com/input/?i=
    Integrate+2*pi*x*exp(-x%5E2%2F(2*s%5E2))%2F(s*sqrt(2*pi))+dx+from+0+to+r
    """
    return 1 - np.exp(-xdivsigma**2/2.)
