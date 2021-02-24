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
from _gauss2d import loglike_gaussians_pixel as loglike_gaussians_pixel_pb


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


zero_double = np.zeros(0)
zeros_double = np.zeros((0, 0))
zeros_uint64 = np.zeros((0, 0), dtype=np.uint64)


def loglike_gaussians_pixel(
        data, sigma_inv, gaussians, x_min=None, x_max=None, y_min=None, y_max=None,
        to_add=False, output=None, residual=None, grad=None, grad_param_map=None, grad_param_factor=None,
        sersic_param_map=None, sersic_param_factor=None, background=None,
):
    """A wrapper for the pybind11 loglike_gaussians_pixel function with
    reasonable default values for allnof the arguments, since most of them
    don't accept None (nullptr).

    Parameters
    ----------
    data : `numpy.ndarray`
        Image data if computing likelihood. Must be 2D.
    sigma_inv : `numpy.ndarray`
        Inverse sigma (sqrt variance) image if computing likelihood. Must have
        the same shape as `data`.
    gaussians : `numpy.ndarray`
        Parameters for N Gaussians: (cenx, ceny, flux, sigma_x, sigma_y, rho);
        optionally add PSF terms (sigma_x, sigma_y, rho).
        Shape must be (6|9, N).
    x_min : float
        The x-coordinate of the left edge of the image.
    x_max : float
        The x-coordinate of the right edge of the image.
    y_min : float
        The y-coordinate of the bottom edge of the image.
    y_max : float
        The y-coordinate of the top edge of the image.
    to_add : bool
        If outputting the model, should it add it to the output image or
        overwrite it?
    output : `np.ndarray`
        The output image for the model. Shape must be identical to data
        if computing likelihood.
    residual : `np.ndarray`
        The output image for the residual. Shape must be identical to data
        if computing likelihood.
    grad : `np.ndarray`
        Array to store likelihood gradient or Jacobian. Shape must match
        gaussian.shape to output gradient, or
        (data.shape[0], data.shape[1], N=gaussians.size) to output Jacobian.
    grad_param_map : `np.ndarray`
        Optional array of integer indices of where each Gaussian parameter
        gradient/Jacobian should be placed into the `grad` array.
        Length must be N to match `gaussians`.
    grad_param_factor : `np.ndarray`
        Multiplicative factor for each parameter computed for `grad`
        (e.g. for multiplying by transformation derivatives).
        Shape must match `gaussians`.
    sersic_param_map : `np.ndarray`
        Indices as in `grad_param_map` but for the Sersic index of
        multi-Gaussian (mixture) profiles.
    sersic_param_factor : `np.ndarray`
        Multiplicative factor for each parameter as in `grad_param_factor`,
        but for the Sersic index of multi-Gaussian (mixture) profiles.
    background : `np.ndarray`
        A constant background level; default zero. Shape must be (1,).

    Returns
    -------
    loglike : `float`
        The log-likelihood of the model, or zero if computing the Jacobian
        (the residuals are output into the array `residual`).
    """
    if output is None:
        output = zeros_double
    if residual is None:
        residual = zeros_double
    if grad is None:
        grad = zeros_double
    if grad_param_map is None:
        grad_param_map = zeros_uint64
    if sersic_param_map is None:
        sersic_param_map = zeros_uint64
    if grad_param_factor is None:
        grad_param_factor = zeros_double
    if sersic_param_factor is None:
        sersic_param_factor = zeros_double
    if x_min is None:
        x_min = 0
    if x_max is None:
        x_max = data.shape[1]
    if y_min is None:
        y_min = 0
    if y_max is None:
        y_max = data.shape[0]
    if background is None:
        background = zero_double
    return loglike_gaussians_pixel_pb(
        data=data, sigma_inv=sigma_inv, gaussians=gaussians, x_min=x_min, x_max=x_max,
        y_min=y_min, y_max=y_max, to_add=to_add, output=output, residual=residual, grad=grad,
        grad_param_map=grad_param_map, grad_param_factor=grad_param_factor,
        sersic_param_map=sersic_param_map, sersic_param_factor=sersic_param_factor, background=background)
