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

import pytest

import lsst.gauss2d as g2d
import numpy as np
import sys

prefix_namespace = "lsst.gauss2d."

check_str_repr = "linux" in sys.platform


def test_Image():
    img = g2d.ImageD(1, 1)
    img.fill(1)
    assert img.get_value(0, 0) == 1
    img += 1
    assert img.get_value(0, 0) == 2
    img.set_value_unchecked(0, 0, -1)
    assert img.get_value_unchecked(0, 0) == -1
    str_coordsys = str(img.coordsys)
    if check_str_repr:
        assert str(img) == f"ImageD(coordsys={str_coordsys}, n_rows=1, n_cols=1)"
        assert repr(img) == (
            f"{prefix_namespace}python.ImageD(coordsys={prefix_namespace}{str_coordsys}, n_rows=1, n_cols=1)"
        )
    n_rows, n_cols = 11, 13
    coordsys = g2d.CoordinateSystem(1.3, -0.7, -11.6, 365.1)
    img = g2d.ImageF(n_rows=n_rows, n_cols=n_cols, coordsys=coordsys)
    if check_str_repr:
        assert str(img) == f"ImageF(coordsys={coordsys}, n_rows={n_rows}, n_cols={n_cols})"


def test_numpy_Image():
    data = np.arange(12).astype(np.int32).reshape((3, 4))
    img = g2d.ImageI(data)

    data[0, 0] = -1
    assert img.get_value(0, 0) == -1
    assert img.get_value_unchecked(0, 2) == 2
    assert data[0, 2] == 2

    assert (img.n_rows, img.n_cols) == data.shape
    assert tuple(img.shape) == data.shape
    assert img.size == data.size

    with pytest.raises(IndexError):
        img.get_value(3, 0)
    with pytest.raises(IndexError):
        img.get_value(0, 4)
    with pytest.raises(TypeError):
        img.get_value(0, -1)
    with pytest.raises(TypeError):
        img.get_value(-1, 0)


def test_ImageArray():
    imgs = [g2d.ImageD(3, 2), g2d.ImageD(3, 2)]
    img_arr = g2d.ImageArrayD(imgs)
    assert img_arr.at(0) is imgs[0]
    assert img_arr.size == len(imgs)
    assert len(img_arr) == len(imgs)
    assert img_arr[0] is imgs[0]

    with pytest.raises(ValueError):
        g2d.ImageArrayD([imgs[0], g2d.ImageD(2, 3)])

    with pytest.raises(ValueError):
        g2d.ImageArrayD([imgs[0], None])

    if check_str_repr:
        assert str(img_arr) == f"ImageArrayD(data=[{str(imgs[0])}, {str(imgs[1])}])"
        assert repr(img_arr) == (
            f"{prefix_namespace}ImageArrayD(data=[{repr(imgs[0])}, {repr(imgs[1])}])"
        )


@pytest.fixture(scope="module")
def gaussians():
    gauss1 = g2d.Gaussian(centroid=None, ellipse=g2d.Ellipse(3, 6, 0))
    gauss2 = g2d.Gaussian(centroid=None, ellipse=g2d.Ellipse(4, 8, 0))
    return gauss1, gauss2


@pytest.fixture(scope="module")
def convolved_gaussian(gaussians):
    gauss_conv = g2d.ConvolvedGaussian(gaussians[0], gaussians[1])
    return gauss_conv


@pytest.fixture(scope="module")
def convolved_gaussians(convolved_gaussian):
    gaussians_conv = g2d.ConvolvedGaussians([
        convolved_gaussian, g2d.ConvolvedGaussian(convolved_gaussian.kernel, convolved_gaussian.source)
    ])
    return gaussians_conv


def test_GaussianEvaluator(convolved_gaussians):
    output = g2d.ImageD(coordsys=g2d.CoordinateSystem(x_min=-5, y_min=-5), n_rows=11, n_cols=9)
    evaluator = g2d.GaussianEvaluatorD(gaussians=convolved_gaussians, output=output)
    result = evaluator.loglike_pixel()

    assert result == 0
    assert np.sum(output.data) > 0

    sigma_inv = g2d.ImageD(coordsys=output.coordsys, data=np.ones_like(output.data))
    evaluator_ll = g2d.GaussianEvaluatorD(gaussians=convolved_gaussians, data=output, sigma_inv=sigma_inv)

    result = evaluator_ll.loglike_pixel()
    assert result == 0

    output += 1e-6
    assert evaluator_ll.loglike_pixel() < 0

    if check_str_repr:
        assert str(evaluator_ll) == (
            f"GaussianEvaluatorD(gaussians={str(convolved_gaussians)}, "
            f"do_extra=0, do_output=0, do_residual=0, has_background=0, "
            f"is_sigma_image=1, backgroundtype=0, get_likelihood=1, "
            f"data={str(output)}, sigma_inv={str(sigma_inv)}, output=None, residual=None, "
            f"grads=None, grad_param_map=None, grad_param_factor=None, "
            f"extra_param_map=None, extra_param_factor=None, grad_extra=None, grad_param_idx=[], "
            f"n_cols={output.n_cols}, n_rows={output.n_rows}, coordsys={str(output.coordsys)})"
        )
        assert repr(evaluator_ll) == (
            f"{prefix_namespace}GaussianEvaluatorD(gaussians={repr(convolved_gaussians)}, "
            f"data={repr(output)}, sigma_inv={repr(sigma_inv)}, output=None, residual=None, "
            f"grads=None, grad_param_map=None, grad_param_factor=None, "
            f"extra_param_map=None, extra_param_factor=None, background=None)"
        )


def test_make_gaussians_pixel(convolved_gaussians):
    output = g2d.make_gaussians_pixel_D(convolved_gaussians, n_rows=5, n_cols=7)
    assert np.sum(output.data) > 0
