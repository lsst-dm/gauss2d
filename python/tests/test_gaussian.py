import pytest

import lsst.gauss2d as g2d

prefix_namespace = "lsst.gauss2d."


def test_Gaussian():
    gauss0 = g2d.Gaussian()
    assert gauss0 == g2d.Gaussian()
    centroid = g2d.Centroid()
    ellipse = g2d.Ellipse()
    gauss1 = g2d.Gaussian(centroid, ellipse)
    gauss2 = g2d.Gaussian(centroid, ellipse, g2d.GaussianIntegralValue(2))
    gauss2.integral_value = gauss1.integral_value
    assert gauss1 == gauss2

    centroid.x = 1
    assert gauss1.centroid is centroid
    kwargs = {
        "str_values_cen": "CentroidValues(x=1.000000e+00, y=0.000000e+00)",
        "str_values_ell": "EllipseValues(sigma_x=0.000000e+00, sigma_y=0.000000e+00, rho=0.000000e+00)",
        "str_integral": "GaussianIntegralValue(value=1.000000e+00)",
    }
    str_gauss_format = (
        "{prefix_namespace}Gaussian(centroid={prefix_namespace}Centroid(data={prefix_namespace}"
        "{str_values_cen}), ellipse={prefix_namespace}Ellipse(data={prefix_namespace}{str_values_ell}), "
        "integral={prefix_namespace}{str_integral})"
    )
    assert str(gauss1) == str_gauss_format.format(prefix_namespace="", **kwargs)
    assert repr(gauss1) == str_gauss_format.format(prefix_namespace=prefix_namespace, **kwargs)


def test_Gaussians():
    gauss1 = g2d.Gaussian()
    gs = g2d.Gaussians([gauss1, gauss1])
    assert gs.at(0) == gs.at(1)
    with pytest.raises(IndexError):
        gs.at(2)
    for idx in range(len(gs)):
        assert gs.at(idx) == gauss1

    str_g1 = str(gauss1)
    str_gs = f"Gaussians(data=[{str_g1}, {str_g1}])"
    assert str(gs) == str_gs


@pytest.fixture(scope="module")
def gaussians():
    gauss1 = g2d.Gaussian(centroid=None, ellipse=g2d.Ellipse(3, 6, 0))
    gauss2 = g2d.Gaussian(centroid=None, ellipse=g2d.Ellipse(4, 8, 0))
    return gauss1, gauss2


@pytest.fixture(scope="module")
def convolved_gaussian(gaussians):
    gauss_conv = g2d.ConvolvedGaussian(gaussians[0], gaussians[1])
    return gauss_conv


def test_ConvolvedGaussian(convolved_gaussian, gaussians):
    assert convolved_gaussian.source is gaussians[0]
    assert convolved_gaussian.kernel is gaussians[1]


def test_ConvolvedGaussians(convolved_gaussian, gaussians):
    conv_gaussian2 = g2d.ConvolvedGaussian(gaussians[1], gaussians[0])
    conv_list = [convolved_gaussian, conv_gaussian2]
    n_conv = len(conv_list)
    conv_gaussians = g2d.ConvolvedGaussians(conv_list)
    assert len(conv_gaussians) == n_conv
    assert conv_gaussians.size == n_conv
    assert conv_gaussians[0] is conv_list [0]
    assert conv_gaussians.at(1) is conv_list[1]
    with pytest.raises(IndexError):
        conv_gaussians[n_conv]
