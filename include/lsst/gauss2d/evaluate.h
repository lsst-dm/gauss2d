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

#ifndef LSST_GAUSS2D_EVALUATE_H
#define LSST_GAUSS2D_EVALUATE_H

#include <iostream>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>

#include <vector>

#include "gaussian.h"
#include "image.h"
#include "to_string.h"

namespace lsst::gauss2d {

static const ConvolvedGaussians GAUSSIANS_NULL{std::nullopt};

enum class BackgroundType : unsigned char {
    none = 0,
    constant = 1,
};
enum class OutputType : unsigned char {
    none = 0,
    overwrite = 1,
    add = 2,
};
enum class GradientType : unsigned char {
    none = 0,
    loglike = 1,
    jacobian = 2,
};

const size_t N_EXTRA_MAP = 2;
const size_t N_EXTRA_FACTOR = 3;
const size_t N_PARAMS_GAUSS2D = 6;

namespace detail {

typedef std::vector<double> vecd;
typedef std::unique_ptr<vecd> vecdptr;

struct ValuesGauss {
    double cen_x;
    double cen_y;
    double L;
    double sigma_x;
    double sigma_y;
    double rho;
};

struct TermsMoments {
    vecdptr x;
    vecdptr x_norm;
    vecdptr xx;
};

inline TermsMoments gaussian_pixel_x_xx(const double cen_x, const double x_min, const double bin_x,
                                        const unsigned int dim_x, const double xx_weight,
                                        const double xy_weight) {
    const double x_init = x_min - cen_x + bin_x / 2.;
    vecdptr x = std::make_unique<vecd>(dim_x);
    vecdptr xx = std::make_unique<vecd>(dim_x);
    vecdptr x_norm = std::make_unique<vecd>(dim_x);
    for (unsigned int i = 0; i < dim_x; i++) {
        double dist = x_init + i * bin_x;
        (*x)[i] = dist;
        (*xx)[i] = dist * dist * xx_weight;
        (*x_norm)[i] = dist * xy_weight;
    }
    return TermsMoments({std::move(x), std::move(x_norm), std::move(xx)});
}

/*
This is some largely unnecessary algebra to derive conversions between the ellipse parameterization and the
covariance matrix parameterization of a bivariate Gaussian.

See e.g. http://mathworld.wolfram.com/BivariateNormalDistribution.html
... and https://www.unige.ch/sciences/astro/files/5413/8971/4090/2_Segransan_StatClassUnige.pdf

tan 2th = 2*rho*sig_x*sig_y/(sig_x^2 - sig_y^2)

sigma_maj^2 = (cos2t*sig_x^2 - sin2t*sig_y^2)/(cos2t-sin2t)
sigma_min^2 = (cos2t*sig_y^2 - sin2t*sig_x^2)/(cos2t-sin2t)

(sigma_maj^2*(cos2t-sin2t) + sin2t*sig_y^2)/cos2t = sig_x^2
-(sigma_min^2*(cos2t-sin2t) - cos2t*sig_y^2)/sin2t = sig_x^2

(sigma_maj^2*(cos2t-sin2t) + sin2t*sig_y^2)/cos2t = -(sigma_min^2(cos2t-sin2t) - cos2t*sig_y^2)/sin2t
sin2t*(sigma_maj^2*(cos2t-sin2t) + sin2t*sig_y^2)c + cos2t*(sigma_min^2*(cos2t-sin2t) - cos2t*sig_y^2) = 0
cos4t*sig_y^2 - sin^4th*sig_y^2 = sin2t*(sigma_maj^2*(cos2t-sin2t)) + cos2t*(sigma_min^2*(cos2t-sin2t))

cos^4x - sin^4x = (cos^2 + sin^2)*(cos^2-sin^2) = cos^2 - sin^2

sig_y^2 = (sin2t*(sigma_maj^2*(cos2t-sin2t)) + cos2t*(sigma_min^2*(cos2t-sin2t)))/(cos2t - sin2t)
       = (sin2t*sigma_maj^2 + cos2t*sigma_min^2)

sigma_maj^2*(cos2t-sin2t) + sin2t*sig_y^2 = cos2t*sig_x^2
sig_x^2 = (sigma_maj^2*(cos2t-sin2t) + sin2t*sig_y^2)/cos2t
       = (sigma_maj^2*(cos2t-sin2t) + sin2t*(sin2t*sigma_maj^2 + cos2t*sigma_min^2))/cos2t
       = (sigma_maj^2*(cos2t-sin2t+sin4t) + sin2tcos2t*sigma_min^2)/cos2t
       = (sigma_maj^2*cos4t + sin2t*cos2t*sigma_min^2)/cos2t
       = (sigma_maj^2*cos2t + sigma_min^2*sin2t)

sig_x^2 - sig_y^2 = sigma_maj^2*cos2t + sigma_min^2*sin2t - sigma_maj^2*sin2t - sigma_min^2*cos2t
                = (sigma_maj^2 - sigma_min^2*)*(cos2t-sin2t)
                = (sigma_maj^2 - sigma_min^2*)*(1-tan2t)*cos2t

rho = tan2th/2/(sig_x*sig_y)*(sig_x^2 - sig_y^2)
    = tanth/(1-tan2t)/(sig_x*sig_y)*(sig_x^2 - sig_y^2)
    = tanth/(1-tan2t)/(sig_x*sig_y)*(sigma_maj^2 - sigma_min^2*)*(1-tan2t)*cos2t
    = tanth/(sig_x*sig_y)*(sigma_maj^2 - sigma_min^2)*cos2t
    = sint*cost/(sig_x*sig_y)*(sigma_maj^2 - sigma_min^2)

Derivation of 2D Gaussian function gradients (Jacobian):

    xmc = x - cen_x

    m = L*bin_x*bin_y/(2.*M_PI*sig_x*sig_y)/sqrt(1-rho*rho) * exp(-(
        + xmc[g]^2/sig_x^2/(2*(1-rho*rho))
        + ymc[g][j]^2/sig_y^2/(2*(1-rho*rho))
        - rho*xmc[g]*ymc[g][j]/(1-rho*rho)/sig_x/sig_y
    ))

    norm_exp = 1./(2*(1-rho*rho))
    weight = L*bin_x*bin_y
    T = Terms
    T.weight = weight/(2.*M_PI*sig_x*sig_y)*sqrt(2.*norm_exp),
    T.xx = norm_exp/sig_x/sig_x,
    T.yy = norm_exp/sig_y/sig_y,
    T.xy = 2.*rho*norm_exp/sig_x/sig_y

    m = T.weight * exp(-(xmc[g]^2*T.xx + ymc[g][j]^2*T.yy - rho*xmc[g]*ymc[g][j]*T.xy))

    rho -> r; sig_x/y -> s;

    weight_xx[g] = T.xx
    dxmc^2/dcen_x = -2*xmc[g] - r*ymc[g][j]*T.xy

    ymc_weighted = r*ymc*T.xy

    dm/dL = m/L
    dm/dcen_x = (2*xmc[g]*weight_xx[g] - ymc_weighted[g][j])*m
    dm/ds = -m/s + m*(2/s*[xy]sqs[g] - xy/s) = m/s*(2*[xy]sqs[g] - (1+xy))
    dm/dr = -m*(r/(1-r^2) + 4*r*(1-r^2)*(xmc_sq_norm[g] + yy_weighted[g][j]) -
        (1/r + 2*r/(1+r)*xmc[g]*ymc[g][j])

    To verify (all on one line):

    https://www.wolframalpha.com/input/?i=differentiate+F*s*t%2F(2*pi*sqrt(1-r%5E2))*
    exp(-((x-m)%5E2*s%5E2+%2B+(y-n)%5E2*t%5E2+-+2*r*(x-m)*(y-n)*s*t)%2F(2*(1-r%5E2)))+wrt+s
*/

typedef std::array<double, N_PARAMS_GAUSS2D> Weights;

// Various multiplicative terms that appread in a Gaussian PDF
struct Terms {
    double weight;
    double xx;
    double yy;
    double xy;
};

Terms terms_from_covar(double weight, const Ellipse& ell);

/**
 *  Storage for terms common to Gaussians for a single pixel.
 **/
class TermsPixel {
public:
    double weight = 0;
    double weight_kernel = 0;
    double xmc = 0;
    double xmc_weighted = 0;
    vecdptr ymc_weighted = nullptr;
    double weight_xx = 0;
    double xmc_sq_norm = 0;
    vecdptr yy_weighted = nullptr;
    // TODO: Give these variables more compelling names

    explicit TermsPixel(double weight_i = 0, double weight_kernel = 0, double xmc_i = 0,
                        double xmc_weighted_i = 0, vecdptr ymc_weighted_i = nullptr, double weight_xx_i = 0,
                        double xmc_sq_norm_i = 0, vecdptr yy_weighted_i = nullptr)
            : weight(weight_i),
              xmc(xmc_i),
              xmc_weighted(xmc_weighted_i),
              ymc_weighted(std::move(ymc_weighted_i)),
              weight_xx(weight_xx_i),
              xmc_sq_norm(xmc_sq_norm_i),
              yy_weighted(std::move(yy_weighted_i)) {}

    // Set values that are specific to a given gaussian
    void set(double weight_, double weight_kernel_, double xmc_, double xx, vecdptr ymc_weighted_,
             vecdptr yy_weighted_) {
        this->weight = weight_;
        this->weight_kernel = weight_kernel_;
        this->xmc = xmc_;
        this->weight_xx = xx;
        this->ymc_weighted = std::move(ymc_weighted_);
        this->yy_weighted = std::move(yy_weighted_);
    }
};

typedef std::vector<TermsPixel> TermsPixelVec;

class TermsGradient {
public:
    vecdptr ymc = nullptr;
    double xx_weight = 0;
    double xy_weight = 0;
    double yy_weight = 0;
    double rho_factor = 0;
    double sig_x_inv = 0;
    double sig_y_inv = 0;
    double rho_xy_factor = 0;
    double sig_x_src_div_conv = 0;
    double sig_y_src_div_conv = 0;
    double drho_c_dsig_x_src = 0;
    double drho_c_dsig_y_src = 0;
    double drho_c_drho_s = 0;
    double xmc_t_xy_weight = 0;

    explicit TermsGradient(vecdptr ymc_i = nullptr, double xx_weight_i = 0, double xy_weight_i = 0,
                           double yy_weight_i = 0, double rho_factor_i = 0, double sig_x_inv_i = 0,
                           double sig_y_inv_i = 0, double rho_xy_factor_i = 0,
                           double sig_x_src_div_conv_i = 0, double sig_y_src_div_conv_i = 0,
                           double drho_c_dsig_x_src_i = 0, double drho_c_dsig_y_src_i = 0,
                           double drho_c_drho_s_i = 0, double xmc_t_xy_weight_i = 0)
            : ymc(std::move(ymc_i)),
              xx_weight(xx_weight_i),
              xy_weight(xy_weight_i),
              yy_weight(yy_weight_i),
              rho_factor(rho_factor_i),
              sig_x_inv(sig_x_inv_i),
              sig_y_inv(sig_y_inv_i),
              rho_xy_factor(rho_xy_factor_i),
              sig_x_src_div_conv(sig_x_src_div_conv_i),
              sig_y_src_div_conv(sig_y_src_div_conv_i),
              drho_c_dsig_x_src(drho_c_dsig_x_src_i),
              drho_c_dsig_y_src(drho_c_dsig_y_src_i),
              drho_c_drho_s(drho_c_drho_s_i),
              xmc_t_xy_weight(xmc_t_xy_weight_i) {}

    void set(const Terms& terms, const Ellipse& ell_src, const Ellipse& ell_psf, const Ellipse& ell,
             vecdptr ymc_) {
        double sig_x = ell.get_sigma_x();
        double sig_y = ell.get_sigma_y();
        double rho = ell.get_rho();

        this->ymc = std::move(ymc_);
        const double sig_xy = sig_x * sig_y;
        this->xx_weight = terms.xx;
        this->xy_weight = terms.xy;
        this->yy_weight = terms.yy;
        this->rho_factor = rho / (1. - rho * rho);
        this->sig_x_inv = 1 / sig_x;
        this->sig_y_inv = 1 / sig_y;
        this->rho_xy_factor = 1. / (1. - rho * rho) / sig_xy;
        /*
    sigma_conv = sqrt(sigma_src^2 + sigma_psf^2)
    https://www.wolframalpha.com/input/?i=d(sqrt(x%5E2%2By%5E2))%2Fdx
    dsigma_conv/dsigma_src = sigma_src/sigma_conv
    df/dsigma_src = df/d_sigma_conv * dsigma_conv/dsigma_src

    rho_conv*sigmaxy_conv = rho_src*sigmaxy_src + rho_psf*sigmaxy_psf
    rho_conv = (rho_src*sigmaxy_src + rho_psf*sigmaxy_psf)/sigmaxy_conv
    drho_conv/drho_src = sigmaxy_src/sigmaxy_conv

    drho_conv/dsigmax_src = rho_src*sigmay_src/sigmay_conv*sigmax_psf^2/sigmax_conv^3 (*sigmax_src/sigmax_src)
                          = rho_conv/sigmax_src*(sigmax_psf/sigmax_conv)^2
        + rho_psf*sigmax_psf*sigmay_psf*sigmax_src/sigmay_conv/sigmax_conv^3

        */
        double sig_x_src = ell_src.get_sigma_x();
        double sig_y_src = ell_src.get_sigma_y();
        double rho_src = ell_src.get_rho();
        double sig_x_psf = ell_psf.get_sigma_x();
        double sig_y_psf = ell_psf.get_sigma_y();
        double rho_psf = ell_psf.get_rho();

        this->sig_x_src_div_conv = sig_x_src / sig_x;
        this->sig_y_src_div_conv = sig_y_src / sig_y;

        double covar_psf_dsig_xy = rho_psf * sig_x_psf * sig_y_psf / sig_xy;
        double drho = 0;
        if (sig_x_psf > 0) {
            double sig_p_ratio = sig_x_psf / sig_x;
            drho = rho_src * this->sig_y_src_div_conv * sig_p_ratio * sig_p_ratio / sig_x
                   - covar_psf_dsig_xy * this->sig_x_src_div_conv / sig_x;
        }
        this->drho_c_dsig_x_src = drho;
        if (sig_y_psf > 0) {
            double sig_p_ratio = sig_y_psf / sig_y;
            drho = rho_src * this->sig_x_src_div_conv * sig_p_ratio * sig_p_ratio / sig_y
                   - covar_psf_dsig_xy * this->sig_y_src_div_conv / sig_y;
        }
        this->drho_c_dsig_y_src = drho;
        this->drho_c_drho_s = sig_x_src * sig_y_src / sig_xy;
    }
};

typedef std::vector<TermsGradient> TermsGradientVec;

template <typename t, class Data, class Indices>
class GradientsExtra : public Object {
public:
    GradientsExtra(const Image<idx_type, Indices>& param_map_in, const Image<t, Data>& param_factor_in,
                   std::shared_ptr<ImageArray<t, Data>> output, size_t n_gauss)
            : _param_map(param_map_in), _param_factor(param_factor_in), _output(output) {
        if(output == nullptr) {
            throw std::invalid_argument("GradientsExtra output must not be a nullptr");
        }
        const auto n_extra_map_rows = param_map_in.get_n_rows();
        const auto n_extra_fac_rows = param_factor_in.get_n_rows();
        const auto n_extra_map_cols = param_map_in.get_n_cols();
        const auto n_extra_fac_cols = param_factor_in.get_n_cols();
        std::string errmsg = "";
        if (n_extra_map_rows != n_gauss) {
            errmsg += "extra_param_map n_rows=" + std::to_string(n_extra_map_rows)
                      + " != n_gauss=" + std::to_string(n_gauss) + ". ";
        }
        if (n_extra_fac_rows != n_gauss) {
            errmsg += "extra_param_factor n_rows=" + std::to_string(n_extra_fac_rows)
                      + " != n_gauss=" + std::to_string(n_gauss) + ". ";
        }
        if (n_extra_map_cols != N_EXTRA_MAP) {
            errmsg += "extra_param_map n_cols=" + std::to_string(n_extra_map_cols) + " != 2. ";
        }
        if (n_extra_fac_cols != N_EXTRA_FACTOR) {
            errmsg += "extra_param_factor n_cols=" + std::to_string(n_extra_fac_cols) + " != 3. ";
        }
        if (!errmsg.empty()) throw std::runtime_error(errmsg);
    }

    GradientsExtra() = delete;

    inline void add_index_to_set(size_t g, std::set<size_t>& set) { set.insert(_param_map.get_value(g, 1)); }

    inline void add(size_t g, size_t dim1, size_t dim2, const ValuesGauss& gradients) {
        // Reset g_check to zero once g rolls back to zero itself
        _g_check *= (g != 0);
        const bool to_add = g == _param_map.get_value(_g_check, 0);
        const auto idx = _param_map.get_value(g, 1);
        double value = gradients.L * _param_factor.get_value(g, 0)
                       + gradients.sigma_x * _param_factor.get_value(g, 1)
                       + gradients.sigma_y * _param_factor.get_value(g, 2);
        (*_output)[idx].add_value_unchecked(dim1, dim2, value);
        _g_check += to_add;
    }

    std::string repr(bool name_keywords = false,
                     std::string_view namespace_separator = Object::CC_NAMESPACE_SEPARATOR) const override {
        bool is_kw = name_keywords;
        std::string rval = (  //
                type_name_str<GradientsExtra<t, Data, Indices>>(false, namespace_separator) + "("
                + (is_kw ? "param_map=" : "") + _param_map.repr(is_kw, namespace_separator) + ", "
                + (is_kw ? "param_factor=" : "") + _param_factor.repr(is_kw, namespace_separator) + ", "
                + (is_kw ? "output=" : "") + _output->repr(is_kw, namespace_separator) + ", "
                + (is_kw ? "n_gauss=" : "") + std::to_string(_param_map.get_n_rows()) + ")"  //
        );
        return rval;
    }
    std::string str() const override {
        std::string rval = (                                                  //
                type_name_str<GradientsExtra<t, Data, Indices>>(true) + "("   //
                + "param_map=" + _param_map.str() + ", "                      //
                + "param_factor=" + _param_factor.str() + ", "                //
                + "output=" + _output->str() + ", "                            //
                + "n_gauss=" + std::to_string(_param_map.get_n_rows()) + ")"  //
        );
        return rval;
    }

private:
    const Image<idx_type, Indices>& _param_map;
    const Image<t, Data>& _param_factor;
    std::shared_ptr<ImageArray<t, Data>> _output;
    size_t _g_check = 0;
};

template <class Data>
inline void gaussian_pixel(Data& image, const double weight, const double cen_x, const double cen_y,
                           const double xx_weight, const double yy_weight, const double xy_weight) {
    const unsigned int dim_x = image.get_n_cols();
    const unsigned int dim_y = image.get_n_rows();
    const auto coordsys = image.get_coordsys();
    const double x_min = coordsys.x_min;
    const double y_min = coordsys.y_min;
    const double bin_x = coordsys.get_dx1() / image.get_n_cols();
    const double bin_y = coordsys.get_dy2() / image.get_n_rows();

    const double weight_pix = weight * bin_x * bin_y;

    const auto moment_terms_y = gaussian_pixel_x_xx(cen_y, y_min, bin_y, dim_y, yy_weight, xy_weight);
    const vecd& y_weighted = *(moment_terms_y.x_norm);
    const vecd& yy_weighted = *(moment_terms_y.xx);

    double x = x_min - cen_x + bin_x / 2.;
    // TODO: Consider a version with cached xy, although it doubles memory usage
    for (unsigned int i = 0; i < dim_x; i++) {
        const double xx_weighted = x * x * xx_weight;
        for (unsigned int j = 0; j < dim_y; j++) {
            image.set_value_unchecked(j, i,
                                      weight_pix * exp(-(xx_weighted + yy_weighted[j] - x * y_weighted[j])));
        }
        x += bin_x;
    }
}

/*
// Evaluate a Gaussian on a grid given the three elements of the symmetric covariance matrix
// Actually, rho is scaled by sig_x and sig_y (i.e. the covariance is rho*sig_x*sig_y)
template <class Data>
void add_gaussian_pixel(const Gaussian & gauss, Data & image)
{
    const double L = gauss.get_integral();
    const Terms terms = terms_from_covar(L, *gauss.ell);

    gaussian_pixel(image, terms.weight, gauss.cen->x, gauss.cen->y, terms.xx, terms.yy, terms.xy);
}

template <class Data>
void add_gaussian_pixel_extra(const Gaussian & gauss, Data & data)
{
    // I don't remember why this isn't just 1/(2*ln(2)) but anyway it isn't
    const double weight_sersic = 0.69314718055994528622676398299518;

    size_t dim_x = data.get_n_cols();
    size_t dim_y = data.get_n_rows();

    const auto coordsys = data.get_coordsys();
    const double bin_x = coordsys.get_dx1();
    const double bin_y = coordsys.get_dy2();

    const auto ell = EllipseMajor(*gauss.ell);
    double r_eff = ell.get_r_major();
    double axrat = ell.get_axrat();
    double ang = ell.get_angle_radians();
    double cen_x = gauss.cen->x;
    double cen_y = gauss.cen->y;

    const double weight = gauss.get_integral()*weight_sersic/(M_PI*axrat)/(r_eff*r_eff)*bin_x*bin_y;
    const double axrat_sq_inv = 1.0/axrat/axrat;

    double x,y;
    const double bin_x_half=bin_x/2.;
    const double bin_y_half=bin_y/2.;

    const double reff_y_inv = sin(ang)*sqrt(weight_sersic)/r_eff;
    const double reff_x_inv = cos(ang)*sqrt(weight_sersic)/r_eff;


    unsigned int i=0,j=0;
    x = coordsys.x_min - cen_x + bin_x_half;
    for(i = 0; i < dim_x; i++)
    {
        y = coordsys.y_min - cen_y + bin_y_half;
        for(j = 0; j < dim_y; j++)
        {
           const double dist_1 = (x*reff_x_inv + y*reff_y_inv);
           const double dist_2 = (x*reff_y_inv - y*reff_x_inv);
           // mat(j,i) = ... is slower, but perhaps allows for Indices with dim_x*dim_y > INT_MAX ?
           data.set_value_unchecked(j, i, weight*exp(-(dist_1*dist_1 + dist_2*dist_2*axrat_sq_inv)));
           y += bin_y;
        }
        x += bin_x;
    }
}
*/
inline void gaussian_pixel_get_jacobian(ValuesGauss& out, const double m, const double m_unweight,
                                        const double xmc_norm, const double ymc_weighted, const double xmc,
                                        const double ymc, const double norms_yy, const double xmc_t_xy_factor,
                                        const double sig_x_inv, const double sig_y_inv, const double xy_norm,
                                        const double xx_norm, const double yy_weighted,
                                        const double rho_div_one_m_rhosq, const double norms_xy_div_rho,
                                        const double dsig_x_conv_src = 1, const double dsig_y_conv_src = 1,
                                        const double drho_c_dsig_x_src = 0,
                                        const double drho_c_dsig_y_src = 0, const double drho_c_drho_s = 1) {
    out.cen_x = m * (2 * xmc_norm - ymc_weighted);
    out.cen_y = m * (2 * ymc * norms_yy - xmc_t_xy_factor);
    out.L = m_unweight;
    const double one_p_xy = 1. + xy_norm;
    // df/dsigma_src = df/d_sigma_conv * dsigma_conv/dsigma_src
    const double dfdrho_c = m
                            * (rho_div_one_m_rhosq * (1 - 2 * (xx_norm + yy_weighted - xy_norm))
                               + xmc * ymc * norms_xy_div_rho);
    out.sigma_x = dfdrho_c * drho_c_dsig_x_src + dsig_x_conv_src * (m * sig_x_inv * (2 * xx_norm - one_p_xy));
    out.sigma_y
            = dfdrho_c * drho_c_dsig_y_src + dsig_y_conv_src * (m * sig_y_inv * (2 * yy_weighted - one_p_xy));
    // drho_conv/drho_src = sigmaxy_src/sigmaxy_conv
    out.rho = dfdrho_c * drho_c_drho_s;
}

inline void gaussian_pixel_get_jacobian_from_terms(ValuesGauss& out, const size_t dim1,
                                                   const TermsPixel& terms_pixel,
                                                   const TermsGradient& terms_grad, double m,
                                                   double m_unweighted, double xy_norm) {
    gaussian_pixel_get_jacobian(
            out, m, m_unweighted * terms_pixel.weight_kernel, terms_pixel.xmc_weighted,
            (*terms_pixel.ymc_weighted)[dim1], terms_pixel.xmc, (*terms_grad.ymc)[dim1], terms_grad.yy_weight,
            terms_grad.xmc_t_xy_weight, terms_grad.sig_x_inv, terms_grad.sig_y_inv, xy_norm,
            terms_pixel.xmc_sq_norm, (*terms_pixel.yy_weighted)[dim1], terms_grad.rho_factor,
            terms_grad.rho_xy_factor, terms_grad.sig_x_src_div_conv, terms_grad.sig_y_src_div_conv,
            terms_grad.drho_c_dsig_x_src, terms_grad.drho_c_dsig_y_src, terms_grad.drho_c_drho_s);
}

template <typename t>
inline void gaussian_pixel_add_values(t& cen_x, t& cen_y, t& L, t& sig_x, t& sig_y, t& rho,
                                      const ValuesGauss& values, const t weight_cen_x = 1,
                                      const t weight_cen_y = 1, const t weight_L = 1,
                                      const t weight_sig_x = 1, const t weight_sig_y = 1,
                                      const t weight_rho = 1) {
    cen_x += weight_cen_x * values.cen_x;
    cen_y += weight_cen_y * values.cen_y;
    L += weight_L * values.L;
    sig_x += weight_sig_x * values.sigma_x;
    sig_y += weight_sig_y * values.sigma_y;
    rho += weight_rho * values.rho;
}

// Computes and stores LL along with dll/dx for all components
template <typename t, class Data, class Indices, bool do_extra>
inline void gaussians_pixel_add_like_grad(
        Image<t, Data> * output, const Image<idx_type, Indices>& grad_param_map,
        const Image<t, Data>& grad_param_factor, const size_t N_GAUSS,
        const std::vector<Weights>& gaussweights, double& chi_pix, double& loglike, const double model,
        double data, double sigma_inv, unsigned int dim1, unsigned int dim2, const TermsPixelVec& terms_pixel,
        const TermsGradientVec& terms_grad, ValuesGauss& gradients,
        const std::shared_ptr<const Image<idx_type, Indices>> extra_param_map,
        const std::shared_ptr<const Image<t, Data>> extra_param_factor) {
    double diff = data - model;
    chi_pix = sigma_inv * diff;
    double diffvar = sigma_inv * chi_pix;
    loglike -= diff * diffvar / 2.;
    for (size_t g = 0; g < N_GAUSS; ++g) {
        const Weights& weights = gaussweights[g];
        /*
        * Derivation:
        *
            ll = sum(-(data-model)^2*varinv/2)
            dll/dx = --2*dmodel/dx*(data-model)*varinv/2
            dmodelsum/dx = d(.. + model[g])/dx = dmodel[g]/dx
            dll/dx = dmodel[g]/dx*diffvar
        */
        gaussian_pixel_get_jacobian_from_terms(gradients, dim1, terms_pixel[g], terms_grad[g],
                                               weights[0] * diffvar, weights[1] * diffvar, weights[2]);
        gaussian_pixel_add_values<t>(
                output->_get_value_unchecked(0, grad_param_map.get_value_unchecked(g, 0)),
                output->_get_value_unchecked(0, grad_param_map.get_value_unchecked(g, 1)),
                output->_get_value_unchecked(0, grad_param_map.get_value_unchecked(g, 2)),
                output->_get_value_unchecked(0, grad_param_map.get_value_unchecked(g, 3)),
                output->_get_value_unchecked(0, grad_param_map.get_value_unchecked(g, 4)),
                output->_get_value_unchecked(0, grad_param_map.get_value_unchecked(g, 5)), gradients,
                grad_param_factor.get_value_unchecked(g, 0), grad_param_factor.get_value_unchecked(g, 1),
                grad_param_factor.get_value_unchecked(g, 2), grad_param_factor.get_value_unchecked(g, 3),
                grad_param_factor.get_value_unchecked(g, 4), grad_param_factor.get_value_unchecked(g, 5));
        if constexpr (do_extra) {
            double value = gradients.L * extra_param_factor->get_value_unchecked(g, 0)
                           + gradients.sigma_x * extra_param_factor->get_value_unchecked(g, 1)
                           + gradients.sigma_y * extra_param_factor->get_value_unchecked(g, 2);
            output->_get_value_unchecked(0, extra_param_map->get_value_unchecked(g, 1)) += value;
        }
    }
}

template <typename T, class Data, class Indices, GradientType gradient_type, bool do_extra>
inline T gaussian_pixel_add_all(size_t g, size_t j, size_t i, double weight, double sigma_inv,
                                const TermsPixelVec& terms_pixel_vec, ImageArray<T, Data>& output_jac_ref,
                                const Image<idx_type, Indices>& grad_param_map,
                                const Image<T, Data>& grad_param_factor, std::vector<Weights>& gradweights,
                                const TermsGradientVec& terms_grad_vec, ValuesGauss& gradients,
                                GradientsExtra<T, Data, Indices>* grad_extra) {
    const TermsPixel& terms_pixel = terms_pixel_vec[g];
    const double xy_norm = terms_pixel.xmc * (*terms_pixel.ymc_weighted)[j];
    const double value_unweight
            = terms_pixel.weight * exp(-(terms_pixel.xmc_sq_norm + (*terms_pixel.yy_weighted)[j] - xy_norm));
    const double value = weight * value_unweight;
    // Computes dmodel/dx for x in [cen_x, cen_y, flux, sigma_x, sigma_y, rho]
    if constexpr (gradient_type == GradientType::loglike) {
        Weights& output = gradweights[g];
        output[0] = value;
        output[1] = value_unweight;
        output[2] = xy_norm;
    } else if constexpr (gradient_type == GradientType::jacobian) {
        gaussian_pixel_get_jacobian_from_terms(gradients, j, terms_pixel, terms_grad_vec[g],
                                               sigma_inv * value, sigma_inv * value_unweight, xy_norm);
        gaussian_pixel_add_values<T>(
                output_jac_ref[grad_param_map.get_value_unchecked(g, 0)]._get_value_unchecked(j, i),
                output_jac_ref[grad_param_map.get_value_unchecked(g, 1)]._get_value_unchecked(j, i),
                output_jac_ref[grad_param_map.get_value_unchecked(g, 2)]._get_value_unchecked(j, i),
                output_jac_ref[grad_param_map.get_value_unchecked(g, 3)]._get_value_unchecked(j, i),
                output_jac_ref[grad_param_map.get_value_unchecked(g, 4)]._get_value_unchecked(j, i),
                output_jac_ref[grad_param_map.get_value_unchecked(g, 5)]._get_value_unchecked(j, i),
                gradients, grad_param_factor.get_value_unchecked(g, 0),
                grad_param_factor.get_value_unchecked(g, 1), grad_param_factor.get_value_unchecked(g, 2),
                grad_param_factor.get_value_unchecked(g, 3), grad_param_factor.get_value_unchecked(g, 4),
                grad_param_factor.get_value_unchecked(g, 5));
        if constexpr (do_extra) grad_extra->add(g, j, i, gradients);
    }
    return value;
}

template <typename Indices>
const std::shared_ptr<const Indices> _param_map_default(size_t n_gaussians,
                                                        size_t n_params = N_PARAMS_GAUSS2D,
                                                        size_t increment = 1) {
    auto param_map = std::make_shared<Indices>(n_gaussians, n_params);
    size_t index = 0;
    for (size_t g = 0; g < n_gaussians; ++g) {
        for (size_t p = 0; p < n_params; ++p) {
            param_map->set_value(g, p, index);
            index += increment;
        }
    }
    return std::const_pointer_cast<const Indices>(param_map);
}

template <typename Data>
const std::shared_ptr<const Data> _param_factor_default(size_t n_gaussians,
                                                        size_t n_params = N_PARAMS_GAUSS2D,
                                                        double value = 1.) {
    auto param_factor = std::make_shared<Data>(n_gaussians, n_params);
    for (size_t g = 0; g < n_gaussians; ++g) {
        for (size_t p = 0; p < n_params; ++p) {
            param_factor->set_value(g, p, value);
        }
    }
    return std::const_pointer_cast<const Data>(param_factor);
}
}  // namespace detail

/**
 * @brief A class that evaluates 2D Gaussians and renders them in images.
 *
 * This class is designed for repeated re-evaluation of Gaussian mixture models
 * by caching references to the inputs and outputs.
 *
 * @tparam T The data type (e.g. float, int)
 * @tparam Data The data array class
 * @tparam Indices The index array class (usually a size_t array)
 */
template <typename T, class Data, class Indices>
class GaussianEvaluator : public Object {
public:
    typedef Image<T, Data> DataT;
    typedef ImageArray<T, Data> ImageArrayT;
    typedef Image<idx_type, Indices> IndicesT;

    GaussianEvaluator(int x = 0, const std::shared_ptr<const ConvolvedGaussians> gaussians = nullptr) {}

    /**
     * @brief Construct a GaussianEvaluator, inferring outputs from inputs.
     *
     * @param gaussians N x 6 matrix of Gaussian parameters [cen_x, cen_y, flux, sigma_x, sigma_y, rho]
     * @param coordsys Coordinate system for all images.
     * @param data 2D input Data matrix.
     * @param sigma_inv 2D inverse sigma (sqrt variance) map of the same size as data.
     * @param output 2D output matrix of the same size as ImageD.
     * @param residual 2D output matrix for residual ((data-model)/sigma) of the same size as data.
     * @param grads Output for gradients. Can either an M x 1 vector or M x Data 3D Jacobian matrix,
     *    where M <= N x 6 to allow for condensing gradients based on grad_param_map.
     * @param grad_param_map Nx6 matrix of indices of grad to add each gradient to. For example, if four
     *    gaussians share the same cen_x, one could set grad_param_map[0:4,0] = 0. All values must have
     *    index < grad.size().
     * @param grad_param_factor Nx6 matrix of multiplicative factors for each gradient term. For example, if a
     *    Gaussian is a sub-component of a multi-Gaussian component with a total flux parameter but fixed
     *    ratios, as in multi-Gaussian Sersic models.
     * @param extra_param_map Nx2 matrix of indices to add to an extra (meta)parameter. The first item is
     *    the index of the Gaussian to add and the second is the index of the metaparameter.
     * @param extra_param_factor Nx3 matrix of multiplicative factors for each extra gradient term. The
     *    factors are ordered L, sigma_x, sigma_y.
     * @param background A background model image. Only 1x1 (constant level) backgrounds are supported.
     */
    GaussianEvaluator(const std::shared_ptr<const ConvolvedGaussians> gaussians,
                      const std::shared_ptr<const Image<T, Data>> data = nullptr,
                      const std::shared_ptr<const Image<T, Data>> sigma_inv = nullptr,
                      const std::shared_ptr<Image<T, Data>> output = nullptr,
                      const std::shared_ptr<Image<T, Data>> residual = nullptr,
                      const std::shared_ptr<ImageArrayT> grads = nullptr,
                      const std::shared_ptr<const Image<idx_type, Indices>> grad_param_map = nullptr,
                      const std::shared_ptr<const Image<T, Data>> grad_param_factor = nullptr,
                      const std::shared_ptr<const Image<idx_type, Indices>> extra_param_map = nullptr,
                      const std::shared_ptr<const Image<T, Data>> extra_param_factor = nullptr,
                      const std::shared_ptr<const Image<T, Data>> background = nullptr)
            : _gaussians(gaussians == nullptr ? GAUSSIANS_NULL : *gaussians),
              _gaussians_ptr(gaussians == nullptr ? nullptr : gaussians),
              _n_gaussians(_gaussians.size()),
              _do_output(output != nullptr),
              _is_sigma_image((sigma_inv != nullptr) && (sigma_inv->size() > 1)),
              _do_residual(residual != nullptr),
              _has_background(background != nullptr),
              _gradienttype((grads != nullptr && grads->size() > 0)
                                    ? (((grads->size() == 1) && ((*grads)[0].get_n_rows() == 1))
                                               ? GradientType::loglike
                                               : GradientType::jacobian)
                                    : GradientType::none),
              _do_extra(
                      extra_param_map != nullptr && extra_param_factor != nullptr
                      && (_gradienttype == GradientType::loglike || _gradienttype == GradientType::jacobian)),
              _backgroundtype(background != nullptr ? (_background->size() == 1 ? BackgroundType::constant
                                                                                : BackgroundType::none)
                                                    : BackgroundType::none),
              _get_likelihood((data != nullptr) && (sigma_inv != nullptr)),
              _data(data),
              _sigma_inv(sigma_inv),
              _output(output),
              _residual(residual),
              _grads(grads),
              _grad_param_map((_gradienttype == GradientType::none)
                                      ? nullptr
                                      : ((grad_param_map != nullptr)
                                                 ? grad_param_map
                                                 : detail::_param_map_default<Indices>(_gaussians.size()))),
              _grad_param_factor(
                      (_gradienttype == GradientType::none)
                              ? nullptr
                              : ((grad_param_factor != nullptr)
                                         ? grad_param_factor
                                         : detail::_param_factor_default<Data>(_gaussians.size()))),
              _extra_param_map((_gradienttype == GradientType::none)
                                       ? nullptr
                                       : ((extra_param_map != nullptr)
                                                  ? extra_param_map
                                                  : detail::_param_map_default<Indices>(_gaussians.size(),
                                                                                        N_EXTRA_MAP, 0))),
              _extra_param_factor((_gradienttype == GradientType::none)
                                          ? nullptr
                                          : ((extra_param_factor != nullptr)
                                                     ? extra_param_factor
                                                     : detail::_param_factor_default<Data>(
                                                               _gaussians.size(), N_EXTRA_FACTOR, 0.))),
              _grad_extra(_do_extra ? std::make_unique<detail::GradientsExtra<T, Data, Indices>>(
                                              *extra_param_map, *extra_param_factor, _grads, _n_gaussians)
                                    : nullptr),
              _grad_param_idx(_grad_param_map == nullptr ? std::vector<size_t>{} : _get_grad_param_idx()),
              _n_cols(_data == nullptr ? (_output == nullptr ? (_gradienttype == GradientType::jacobian
                                                                        ? (*_grads)[0].get_n_cols()
                                                                        : 0)
                                                             : _output->get_n_cols())
                                       : _data->get_n_cols()),
              _n_rows(_data == nullptr ? (_output == nullptr ? (_gradienttype == GradientType::jacobian
                                                                        ? (*_grads)[0].get_n_rows()
                                                                        : 0)
                                                             : _output->get_n_rows())
                                       : _data->get_n_rows()),
              _size(_n_cols * _n_rows),
              _coordsys(_data == nullptr ? (_output == nullptr ? (_gradienttype == GradientType::jacobian
                                                                          ? (*_grads)[0].get_coordsys()
                                                                          : COORDS_DEFAULT)
                                                               : _output->get_coordsys())
                                         : _data->get_coordsys()) {
        if (gaussians == nullptr) throw std::runtime_error("Gaussians can't be null");
        if (_has_background && background->size() > 1)
            throw std::runtime_error("Background model size can't be > 1");
        if (!(_n_cols > 0 && _n_rows > 0))
            throw std::runtime_error("Data can't have n_rows/cols == 0; "
                                     + std::to_string(output == nullptr));
        if (!_do_extra && ((extra_param_map != nullptr) || (extra_param_factor != nullptr))) {
            throw std::runtime_error(
                    "Must pass all of extra_param_map, extra_param_factor and grads"
                    " to compute extra gradients");
        }
        if (_get_likelihood) {
            // The case of constant variance per pixel
            if (_is_sigma_image) {
                if (_n_cols != _sigma_inv->get_n_cols() || _n_rows != _sigma_inv->get_n_rows()) {
                    throw std::runtime_error("Data matrix dimensions [" + std::to_string(_n_cols) + ','
                                             + std::to_string(_n_rows)
                                             + "] don't match inverse variance dimensions ["
                                             + std::to_string(_sigma_inv->get_n_cols()) + ','
                                             + std::to_string(_sigma_inv->get_n_rows()) + ']');
                }
                if (_sigma_inv->get_coordsys() != _data->get_coordsys()) {
                    throw std::runtime_error("sigma_inv coordsys=" + _sigma_inv->get_coordsys().str()
                                             + "!= coordsys=" + get_coordsys().str());
                }
            }
        } else if ((data != nullptr) || (sigma_inv != nullptr)) {
            throw std::runtime_error("Passed only one non-null data/sigma_inv");
        }
        // TODO: Add more coordsys compatibility checks
        if (_gradienttype != GradientType::none) {
            const size_t n_gpm = _grad_param_map->get_n_rows();
            const size_t n_gpf = _grad_param_factor->get_n_rows();

            if ((n_gpm != _n_gaussians) || (n_gpf != _n_gaussians)) {
                throw std::invalid_argument("nrows for grad_param_map,factor=" + std::to_string(n_gpm) + ","
                                            + std::to_string(n_gpf)
                                            + " != n_gaussians=" + std::to_string(_n_gaussians));
            }

            if ((_grad_param_map->get_n_cols() != N_PARAMS_GAUSS2D)
                || (_grad_param_factor->get_n_cols() != N_PARAMS_GAUSS2D)) {
                throw std::invalid_argument("n_cols for grad_param_map,factor=" + std::to_string(n_gpm) + ","
                                            + std::to_string(n_gpf));
            }

            const size_t n_grad = (_gradienttype == GradientType::jacobian) ? _grads->size()
                                                                            : (_grads->at(0).get_n_cols());
            size_t idx_g_max = 0;
            size_t idx_e_gauss_max = 0;
            size_t idx_e_param_max = 0;
            for (size_t g = 0; g < _n_gaussians; ++g) {
                for (size_t p = 0; p < N_PARAMS_GAUSS2D; ++p) {
                    auto value = _grad_param_map->get_value_unchecked(g, p);
                    if (value > idx_g_max) idx_g_max = value;
                }
                auto value = _extra_param_map->get_value_unchecked(g, 0);
                if (value > idx_e_gauss_max) idx_e_gauss_max = value;
                value = _extra_param_map->get_value_unchecked(g, 1);
                if (value > idx_e_param_max) idx_e_param_max = value;
            }
            if (!((idx_e_gauss_max < _n_gaussians) && (idx_e_param_max < n_grad))) {
                throw std::invalid_argument(" max extra_param_map[0]=" + std::to_string(idx_e_gauss_max)
                                            + " !< n_gaussians=" + std::to_string(_n_gaussians) + " and/or "
                                            + " max extra_param_map[1]=" + std::to_string(idx_e_param_max)
                                            + " !< n_grad=" + std::to_string(n_grad));
            }
            if (!(idx_g_max < n_grad)) {
                throw std::invalid_argument("max grad_param_map,extra_param_map[1]="
                                            + std::to_string(idx_g_max)
                                            + " !< n_grad=" + std::to_string(n_grad));
            }
        }
        if (_gradienttype == GradientType::loglike) {
            if (!_get_likelihood) {
                throw std::runtime_error(
                        "Can't compute likelihood gradient without computing likelihood;"
                        " did you pass data and sigma_inv arrays?");
            }
        } else if (_gradienttype == GradientType::jacobian) {
            // This should never happen because the ImageArray constructor requires identical dimensions, but
            // anyway
            size_t idx_jac = 0;
            for (const auto& jac_img : *_grads) {
                if (!(jac_img->get_n_rows() == _n_rows && jac_img->get_n_cols() == _n_cols)) {
                    throw std::runtime_error("grads[0] matrix dimensions [" + std::to_string(_n_cols) + ','
                                             + std::to_string(_n_rows)
                                             + "] don't match Jacobian matrix (grads["
                                             + std::to_string(idx_jac) + "]) dimensions ["
                                             + std::to_string(jac_img->get_n_cols()) + ','
                                             + std::to_string(jac_img->get_n_rows()) + ']');
                }
            }
        }
    }
    ~GaussianEvaluator() override = default;

    Data const& IMAGE_NULL_CONST() const { return this->_data_null_const; }
    Indices const& INDICES_NULL_CONST() const { return this->_indices_null_const; }
    ImageArray<T, Data> const& IMAGEARRAY_NULL_CONST() const { return this->_data_array_null_const; }

    const CoordinateSystem& get_coordsys() const { return _coordsys; }
    size_t get_n_cols() const { return (_n_cols); }
    size_t get_n_rows() const { return (_n_rows); }
    size_t get_size() const { return (_size); }

    /**
     * Compute the model and/or log-likelihood and/or gradient (d(log-likelihood)/dx)
     * and/or Jacobian (dmodel/dx) for a Gaussian mixture model.
     *
     * This function calls a series of templated functions with explicit instantiations. This is solely for
     * the purpose of avoiding having to manually write a series of nested conditionals. My hope is that the
     * templating will insert no-op functions wherever there's nothing to do instead of a needless branch
     * inside each pixel's loop, and that the compiler will actually inline everything for maximum
     * performance. Whether that actually happens or not is anyone's guess.
     *
     * TODO: Consider override to compute LL and Jacobian, even if it's only useful for debug purposes.
     *
     * @return The log likelihood.
     */
    double loglike_pixel(bool to_add = false) {
        const OutputType output_type
                = this->_do_output ? (to_add ? OutputType::add : OutputType::overwrite) : OutputType::none;
        return output_type == OutputType::overwrite
                       ? loglike_gaussians_pixel_output<OutputType::overwrite>()
                       : (output_type == OutputType::add
                                  ? loglike_gaussians_pixel_output<OutputType::add>()
                                  : loglike_gaussians_pixel_output<OutputType::none>());
    }

    std::string repr(bool name_keywords = false,
                     std::string_view namespace_separator = Object::CC_NAMESPACE_SEPARATOR) const override {
        std::string null = namespace_separator == Object::CC_NAMESPACE_SEPARATOR ? "nullptr" : "None";
        // These are sadly just to keep formatting aligned nicely
        auto is_kw = name_keywords;
        auto name_sep = namespace_separator;
        std::string c = ", ";
        std::string rval = (  // I really want this on a separate line; I'd also prefer single indent, but...
                type_name_str<GaussianEvaluator>(false, name_sep) + "("  // make clang-format align nicely
                + (is_kw ? "gaussians=" : "") + repr_ptr(_gaussians_ptr, is_kw, name_sep) + ", "  //
                + (is_kw ? "data=" : "") + repr_ptr(_data, is_kw, name_sep) + ", "                //
                + (is_kw ? "sigma_inv=" : "") + repr_ptr(_sigma_inv, is_kw, name_sep) + ", "      //
                + (is_kw ? "output=" : "") + repr_ptr(_output, is_kw, name_sep) + ", "            //
                + (is_kw ? "residual=" : "") + repr_ptr(_residual, is_kw, name_sep) + ", "        //
                + (is_kw ? "grads=" : "") + repr_ptr(_grads, is_kw, name_sep) + ", "              //
                + (is_kw ? "grad_param_map=" : "") + repr_ptr(_grad_param_map, is_kw, name_sep) + ", "
                + (is_kw ? "grad_param_factor=" : "") + repr_ptr(_grad_param_factor, is_kw, name_sep) + c
                + (is_kw ? "extra_param_map=" : "") + repr_ptr(_extra_param_map, is_kw, name_sep) + ", "
                + (is_kw ? "extra_param_factor=" : "") + repr_ptr(_extra_param_factor, is_kw, name_sep) + c
                + (is_kw ? "background=" : "") + repr_ptr(_background, is_kw, name_sep) + ")");
        return rval;
    }
    std::string str() const override {
        std::string c = ", ";
        // The comments force clang-format to keep to a more readable one var per line
        std::string rval = (                                                                 //
                type_name_str<GaussianEvaluator>(true) + "("                                 //
                + "gaussians=" + str_ptr(_gaussians_ptr) + c                                 //
                + "do_extra=" + std::to_string(_do_extra) + c                                //
                + "do_output=" + std::to_string(_do_output) + c                              //
                + "do_residual=" + std::to_string(_do_residual) + c                          //
                + "has_background=" + std::to_string(_has_background) + c                    //
                + "is_sigma_image=" + std::to_string(_is_sigma_image) + c                    //
                + "backgroundtype=" + std::to_string(static_cast<int>(_backgroundtype)) + c  //
                + "get_likelihood=" + std::to_string(_get_likelihood) + c                    //
                + "data=" + str_ptr(_data) + c                                               //
                + "sigma_inv=" + str_ptr(_sigma_inv) + c                                     //
                + "output=" + str_ptr(_output) + c                                           //
                + "residual=" + str_ptr(_residual) + c                                       //
                + "grads=" + str_ptr(_grads) + c                                             //
                + "grad_param_map=" + str_ptr(_grad_param_map) + c                           //
                + "grad_param_factor=" + str_ptr(_grad_param_factor) + c                     //
                + "extra_param_map=" + str_ptr(_extra_param_map) + c                         //
                + "extra_param_factor=" + str_ptr(_extra_param_factor) + c                   //
                + "grad_extra=" + str_ptr(_grad_extra ? _grad_extra.get() : nullptr) + c     //
                + "grad_param_idx=" + to_string_iter(_grad_param_idx) + c                    //
                + "n_cols=" + std::to_string(_n_cols) + c                                    //
                + "n_rows=" + std::to_string(_n_rows) + c                                    //
                + "coordsys=" + _coordsys.str()                                              //
                + ")"                                                                        //
        );
        return rval;
    }

private:
    // These are private because they originally return static objects that
    // were (and are) not meant to be mutated, although they probably can be
    ImageArrayT& IMAGEARRAY_NULL() const;

    const ConvolvedGaussians& _gaussians;
    const std::shared_ptr<const ConvolvedGaussians> _gaussians_ptr;
    const size_t _n_gaussians;

    const bool _do_output;
    const bool _is_sigma_image;
    const bool _do_residual;
    const bool _has_background;
    const GradientType _gradienttype;
    const bool _do_extra;
    const BackgroundType _backgroundtype;
    const bool _get_likelihood;

    const std::shared_ptr<const DataT> _data;
    const std::shared_ptr<const DataT> _sigma_inv;
    const std::shared_ptr<DataT> _output;
    const std::shared_ptr<DataT> _residual;
    const std::shared_ptr<ImageArrayT> _grads;
    const std::shared_ptr<const IndicesT> _grad_param_map;
    const std::shared_ptr<const DataT> _grad_param_factor;
    const std::shared_ptr<const IndicesT> _extra_param_map;
    const std::shared_ptr<const DataT> _extra_param_factor;
    const std::unique_ptr<detail::GradientsExtra<T, Data, Indices>> _grad_extra;
    const std::shared_ptr<const DataT> _background;
    const std::vector<size_t> _grad_param_idx;
    const size_t _n_cols;
    const size_t _n_rows;
    const size_t _size;
    const CoordinateSystem& _coordsys;

    // TODO: Make static if possible - static PyImage does not work on
    // Python 3.12; see https://rubinobs.atlassian.net/browse/DM-45445
    const Data _data_null_const{0, 0};
    const ImageArray<T, Data> _data_array_null_const{nullptr};
    const Indices _indices_null_const{0, 0};

    std::vector<size_t> _get_grad_param_idx() const {
        std::vector<size_t> grad_param_idx;
        std::set<size_t> grad_param_idx_uniq;
        for (size_t g = 0; g < _n_gaussians; ++g) {
            for (size_t p = 0; p < N_PARAMS_GAUSS2D; ++p) {
                grad_param_idx_uniq.insert(_grad_param_map->get_value(g, p));
            }
            if (_do_extra) _grad_extra->add_index_to_set(g, grad_param_idx_uniq);
        }
        std::copy(grad_param_idx_uniq.begin(), grad_param_idx_uniq.end(), std::back_inserter(grad_param_idx));
        return grad_param_idx;
    }

    // Compute Gaussian mixtures with the option to write output and/or evaluate the log likelihood
    // TODO: Reconsider whether there's a better way to do this
    // The template arguments ensure that there is no performance penalty to any of the versions of this
    // function. However, some messy validation is required as a result.
    template <OutputType output_type, bool getlikelihood, BackgroundType background_type, bool do_residual,
              GradientType gradient_type, bool do_extra>
    double gaussians_pixel_template() {
        constexpr bool writeoutput = output_type != OutputType::none;
        constexpr bool do_gradient = gradient_type != GradientType::none;
        constexpr bool is_loglike = gradient_type == GradientType::loglike;

        double background_flat = 0;
        if (background_type == BackgroundType::constant) {
            background_flat += _background->get_value_unchecked(0, 0);
        }
        const size_t n_gaussians = _gaussians.size();
        auto data_null = Data{0, 0};
        auto array_null = ImageArray<T, Data>{nullptr};

        // TODO: Return to this>IMAGE_NULL(), etc. when DM-45445 is fixed
        DataT * output_grad = is_loglike ? &(_grads->at(0)) : nullptr;
        ImageArrayT& output_jac_ref = gradient_type == GradientType::jacobian ? (*_grads) : array_null;
        size_t grad_param_idx_size = _grad_param_idx.size();

        DataT& outputref = writeoutput ? (*_output) : data_null;
        DataT& residual_ref = do_residual ? (*_residual) : data_null;
        const IndicesT& grad_param_map_ref = do_gradient ? (*_grad_param_map) : INDICES_NULL_CONST();
        const DataT& grad_param_factor_ref = do_gradient ? (*_grad_param_factor) : IMAGE_NULL_CONST();

        if (writeoutput) {
            if (_n_cols != _output->get_n_cols() || _n_rows != _output->get_n_rows()) {
                throw std::runtime_error(
                        "Data matrix dimensions [" + std::to_string(_n_cols) + ',' + std::to_string(_n_rows)
                        + "] don't match _output matrix dimensions [" + std::to_string(_output->get_n_cols())
                        + ',' + std::to_string(_output->get_n_rows()) + ']');
            }
        }

        const auto& coordsys = this->get_coordsys();
        const double x_min = coordsys.get_x_min();
        const double y_min = coordsys.get_y_min();
        const double bin_x = coordsys.get_dx1();
        const double bin_y = coordsys.get_dy2();
        const double bin_x_half = bin_x / 2.;

        detail::TermsPixelVec terms_pixel(n_gaussians);
        // These are to store pre-computed values for gradients and are unused otherwise
        const size_t ngaussgrad = n_gaussians * (do_gradient);
        detail::TermsGradientVec terms_grad(ngaussgrad);
        std::vector<detail::Weights> weights_grad(n_gaussians * (gradient_type == GradientType::loglike));

        std::vector<double> weights_conv(n_gaussians);

        auto* grad_extra_ptr = do_extra ? _grad_extra.get() : nullptr;

        for (size_t g = 0; g < n_gaussians; ++g) {
            const auto& src = _gaussians[g].get_source();
            const auto& kernel = _gaussians[g].get_kernel();
            const auto weight_src = src.get_integral_value();
            const auto weight_kernel = kernel.get_integral_value();
            // TODO: intelligently skip if weight is zero
            weights_conv[g] = weight_src * weight_kernel;
            const double cen_x = src.get_centroid_const().get_x();
            const double cen_y = src.get_centroid_const().get_y();
            const auto ell_psf = kernel.get_ellipse_const();
            const auto ell_src = src.get_ellipse_const();
            const Covariance cov_psf = Covariance(ell_psf);
            const Covariance cov_src = Covariance(ell_src);
            try {
                const auto cov_conv = cov_src.make_convolution(cov_psf);
                const auto ell_conv = Ellipse(*cov_conv);

                // Deliberately omit weights for now
                const detail::Terms terms = detail::terms_from_covar(bin_x * bin_y, ell_conv);
                auto yvals = detail::gaussian_pixel_x_xx(cen_y, y_min, bin_y, _n_rows, terms.yy, terms.xy);

                terms_pixel[g].set(terms.weight, weight_kernel, x_min - cen_x + bin_x_half, terms.xx,
                                   std::move(yvals.x_norm), std::move(yvals.xx));
                if (do_gradient) {
                    terms_grad[g].set(terms, ell_src, ell_psf, ell_conv, std::move(yvals.x));
                }
            } catch (const std::invalid_argument& err) {
                throw std::runtime_error("Failed to convolve cov_src=" + cov_src.str() + " with "
                                         + cov_psf.str() + " due to: " + err.what());
            }
        }
        const DataT& data_ref = getlikelihood ? *_data : IMAGE_NULL_CONST();
        const DataT& sigma_inv_ref = getlikelihood ? *_sigma_inv : IMAGE_NULL_CONST();
        double loglike = 0;
        double model = 0;
        double data_pix = 0;
        double sigma_inv_pix = 0;
        double chi_pix = 0;
        detail::ValuesGauss gradients = {0, 0, 0, 0, 0, 0};
        // TODO: Consider a version with cached xy, although it doubles memory usage
        for (unsigned int i = 0; i < _n_cols; i++) {
            for (size_t g = 0; g < n_gaussians; ++g) {
                terms_pixel[g].xmc_weighted = terms_pixel[g].xmc * terms_pixel[g].weight_xx;
                terms_pixel[g].xmc_sq_norm = terms_pixel[g].xmc_weighted * terms_pixel[g].xmc;
                if constexpr (gradient_type != GradientType::none) {
                    terms_grad[g].xmc_t_xy_weight = terms_pixel[g].xmc * terms_grad[g].xy_weight;
                }
            }
            for (unsigned int j = 0; j < _n_rows; j++) {
                model = background_flat;
                if constexpr (gradient_type == GradientType::jacobian) {
                    for (size_t k = 0; k < grad_param_idx_size; ++k) {
                        output_jac_ref[_grad_param_idx[k]].set_value_unchecked(j, i, 0);
                    }
                }
                if constexpr (getlikelihood || (gradient_type == GradientType::jacobian)) {
                    data_pix = data_ref.get_value_unchecked(j, i);
                    sigma_inv_pix
                            = sigma_inv_ref.get_value_unchecked(j * _is_sigma_image, i * _is_sigma_image);
                }
                for (size_t g = 0; g < n_gaussians; ++g) {
                    model += detail::gaussian_pixel_add_all<T, Data, Indices, gradient_type, do_extra>(
                            g, j, i, weights_conv[g], sigma_inv_pix, terms_pixel, output_jac_ref,
                            grad_param_map_ref, grad_param_factor_ref, weights_grad, terms_grad, gradients,
                            grad_extra_ptr);
                }
                if constexpr (output_type == OutputType::overwrite) {
                    outputref.set_value_unchecked(j, i, model);
                } else if constexpr (output_type == OutputType::add) {
                    outputref.add_value_unchecked(j, i, model);
                }
                if constexpr (getlikelihood) {
                    static_assert(getlikelihood);
                    if constexpr ((gradient_type == GradientType::none)
                                  || (gradient_type == GradientType::jacobian)) {
                        chi_pix = (data_pix - model) * sigma_inv_pix;
                        loglike -= chi_pix * chi_pix / 2.;
                    }
                    // gaussians_pixel_add_like_grad adds to the loglike to avoid redundant calculations
                    else if constexpr (gradient_type == GradientType::loglike) {
                        detail::gaussians_pixel_add_like_grad<T, Data, Indices, do_extra>(
                                output_grad, grad_param_map_ref, grad_param_factor_ref, n_gaussians,
                                weights_grad, chi_pix, loglike, model, data_pix, sigma_inv_pix, j, i,
                                terms_pixel, terms_grad, gradients, _extra_param_map, _extra_param_factor);
                    }
                }
                if constexpr (do_residual) residual_ref.set_value_unchecked(j, i, chi_pix);
            }
            for (size_t g = 0; g < n_gaussians; ++g) terms_pixel[g].xmc += bin_x;
        }
        return loglike;
    }

    // See loglike_gaussians_pixel for docs. This is just for explicit template instantiation.
    template <OutputType output_type, bool get_likelihood, BackgroundType background_type, bool do_residual,
              bool do_extra>
    double loglike_gaussians_pixel_getlike() {
        // Simplified version for testing compilation (will not generally work at runtime)
        /*
            return gaussians_pixel_template<output_type, get_likelihood, background_type, do_residual,
                        GradientType::jacobian, do_extra>();
        */
        return (_gradienttype == GradientType::loglike)
                       ? gaussians_pixel_template<output_type, get_likelihood, background_type, do_residual,
                                                  GradientType::loglike, do_extra>()
                       : (_gradienttype == GradientType::jacobian
                                  ? gaussians_pixel_template<output_type, get_likelihood, background_type,
                                                             do_residual, GradientType::jacobian, do_extra>()
                                  : gaussians_pixel_template<output_type, get_likelihood, background_type,
                                                             do_residual, GradientType::none, do_extra>());
    }

    // See loglike_gaussians_pixel for docs. This is just for explicit template instantiation.
    template <OutputType output_type, bool get_likelihood, BackgroundType background_type>
    double loglike_gaussians_pixel_extra() {
        // Simplified version for testing compilation (will not generally work at runtime)
        /*
        return loglike_gaussians_pixel_getlike<output_type, get_likelihood, background_type, true, true>();
        */
        return _do_residual ? (_do_extra ? loglike_gaussians_pixel_getlike<output_type, get_likelihood,
                                                                           background_type, true, true>()
                                         : loglike_gaussians_pixel_getlike<output_type, get_likelihood,
                                                                           background_type, true, false>())
                            : (_do_extra ? loglike_gaussians_pixel_getlike<output_type, get_likelihood,
                                                                           background_type, false, true>()
                                         : loglike_gaussians_pixel_getlike<output_type, get_likelihood,
                                                                           background_type, false, false>());
    }

    // See loglike_gaussians_pixel for docs. This is just for explicit template instantiation.
    template <OutputType output_type>
    double loglike_gaussians_pixel_output() {
        // Simplified version for testing compilation (will not generally work at runtime)
        /*
        return loglike_gaussians_pixel_extra<output_type, true, BackgroundType::constant>();
        */
        return _get_likelihood
                       ? (_has_background
                                  ? loglike_gaussians_pixel_extra<output_type, true,
                                                                  BackgroundType::constant>()
                                  : loglike_gaussians_pixel_extra<output_type, true, BackgroundType::none>())
                       : (_has_background ? loglike_gaussians_pixel_extra<output_type, false,
                                                                          BackgroundType::constant>()
                                          : loglike_gaussians_pixel_extra<output_type, false,
                                                                          BackgroundType::none>());
    }
};

/**
 * @brief Add gaussians to an image, creating one if needed
 *
 * @tparam t The data type (e.g. float, int)
 * @tparam Data The data array class
 * @tparam Indices The index array class (usually a size_t array)
 * @param gaussians The gaussians to add
 * @param output The image to add gaussians to. If null, it will be initialized
 * @param n_rows The number of rows, if creating a new image
 * @param n_cols The number of columns, if creating a new image
 * @param coordsys The image's coordinate system
 * @param to_add Whether to add to an existing image. throws if output is nullptr.
 * @return std::shared_ptr<Data> The output pointer, assigned if originally null
 */
template <typename T, class Data, class Indices>
std::shared_ptr<Data> make_gaussians_pixel(const std::shared_ptr<const ConvolvedGaussians> gaussians,
                                           std::shared_ptr<Data> output = nullptr,
                                           const unsigned int n_rows = 0, const unsigned int n_cols = 0,
                                           const std::shared_ptr<const CoordinateSystem> coordsys = nullptr,
                                           bool to_add = false) {
    if (output == nullptr) {
        if (to_add) {
            throw std::invalid_argument("Cannot set to_add if output is nullptr");
        }
        output = std::make_shared<Data>(n_rows, n_cols, nullptr, coordsys);
    }
    auto evaluator
            = std::make_unique<GaussianEvaluator<T, Data, Indices>>(gaussians, nullptr, nullptr, output);
    evaluator->loglike_pixel(to_add);
    return output;
}
}  // namespace lsst::gauss2d
#endif
