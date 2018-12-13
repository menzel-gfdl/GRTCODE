#include <math.h>
#include <stdint.h>
#include "debug.h"
#include "floating_point_type.h"
#include "kernels.h"
#include "kernel_utils.h"
#include "line_shape.h"
#include "RFM_voigt.h"
#include "spectral_bin.h"
#include "tips2017.h"


/** @brief Calculate integrated number densities.
    @return SUCCESS or an error code.*/
int calc_number_densities(int const num_layers, /*Number of atmospheric layers.*/
                          fp_t const * const p, /*Pressure [atm] (levels).*/
                          fp_t * const n /*Integrated number densities
                                           [cm^-2] (layers).*/
                         )
{
    fp_t const c = 2.147822334314468e+25; /*[1/(cm^2*atm)]*/
    int i;
#pragma omp parallel for default(none) private(i)
    for (i=0;i<num_layers;++i)
    {
        fp_t dp = p[i] - p[i+1];
        dp = dp >= 0.f ? dp : dp*-1.f;
        n[i] = c*dp;
    }
    return SUCCESS;
}


/** @brief Calculate layer pressures and temperatures.
    @return SUCCESS or an error code.*/
int calc_pressures_and_temperatures(int const num_layers, /*Number of atmospheric
                                                            layers.*/
                                    fp_t const * const p, /*Pressure [atm] (levels).*/
                                    fp_t const * const t, /*Temperature [K] (levels).*/
                                    fp_t * const pavg, /*Pressure [atm] (layers).*/
                                    fp_t * const tavg /*Pressure [atm] (layers).*/
                                   )
{
    int i;
#pragma omp parallel for default(none) private(i)
    for (i=0;i<num_layers;++i)
    {
        pavg[i] = 0.5f*(p[i] + p[i+1]);
        tavg[i] = 0.5f*(t[i] + t[i+1]);
    }
    return SUCCESS;
}


/** @brief Calculate partial pressures and number densities.
    @return SUCCESS or an error code.*/
int calc_partial_pressures_and_number_densities(int const num_layers, /*Number of
                                                                        atmospheric layers.*/
                                                fp_t const * const p, /*Pressure [atm] (levels).*/
                                                fp_t const * const x, /*Abundance (levels).*/
                                                fp_t const * const n, /*Integrated number densities
                                                                        [cm^-2] (layers).*/
                                                fp_t * const ps, /*Partial pressure [atm]
                                                                   (layers).*/
                                                fp_t * const ns /*Integrated molecular number
                                                                  densities [cm^-2] (layers).*/
                                               )
{
    fp_t const third = 1.f/3.f;
    fp_t const sixth = 1.f/6.f;
    int i;
#pragma omp parallel for default(none) private(i)
    for (i=0;i<num_layers;++i)
    {
        ps[i] = third*(x[i]*p[i] + x[i+1]*p[i+1])
                + sixth*(x[i]*p[i+1] + x[i+1]*p[i]);
        ns[i] = n[i]*0.5f*(x[i] + x[i+1]);
    }
    return SUCCESS;
}


/** @brief Calculate pressure-shifted line center positions.
    @return SUCCESS or an error code.*/
int calc_line_centers(uint64_t const num_lines, /*Number of molecular lines.*/
                      int const num_layers, /*Number of atmospheric layers.*/
                      fp_t const * const v0, /*Unshifted line center
                                               positions [1/cm] (lines).*/
                      fp_t const * const delta, /*Air-broadened pressure
                                                  shift [1/(cm*atm)] (lines).*/
                      fp_t const * const p, /*Pressure [atm] (layers).*/
                      fp_t * const vnn /*Pressure-shifted line center
                                         positions [1/cm] (layers,lines).*/
                     )
{
    int i;
    uint64_t j;
#pragma omp parallel for collapse(2) default(none) private(i,j)
    for (i=0;i<num_layers;++i)
    {
        for (j=0;j<num_lines;++j)
        {
            vnn[i*num_lines+j] = v0[j] + delta[j]*p[i];
        }
    }
    return SUCCESS;
}


/** @brief Calculate the total partition functions for each isotopologue.*/
int calc_partition_functions(int const num_layers, /*Number of atmospheric layers.*/
                             int const mol_id, /*Molecule id.*/
                             int const num_iso, /*Number of molecular isotopologues.*/
                             fp_t const * const t, /*Temperature [K] (layers).*/
                             fp_t * const q /*Total partition function.*/
                            )
{
    int i;
    int j;
#pragma omp parallel for collapse(2) default(none) private(i,j)
    for (i=0;i<num_layers;++i)
    {
        for (j=0;j<num_iso;++j)
        {
            q[i*num_iso+j] = 1.f/Q(mol_id,t[i],j+1);
        }
    }
    return SUCCESS;
}


/** @brief Calculate temperature-corrected line intensities.
    @return SUCCESS or an error code.*/
int calc_line_strengths(uint64_t const num_lines, /*Number of molecular lines.*/
                        int const num_layers, /*Number of atmospheric layers.*/
                        int const num_iso, /*Number of molecular
                                             isotopologues.*/
                        int const * const iso, /*Isotopologue id (lines).*/
                        fp_t const * const s0, /*Uncorrected line strengths
                                                 [1/cm] (lines).*/
                        fp_t const * const vnn, /*Line center position [1/cm]
                                                  (lines).*/
                        fp_t const * const en, /*Lower state energies [1/cm]
                                                 (lines).*/
                        fp_t const * const t, /*Temperature [K] (layers).*/
                        fp_t const * const q, /*Total partition function.*/
                        fp_t * const snn /*Temperature-corrected line
                                           strengths [1/cm] (layers,lines).*/
                       )
{
    fp_t const c2 = -1.4387686f;
    int i;
    uint64_t j;
#pragma omp parallel for collapse(2) default(none) private(i,j)
    for (i=0;i<num_layers;++i)
    {
        for (j=0;j<num_lines;++j)
        {
            snn[i*num_lines+j] = s0[j]*EXP(c2*en[j]/t[i])*
                                 (1.f - EXP(c2*vnn[j]/t[i]))*
                                 q[i*num_iso+iso[j]-1];
        }
    }
    return SUCCESS;
}


/** @brief Calculate lorentz halfwidths.
    @return SUCCESS or an error code.*/
int calc_lorentz_hw(uint64_t const num_lines, /*Number of molecular lines.*/
                    int const num_layers, /*Number of atmospheric layers.*/
                    fp_t const * const n, /*Coefficient of temperature
                                            dependence of air-broadened
                                            halfwidths (lines).*/
                    fp_t const * const yair, /*Air-broadended halfwidths
                                               [1/(cm*atm)] at 296K and 1atm.*/
                    fp_t const * const yself, /*Self-broadended halfwidths
                                               [1/(cm*atm)] at 296K and 1atm.*/
                    fp_t const * const t, /*Temperature [K] (layers).*/
                    fp_t const * const p, /*Pressure [atm] (layers).*/
                    fp_t const * const ps, /*Partial pressure [atm] (layers.*/
                    fp_t * const gamma /*Temperature and pressure corrected
                                         lorentz halfwidths [1/cm]
                                         (layers,lines).*/
                   )
{
    fp_t const tref = 296.f; /*[K]*/
    int i;
    uint64_t j;
#pragma omp parallel for collapse(2) default(none) private(i,j)
    for (i=0;i<num_layers;++i)
    {
        for (j=0;j<num_lines;++j)
        {
            gamma[i*num_lines+j] = POW(tref/t[i],n[j])*(yair[j]*(p[i]-ps[i]) +
                                   yself[j]*ps[i]);
        }
    }
    return SUCCESS;
}


/** @brief Calculate doppler halfwidths.
    @return SUCCESS or an error code.*/
int calc_doppler_hw(uint64_t const num_lines, /*Number of molecular lines.*/
                    int const num_layers, /*Number of atmospheric layers.*/
                    fp_t const m, /*Molecular mass [g].*/
                    fp_t const * const vnn, /*Pressure-shifted line center
                                              position [1/cm] (layers,lines).*/
                    fp_t const * const t, /*Temperature [K] (layers).*/
                    fp_t * const alpha /*Doppler halfwidths [1/cm]
                                         (layers,lines).*/
                   )
{
    fp_t const sqrt_ln2 = 0.83255461115f;
    fp_t const kb = 1.380658E-16; /*[erg/K]*/
    fp_t const c = 2.99792458E10; /*[cm/s]*/
    int i;
    uint64_t j;
#pragma omp parallel for collapse(2) default(none) private(i,j)
    for (i=0;i<num_layers;++i)
    {
        for (j=0;j<num_lines;++j)
        {
            alpha[i*num_lines+j] = sqrt_ln2*vnn[i*num_lines+j]*
                                   SQRT((2.f*kb*t[i])/(m*c*c));
        }
    }
    return SUCCESS;
}


/** @brief Sort the line parameters in order of line center wavenumber.*/
int sort_lines(uint64_t const num_lines, /*Number of molecular lines.*/
               int const num_layers, /*Number of molecular_lines.*/
               fp_t * const vnn, /*Pressure-shifted line
                                   center positions [1/cm].
                                   (layers,lines).*/
               fp_t * const snn, /*Line strength [1/cm]
                                   (layers,lines).*/
               fp_t * const gamma, /*Lorentz halfwidth [1/cm]
                                     (layers,lines).*/
               fp_t * const alpha /*Doppler halfwidth [1/cm]
                                    (layers,lines).*/
              )
{
    int k;
#pragma omp parallel for default(none) private(k)
    for (k=0;k<num_layers;++k)
    {
        fp_t *v = &(vnn[k*num_lines]);
        fp_t *s = &(snn[k*num_lines]);
        fp_t *g = &(gamma[k*num_lines]);
        fp_t *a = &(alpha[k*num_lines]);
        uint64_t i;
        for (i=1;i<num_lines;++i)
        {
            fp_t value[4];
            value[0] = v[i];
            value[1] = s[i];
            value[2] = g[i];
            value[3] = a[i];
            uint64_t j = i;
            while (j > 0 && v[j-1] > value[0])
            {
                v[j] = v[j-1];
                s[j] = s[j-1];
                g[j] = g[j-1];
                a[j] = a[j-1];
                j--;
            }
            if (i != j)
            {
                v[j] = value[0];
                s[j] = value[1];
                g[j] = value[2];
                a[j] = value[3];
            }
        }
    }
    return SUCCESS;
}


/** @brief Calculate optical depths.
    @return SUCCESS or an error code.*/
int calc_optical_depth_bin_sweep(uint64_t const num_lines, /*Number of molecular lines.*/
                                 int const num_layers, /*Number of atmospheric layers.*/
                                 fp_t * const vnn, /*Pressure-shifted line
                                                     center positions [1/cm].
                                                     (layers,lines).*/
                                 fp_t * const snn, /*Line strength [1/cm]
                                                     (layers,lines).*/
                                 fp_t * const gamma, /*Lorentz halfwidth [1/cm]
                                                       (layers,lines).*/
                                 fp_t * const alpha, /*Doppler halfwidth [1/cm]
                                                       (layers,lines).*/
                                 fp_t const * const n, /*Integrated number density
                                                         [cm^-2] (layers).*/
                                 SpectralBins_t bins, /*Spectral bins.*/
                                 fp_t * const tau /*Optical depth (layer,wavenumber).*/
                                )
{
    int i;
    uint64_t j;
#pragma omp parallel for collapse(2) default(none) shared(bins) private(i,j)
    for (i=0;i<num_layers;++i)
    {
        for (j=0;j<bins.n;++j)
        {
            fp_t *v = &(vnn[i*num_lines]);
            fp_t *s = &(snn[i*num_lines]);
            fp_t *g = &(gamma[i*num_lines]);
            fp_t *a = &(alpha[i*num_lines]);
            uint64_t nbin_local = 1;
            uint64_t nbin_remote = 25;

            /*Find the "local" lines.*/
            uint64_t nbin = nbin_local;
            fp_t leftw = j > nbin ? bins.w[NIP*(j-nbin)] : bins.w[0];
            fp_t rightw = j >= (bins.n-1)-nbin ? bins.w[NIP*bins.n-1]
                          : bins.w[NIP*(j+nbin+1)-1];

            uint64_t left;
            uint64_t right;
            if (leftw <= v[num_lines-1] && rightw >= v[0])
            {
                uint64_t tmp;
                bracket(num_lines,
                        v,
                        leftw,
                        &left,
                        &tmp);
                bracket(num_lines - left,
                        &(v[left]),
                        rightw,
                        &tmp,
                        &right);
                right += left;

                LineShapeInputs_t in;
                in.w = bins.w0 + bins.l[j]*bins.wres;
                in.num_wpoints = bins.r[j] - bins.l[j] + 1;
                in.wres = bins.wres;
                fp_t t[in.num_wpoints];
                uint64_t k;
                for (k=left;k<=right;++k)
                {
                    in.line_center = v[k];
                    in.lorentz_hwhm = g[k];
                    in.doppler_hwhm = a[k];
                    rfm_voigt_line_shape(in,
                                         t);
                    uint64_t l;
                    for (l=bins.l[j];l<=bins.r[j];++l)
                    {
                        tau[i*bins.num_wpoints+l] += s[k]*n[i]*t[l-bins.l[j]];
                    }
                }
            }
            else if (leftw > v[num_lines-1])
            {
                left = num_lines;
            }
            else
            {
                right = (uint64_t)(-1);
            }

            /*Find the "remote" lines.*/
            nbin = nbin_remote;
            fp_t leftw_r = j > nbin ? bins.w[NIP*(j-nbin)] : bins.w[0];
            if (leftw >= v[0] && leftw_r <= v[num_lines-1])
            {
                uint64_t left_r;
                uint64_t tmp;
                bracket(left,
                        v,
                        leftw_r,
                        &left_r,
                        &tmp);
                LineShapeInputs_t in;
                in.w = bins.w[j*NIP];
                in.num_wpoints = NIP;
                in.wres = bins.w[j*NIP+1] - in.w;
                fp_t t_r[NIP];
                uint64_t k;
                for (k=left_r;k<left;++k)
                {
                    in.line_center = v[k];
                    in.lorentz_hwhm = g[k];
                    in.doppler_hwhm = a[k];
                    rfm_voigt_line_shape(in,
                                         t_r);
                    uint64_t l;
                    for (l=0;l<NIP;++l)
                    {
                        uint64_t offset = i*bins.n*NIP + j*NIP + l;
                        bins.tau[offset] += s[k]*n[i]*t_r[l];
                    }
                }
            }

            fp_t rightw_r = j >= (bins.n-1)-nbin ? bins.w[NIP*bins.n-1]
                            : bins.w[NIP*(j+nbin+1)-1];
            if (rightw <= v[num_lines-1] && rightw_r >= v[0])
            {
                uint64_t f = 0;
                if (right == (uint64_t)(-1))
                {
                    f = 1;
                }
                uint64_t right_r;
                uint64_t tmp;
                bracket(num_lines - (right + f),
                        &(v[right+f]),
                        rightw_r,
                        &tmp,
                        &right_r);
                right_r += right + f;
                LineShapeInputs_t in;
                in.w = bins.w[j*NIP];
                in.num_wpoints = NIP;
                in.wres = bins.w[j*NIP+1] - in.w;
                fp_t t_r[NIP];
                uint64_t k;
                for (k=right+1;k<=right_r;++k)
                {
                    in.line_center = v[k];
                    in.lorentz_hwhm = g[k];
                    in.doppler_hwhm = a[k];
                    rfm_voigt_line_shape(in,
                                         t_r);
                    uint64_t l;
                    for (l=0;l<NIP;++l)
                    {
                        uint64_t offset = i*bins.n*NIP + j*NIP + l;
                        bins.tau[offset] += s[k]*n[i]*t_r[l];
                    }
                }
            }
        }
    }
    return SUCCESS;
}


/** @brief Calculate optical depths.
    @return SUCCESS or an error code.*/
int calc_optical_depth_line_sweep(uint64_t const num_lines, /*Number of molecular lines.*/
                                  int const num_layers, /*Number of atmospheric layers.*/
                                  fp_t * const vnn, /*Pressure-shifted line
                                                      center positions [1/cm].
                                                      (layers,lines).*/
                                  fp_t * const snn, /*Line strength [1/cm]
                                                      (layers,lines).*/
                                  fp_t * const gamma, /*Lorentz halfwidth [1/cm]
                                                        (layers,lines).*/
                                  fp_t * const alpha, /*Doppler halfwidth [1/cm]
                                                        (layers,lines).*/
                                  fp_t const * const n, /*Integrated number density
                                                          [cm^-2] (layers).*/
                                  SpectralBins_t bins, /*Spectral bins.*/
                                  fp_t * const tau /*Optical depth (layer,wavenumber).*/
                                 )
{
    fp_t const bin_width = bins.wres*bins.ppb;
    int i;
    uint64_t j;
#pragma omp parallel for collapse(2) default(none) shared(bins) private(i,j)
    for (i=0;i<num_layers;++i)
    {
        for (j=0;j<num_lines;++j)
        {
            uint64_t o = i*num_lines + j;
            LineShapeInputs_t in;
            in.line_center = vnn[o];
            in.lorentz_hwhm = gamma[o];
            in.doppler_hwhm = alpha[o];
            in.wres = bins.wres;

            /*Local lines.*/
            fp_t wcutoff = 1.5f;
            fp_t leftw = in.line_center - wcutoff;
            if (leftw < bins.w0)
            {
                leftw = bins.w0;
            }
            uint64_t left = floor((leftw-bins.w0)/bin_width);

            fp_t rightw = in.line_center + wcutoff;
            fp_t maxw = bins.w0 + bins.num_wpoints*bins.wres;
            if (rightw > maxw)
            {
                rightw = maxw;
            }
            uint64_t right = floor((rightw-bins.w0)/bin_width);

            uint64_t k;
            for (k=left;k<=right;++k)
            {
                in.w = bins.w0 + bins.l[k]*bins.wres;
                in.num_wpoints = bins.r[k] - bins.l[k] + 1;
                fp_t t[in.num_wpoints];
                rfm_voigt_line_shape(in,
                                     t);
                uint64_t l;
                for (l=bins.l[k];l<=bins.r[k];++l)
                {
#pragma omp atomic update
                    tau[i*bins.num_wpoints+l] += snn[o]*n[i]*t[l-bins.l[k]];
                }
            }

            /*Remote lines.*/
            wcutoff = 25.f;
            leftw = in.line_center - wcutoff;
            if (leftw < bins.w0)
            {
                leftw = bins.w0;
            }
            uint64_t left_r = floor((leftw-bins.w0)/bin_width);
            in.num_wpoints = NIP;
            fp_t t_r[NIP];
            for (k=left_r;k<left;++k)
            {
                in.w = bins.w[k*NIP];
                in.wres = bins.w[k*NIP+1] - in.w;
                rfm_voigt_line_shape(in,
                                     t_r);
                uint64_t l;
                for (l=0;l<NIP;++l)
                {
                    uint64_t offset = i*bins.n*NIP + k*NIP + l;
#pragma omp atomic update
                    bins.tau[offset] += snn[o]*n[i]*t_r[l];
                }
            }

            rightw = in.line_center + wcutoff;
            if (rightw > maxw)
            {
                rightw = maxw;
            }
            uint64_t right_r = floor((rightw-bins.w0)/bin_width);
            for (k=right+1;k<=right_r;++k)
            {
                in.w = bins.w[k*NIP];
                in.wres = bins.w[k*NIP+1] - in.w;
                rfm_voigt_line_shape(in,
                                     t_r);
                uint64_t l;
                for (l=0;l<NIP;++l)
                {
                    uint64_t offset = i*bins.n*NIP + k*NIP + l;
#pragma omp atomic update
                    bins.tau[offset] += snn[o]*n[i]*t_r[l];
                }
            }
        }
    }
    return SUCCESS;
}


/** @brief Calculate optical depths by sampling all lines.
    @return SUCCESS or an error code.*/
int calc_optical_depth_line_sample(uint64_t const num_lines, /*Number of molecular lines.*/
                                   int const num_layers, /*Number of atmospheric layers.*/
                                   fp_t * const vnn, /*Pressure-shifted line
                                                       center positions [1/cm].
                                                       (layers,lines).*/
                                   fp_t * const snn, /*Line strength [1/cm]
                                                       (layers,lines).*/
                                   fp_t * const gamma, /*Lorentz halfwidth [1/cm]
                                                         (layers,lines).*/
                                   fp_t * const alpha, /*Doppler halfwidth [1/cm]
                                                         (layers,lines).*/
                                   fp_t const * const n, /*Integrated number density
                                                           [cm^-2] (layers).*/
                                   SpectralBins_t const bins, /*Spectral bins.*/
                                   fp_t * const tau /*Optical depth (layer,wavenumber).*/
                                  )
{
    uint64_t const fsteps = ceil(25.f/bins.wres);
    int i;
    uint64_t j;
#pragma omp parallel for collapse(2) default(none) private(i,j)
    for (i=0;i<num_layers;++i)
    {
        for (j=0;j<num_lines;++j)
        {
            uint64_t loffset = i*num_lines + j;
            LineShapeInputs_t in;
            in.line_center = vnn[loffset];
            in.lorentz_hwhm = gamma[loffset];
            in.doppler_hwhm = alpha[loffset];
            in.wres = bins.wres;
            uint64_t fcenterid = floor((2*((in.line_center-bins.w0)/
                                       bins.wres)+1)/2);
            if (fcenterid < bins.num_wpoints)
            {
                uint64_t s = (int64_t)(fcenterid-fsteps) < 0 ? 0 :
                             fcenterid-fsteps;
                uint64_t e = fcenterid+fsteps >= bins.num_wpoints ?
                             bins.num_wpoints-1 : fcenterid+fsteps;
                in.w = s*bins.wres + bins.w0;
                in.num_wpoints = e - s + 1;
                fp_t t[in.num_wpoints];
                rfm_voigt_line_shape(in,
                                     t);
                uint64_t f;
                for (f=s;f<=e;++f)
                {
#pragma omp atomic update
                    tau[i*bins.num_wpoints+f] += snn[loffset]*n[i]*t[f-s];
                }
            }
        }
    }
    return SUCCESS;
}


/** @brief Calculate the optical depth contribution of the water vapor
           continuum.*/
int calc_water_vapor_ctm_optical_depth(uint64_t const num_wpoints,
                                       int const num_layers,
                                       fp_t * const tau,
                                       fp_t const * const CS,
                                       fp_t const * const T,
                                       fp_t const * const Ps,
                                       fp_t const * const N,
                                       fp_t const * const T0,
                                       fp_t const * const CF,
                                       fp_t const * const P,
                                       fp_t const * const T0F
                                      )
{
    fp_t const tref = 296.f;
    int i;
    uint64_t j;
#pragma omp parallel for collapse(2) default(none) private(i,j)
    for (i=0;i<num_layers;++i)
    {
        for (j=0;j<num_wpoints;++j)
        {
            tau[i*num_wpoints+j] += N[i]*(tref/T[i])*((CS[j]*Ps[i]*
                                    EXP(T0[j]*(tref-T[i]))) +
                                    (CF[j]*(P[i]-Ps[i])*
                                    EXP(T0F[j]*(tref-T[i]))));
        }
    }
    return SUCCESS;
}


/** @brief Calculate the optical depth contribution of the ozone continuum.*/
int calc_ozone_ctm_optical_depth(uint64_t const num_wpoints,
                                 int const num_layers,
                                 fp_t const * const cross_section,
                                 fp_t const * const N,
                                 fp_t * const tau
                                )
{
    int i;
    uint64_t j;
#pragma omp parallel for collapse(2) default(none) private(i,j)
    for (i=0;i<num_layers;++i)
    {
        for (j=0;j<num_wpoints;++j)
        {
            tau[i*num_wpoints+j] += N[i]*cross_section[j];
        }
    }
    return SUCCESS;
}


/** @brief Do a quadratic interpolation of line wing values in each bin
           (except for the last one.).*/
int interpolate(SpectralBins_t const bins,
                fp_t * const tau
               )
{
    not_null(tau);
    if (bins.do_interp)
    {
        int i;
        uint64_t j;
#pragma omp parallel for collapse(2) default(none) private(i,j)
        for (i=0;i<bins.num_layers;++i)
        {
            for (j=0;j<bins.n-1;++j)
            {
                fp_t *t = &(tau[i*bins.num_wpoints]);
                fp_t const *x = &(bins.w[j*NIP]);
                fp_t const *y = &(bins.tau[i*bins.isize + j*NIP]);
                bin_quad_interp(x,
                                y,
                                bins.l[j],
                                bins.r[j],
                                bins.w0,
                                bins.wres,
                                t);
            }
        }
    }
    else
    {
        int i;
        uint64_t j;
#pragma omp parallel for collapse(2) default(none) private(i,j)
        for (i=0;i<bins.num_layers;++i)
        {
            for (j=0;j<bins.n-1;++j)
            {
                fp_t *t = &(tau[i*bins.num_wpoints]);
                fp_t const *y = &(bins.tau[i*bins.isize + j*NIP]);
                bin_no_interp(bins.l[j],
                              bins.r[j],
                              y,
                              t);
            }
        }
    }
    return SUCCESS;
}


/** @brief Do a quadratic interpolation of line wing values in each bin
           (except for the last one.).*/
int interpolate_last_bin(SpectralBins_t const bins,
                         fp_t * const tau
                        )
{
    not_null(tau);
    uint64_t const j = bins.n - 1;
    if (bins.do_last_interp)
    {
        int i;
#pragma omp parallel for default(none) private(i)
        for (i=0;i<bins.num_layers;++i)
        {
            fp_t *t = &(tau[i*bins.num_wpoints]);
            fp_t const *x = &(bins.w[j*NIP]);
            fp_t const *y = &(bins.tau[i*bins.isize + j*NIP]);
            bin_quad_interp(x,
                            y,
                            bins.l[j],
                            bins.r[j],
                            bins.w0,
                            bins.wres,
                            t);
        }
    }
    else
    {
        int i;
#pragma omp parallel for default(none) private(i)
        for (i=0;i<bins.num_layers;++i)
        {
            fp_t *t = &(tau[i*bins.num_wpoints]);
            fp_t const *y = &(bins.tau[i*bins.isize + j*NIP]);
            bin_no_interp(bins.l[j],
                          bins.r[j],
                          y,
                          t);
        }
    }
    return SUCCESS;
}
