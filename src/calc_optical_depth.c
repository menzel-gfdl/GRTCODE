#include <math.h>
#include <stdint.h>
#include "calc_optical_depth.h"
#include "debug.h"
#include "floating_point_type.h"
#include "line_shape.h"
#include "RFM_voigt.h"
#include "spectral_bin.h"
#include "tips2017.h"


static int sort(uint64_t const num_lines,
                fp_t * const vnn,
                fp_t * const snn,
                fp_t * const gamma,
                fp_t * const alpha)
{
    uint64_t i;
    for (i=0;i<num_lines;++i)
    {
        fp_t value[4];
        value[0] = vnn[i];
        value[1] = snn[i];
        value[2] = gamma[i];
        value[3] = alpha[i];
        uint64_t j = i;
        while (j > 0 && vnn[i-1] > value[0])
        {
            vnn[j] = vnn[j-1];
            snn[j] = snn[j-1];
            gamma[j] = gamma[j-1];
            alpha[j] = alpha[j-1];
            j--;
        }
        if (i != j)
        {
            vnn[j] = value[0];
            snn[j] = value[1];
            gamma[j] = value[2];
            alpha[j] = value[3];
        }
    }
    return SUCCESS;
}


static int bracket(uint64_t const array_size,
                   fp_t const * const array,
                   fp_t const val,
                   uint64_t * const left,
                   uint64_t * const right)
{
    uint64_t l = 0;
    uint64_t r = array_size - 1;
    if (val < array[l] || val > array[r])
    {
        *left = l;
        *right = r;
/*
        fatal(RANGE_ERR,
              "val %e is not in input array.",
              val);
*/
    }
    if (array[l] == val)
    {
        r = l;
    }
    else if (array[r] == val)
    {
        l = r;
    }
    else
    {
        while (r-l > 1)
        {
            uint64_t mid = l + (r-l)/2;
            if (array[mid] == val)
            {
                l = mid;
                r = mid;
                break;
            }
            else if (val > array[mid])
            {
                l = mid;
            }
            else
            {
                r = mid;
            }
            if (l > r || l == r)
            {
                fatal(VALUE_ERR,
                      "Something went wrong (l=%zu,r=%zu).",
                      l,
                      r);
            }
        }
    }
    *left = l;
    *right = r;
    return SUCCESS;
}


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


/** @brief Calculate temperature-corrected line intensities.
    @return SUCCESS or an error code.*/
int calc_line_strengths(uint64_t const num_lines, /*Number of molecular lines.*/
                        int const num_layers, /*Number of atmospheric layers.*/
                        int const mol_id, /*Molecule id.*/
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
                        fp_t * const snn /*Temperature-corrected line
                                           strengths [1/cm] (layers,lines).*/
                       )
{
    fp_t const c2 = -1.4387686f;
    fp_t q[num_layers*num_iso];
    int i;
    uint64_t j;

#pragma omp parallel for collapse(2) default(none) shared(q) private(i,j)
    for (i=0;i<num_layers;++i)
    {
        for (j=0;j<num_iso;++j)
        {
            q[i*num_iso+j] = 1.f/Q(mol_id,t[i],j+1);
        }
    }

#pragma omp parallel for collapse(2) default(none) shared(q) private(i,j)
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


/** @brief Calculate optical depths.
    @return SUCCESS or an error code.*/
int calc_optical_depth(uint64_t const num_lines, /*Number of molecular lines.*/
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
                       SpectralBins_t * const bins, /*Spectral bins.*/
                       fp_t * const tau /*Optical depth (layer,wavenumber).*/
                      )
{
    int i;
    for (i=0;i<num_layers;++i)
    {
        fp_t *v = &(vnn[i*num_lines]);
        fp_t *s = &(snn[i*num_lines]);
        fp_t *g = &(gamma[i*num_lines]);
        fp_t *a = &(alpha[i*num_lines]);
        check(sort(num_lines,
                   v,
                   s,
                   g,
                   a));

        fp_t t[bins->num_wpoints];
        uint64_t j;
#pragma omp parallel for default(none) private(j) shared(v,s,g,a,t,i)
        for (j=0;j<bins->n;++j)
        {
            uint64_t nbin_local = 1;
            uint64_t nbin_remote = 25;

            /*Find the "local" lines.*/
            uint64_t nbin = nbin_local;
            fp_t leftw = j > nbin ? bins->w[NIP*(j-nbin)] : bins->w[0];
            fp_t rightw = j >= (bins->n-1)-nbin ? bins->w[NIP*bins->n-1]
                          : bins->w[NIP*(j+nbin+1)-1];

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

                uint64_t k;
                for (k=left;k<=right;++k)
                {
                    LineShapeInputs_t in;
                    in.w = bins->w0 + bins->l[j]*bins->wres;
                    in.num_wpoints = bins->r[j] - bins->l[j] + 1;
                    in.wres = bins->wres;
                    in.line_center = v[k];
                    in.lorentz_hwhm = g[k];
                    in.doppler_hwhm = a[k];
                    rfm_voigt_line_shape(in,
                                         &t[bins->l[j]]);
                    uint64_t l;
                    for (l=bins->l[j];l<=bins->r[j];++l)
                    {
                        tau[i*bins->num_wpoints+l] += s[k]*n[i]*t[l];
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
            fp_t leftw_r = j > nbin ? bins->w[NIP*(j-nbin)] : bins->w[0];
            if (leftw >= v[0] && leftw_r <= v[num_lines-1])
            {
                uint64_t left_r;
                uint64_t tmp;
                bracket(left,
                        v,
                        leftw_r,
                        &left_r,
                        &tmp);
                uint64_t k;
                for (k=left_r;k<left;++k)
                {
                    LineShapeInputs_t in;
                    in.w = bins->w[j*NIP];
                    in.num_wpoints = NIP;
                    in.wres = bins->w[j*NIP+1] - in.w;
                    in.line_center = v[k];
                    in.lorentz_hwhm = g[k];
                    in.doppler_hwhm = a[k];
                    fp_t t_r[NIP];
                    rfm_voigt_line_shape(in,
                                         t_r);
                    uint64_t l;
                    for (l=0;l<NIP;++l)
                    {
                        uint64_t offset = i*bins->n*NIP + j*NIP + l;
                        bins->tau[offset] += s[k]*n[i]*t_r[l];
                    }
                }
            }

            fp_t rightw_r = j >= (bins->n-1)-nbin ? bins->w[NIP*bins->n-1]
                            : bins->w[NIP*(j+nbin+1)-1];
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
                uint64_t k;
                for (k=right+1;k<=right_r;++k)
                {
                    LineShapeInputs_t in;
                    in.w = bins->w[j*NIP];
                    in.num_wpoints = NIP;
                    in.wres = bins->w[j*NIP+1] - in.w;
                    in.line_center = v[k];
                    in.lorentz_hwhm = g[k];
                    in.doppler_hwhm = a[k];
                    fp_t t_r[NIP];
                    rfm_voigt_line_shape(in,
                                         t_r);
                    uint64_t l;
                    for (l=0;l<NIP;++l)
                    {
                        uint64_t offset = i*bins->n*NIP + j*NIP + l;
                        bins->tau[offset] += s[k]*n[i]*t_r[l];
                    }
                }
            }
        }
    }
    return SUCCESS;
}


/** @brief Calculate optical depths.
    @return SUCCESS or an error code.*/
int calc_optical_depth_2(uint64_t const num_lines, /*Number of molecular lines.*/
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
                         SpectralBins_t * const bins, /*Spectral bins.*/
                         fp_t * const tau /*Optical depth (layer,wavenumber).*/
                        )
{
#ifdef NOTDONE
    int i;
    uint64_t j;

#pragma omp parallel for collapse(2) default(none) private(i,j)
    for (i=0;i<num_layers;++i)
    {
        for (j=0;j<num_lines;++j)
        {
            uint64_t o = i*num_lines + j;
            LineShapeInputs_t in;
            in.line_center = vnn[o];
            in.lorentz_hwhm = gamma[o];
            in.doppler_hwhm = alpha[o];

            /*Local lines.*/
            fp_t wcutoff = 3.f;
            fp_t leftw = in.line_center - wcutoff;
            if (leftw < bins->w0)
            {
                leftw = bins->w0;
            }
            uint64_t left = floor((leftw-bins->w0)/bins->width);

            fp_t rightw = in.line_center + wcutoff;
            fp_t maxw = bins->w0 + bins->num_wpoints*bins->wres;
            if (rightw > maxw)
            {
                rightw = maxw;
            }
            uint64_t right = ceil((rightw-bins->w0)/bins->width);

            uint64_t k;
            for (k=left;k<=right;++k)
            {
                uint64_t l;
                for (l=bins->l[j];l<=bins->r[j];++l)
                {
                    in.w = bins->w0 + l*bins->wres;
#pragma omp atomic update
                    tau[i*bins->num_wpoints+l] += snn[o]*n[i]*
                                                  rfm_voigt_line_shape(in);
                }
            }

            /*Remote lines.*/
            wcutoff = 25.f;
            leftw = in.line_center - wcutoff;
            if (leftw < bins->w0)
            {
                leftw = bins->w0;
            }
            uint64_t left_r = floor((leftw-bins->w0)/bins->width);
            for (k=left_r;k<left;++k)
            {
                uint64_t l;
                for (l=0;l<NIP;++l)
                {
                    uint64_t offset = i*bins->n*NIP + j*NIP + l;
                    in.w = bins->w[j*NIP+l];
#pragma omp atomic update
                    bins->tau[offset] += snn[o]*n[i]*
                                         rfm_voigt_line_shape(in);
                }
            }

            rightw = in.line_center + wcutoff;
            if (rightw > maxw)
            {
                rightw = maxw;
            }
            uint64_t right_r = ceil((rightw-bins->w0)/bins->width);
            for (k=right+1;k<=right_r;++k)
            {
                uint64_t l;
                for (l=0;l<NIP;++l)
                {
                    uint64_t offset = i*bins->n*NIP + j*NIP + l;
                    in.w = bins->w[j*NIP+l];
#pragma omp atomic update
                    bins->tau[offset] += snn[o]*n[i]*
                                         rfm_voigt_line_shape(in);
                }
            }
        }
    }
#endif
    return SUCCESS;
}


/** @brief Calculate optical depths.
    @return SUCCESS or an error code.*/
int calc_optical_depth_old(uint64_t const num_lines, /*Number of molecular lines.*/
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
                           SpectralBins_t * const bins, /*Spectral bins.*/
                           fp_t * const tau /*Optical depth (layer,wavenumber).*/
                          )
{
#ifdef NOTDONE
    int const fsteps = ceil(25.f/bins->wres);
    int lyr;
    uint64_t ltid;

#pragma omp parallel for schedule(static) collapse(2) default(none) private(lyr,ltid)
    for (lyr=0;lyr<num_layers;++lyr)
    {
        for (ltid=0;ltid<num_lines;++ltid)
        {
            uint64_t loffset = lyr*num_lines + ltid;
            LineShapeInputs_t in;
            in.line_center = vnn[loffset];
            in.lorentz_hwhm = gamma[loffset];
            in.doppler_hwhm = alpha[loffset];
            unsigned int fcenterid = (2*((in.line_center-bins->w0)/bins->wres)+1)/2;
            if (fcenterid < bins->num_wpoints)
            {
                int ftid;
                for (ftid=((int)fcenterid)-fsteps;ftid<=((int)fcenterid);++ftid)
                {
                    if (ftid >= 0)
                    {
                        in.w = ftid*bins->wres + bins->w0;
#pragma omp atomic update
                        tau[lyr*bins->num_wpoints+ftid] += snn[loffset]*n[lyr]*
                                                           rfm_voigt_line_shape(in);
                    }
                }
                for (ftid=((int)fcenterid)+fsteps;ftid>((int)fcenterid);--ftid)
                {
                    if (ftid < ((int)bins->num_wpoints))
                    {
                        in.w = ftid*bins->wres + bins->w0;
#pragma omp atomic update
                        tau[lyr*bins->num_wpoints+ftid] += snn[loffset]*n[lyr]*
                                                           rfm_voigt_line_shape(in);
                    }
                }
            }
        }
    }
#endif
    return SUCCESS;
}
