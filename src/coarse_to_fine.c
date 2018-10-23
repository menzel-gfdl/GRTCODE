#include <math.h>
#include <omp.h>
#include <stdint.h>
#include "floating_point_type.h"


#include <stdio.h>

int coarse_to_fine(int const num_layers,
                   double const w0,
                   double const wres_fine,
                   double const wres_coarse,
                   uint64_t const num_wpoints_fine,
                   uint64_t const num_wpoints_coarse,
                   fp_t * const tau_fine,
                   fp_t const * const tau_coarse)
{
    int lyr;
    uint64_t wid;

/*
#pragma omp parallel for schedule(static) \
                         collapse(2) \
                         default(none) \
                         private(lyr,wid)
*/
    for (lyr=0;lyr<num_layers;++lyr)
    {
        for (wid=0;wid<num_wpoints_coarse-2;wid=wid+2)
        {
            fp_t x[3];
            x[0] = w0 + wid*wres_coarse;
            x[1] = w0 + (wid+1)*wres_coarse;
            x[2] = w0 + (wid+2)*wres_coarse;
            fp_t y[3];
            y[0] = tau_coarse[lyr*num_wpoints_coarse+wid];
            y[1] = tau_coarse[lyr*num_wpoints_coarse+wid+1];
            y[2] = tau_coarse[lyr*num_wpoints_coarse+wid+2];
            uint64_t left = roundf((wid*wres_coarse)/wres_fine);
            if (left >= num_wpoints_fine)
            {
                left = num_wpoints_fine - 1;
            }
            uint64_t right = roundf(((wid+2)*wres_coarse)/wres_fine);
            if (right >= num_wpoints_fine)
            {
                right = num_wpoints_fine - 1;
            }
            uint64_t i;
            for (i=left;i<=right;++i)
            {
                fp_t w = w0 + i*wres_fine;
                fp_t t = (w-x[1])*(w-x[2])*y[0]/((x[0]-x[1])*(x[0]-x[2])) +
                         (w-x[0])*(w-x[2])*y[1]/((x[1]-x[0])*(x[1]-x[2])) +
                         (w-x[0])*(w-x[1])*y[2]/((x[2]-x[0])*(x[2]-x[1]));
                if (t < 0.f)
                {
                    t = 0.f;
                }
                tau_fine[lyr*num_wpoints_fine+i] += t;
            }
        }
    }
    return 0;
}









#ifdef FOO
int coarse_to_fine(int const num_layers,
                   double const w0,
                   double const wres_fine,
                   double const wres_coarse,
                   uint64_t const num_wpoints_fine,
                   uint64_t const num_wpoints_coarse,
                   fp_t * const tau_fine,
                   fp_t const * const tau_coarse)
{
    int lyr;
    uint64_t wid;

/*
#pragma omp parallel for schedule(static) \
                         collapse(2) \
                         default(none) \
                         private(lyr,wid)
*/
    for (lyr=0;lyr<num_layers;++lyr)
    {
        for (wid=0;wid<num_wpoints_coarse-2;wid=wid+2)
        {
            fp_t x[3];
            x[0] = w0 + wid*wres_coarse;
            x[1] = w0 + (wid+1)*wres_coarse;
            x[2] = w0 + (wid+2)*wres_coarse;
            fp_t y[3];
            y[0] = tau_coarse[lyr*num_wpoints_coarse+wid];
            fp_t scale = 1.f/y[0];
            y[0] = 1.f;
            y[1] = tau_coarse[lyr*num_wpoints_coarse+wid+1]*scale;
            y[2] = tau_coarse[lyr*num_wpoints_coarse+wid+2]*scale;
            fp_t A;
            fp_t B;
            if ((y[1] > y[0] && y[1] > y[2]) || (y[1] < y[0] && y[1] < y[2]))
            {
                B = (y[0]-y[2])/(x[0]-x[2] -(x[0]*x[0]-x[2]*x[2])/(2.f*x[1]));
                A = -1.f*B/(2.f*x[1]);
            }
            else if (y[0] < y[2])
            {
                B = (y[1]-y[2])/(x[1]-x[2] -(x[1]*x[1]-x[2]*x[2])/(2.f*x[0]));
                A = -1.f*B/(2.f*x[0]);
            }
            else
            {
                B = (y[1]-y[0])/(x[1]-x[0] -(x[1]*x[1]-x[0]*x[0])/(2.f*x[2]));
                A = -1.f*B/(2.f*x[2]);
            }
            fp_t C = y[0] - A*x[0]*x[0] - B*x[0];

            if (isnan((double)B))
            {
                fprintf(stderr,"B: %e %e %e %e %e %e\n",x[0],x[1],x[2],
                        y[0],y[1],y[2]);
            }
            if (isnan((double)A))
            {
                fprintf(stderr,"A: %e %e %e %e %e %e\n",x[0],x[1],x[2],
                        y[0],y[1],y[2]);
            }
            if (isnan((double)C))
            {
                fprintf(stderr,"C: %e %e %e %e %e %e %e %e %e\n",
                        C,A,B,
                        x[0],x[1],x[2],
                        y[0],y[1],y[2]);
            }


            uint64_t left = roundf((wid*wres_coarse)/wres_fine);
            if (left >= num_wpoints_fine)
            {
                left = num_wpoints_fine - 1;
            }
            uint64_t right = roundf(((wid+2)*wres_coarse)/wres_fine);
            if (right >= num_wpoints_fine)
            {
                right = num_wpoints_fine - 1;
            }
            uint64_t i;
            for (i=left;i<=right;++i)
            {
                fp_t w = w0 + i*wres_fine;
                tau_fine[lyr*num_wpoints_fine+i] += (A*w*w + B*w + C)/scale;
                if (isnan((double)(tau_fine[lyr*num_wpoints_fine+i])))
                {
                    fprintf(stderr,"tau: %e %e %e %e %e %e %e\n",
                            tau_fine[lyr*num_wpoints_fine+i],
                            x[0],x[1],x[2],
                            y[0],y[1],y[2]);
                }
            }
/*

            fp_t y0 = tau_coarse[lyr*num_wpoints_coarse+wid];
            fp_t y1 = tau_coarse[lyr*num_wpoints_coarse+wid+1];
            fp_t y2 = tau_coarse[lyr*num_wpoints_coarse+wid+2];
            fp_t a = 2.f*(y0 + y2 - 2.f*y1);
            fp_t b = y2 - y1;
            tau_fine[lyr*num_wpoints_fine+i] += y0 + b*w + a*w*(w-1.f);
*/
        }
    }
    return 0;
}
#endif
