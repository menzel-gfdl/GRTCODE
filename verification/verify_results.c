#include <stdlib.h>
#include <stdio.h>
#include "parseOpticalDepthOutput.h"

int main(void)
{
    /*Local variables*/
    char grt_filename[] = "co2_with_fix.nc"; /*File containing GRTcode output.*/
    OpticalDepth_t grt_optical_depth;        /*Optical depth values from GRTcode.*/
    char rfm_filename[] = "co2_rfm_ref.spc"; /*File containing RFM output.*/
    OpticalDepth_t rfm_optical_depth;        /*Optical depth values from RFM.*/
    size_t time_end;                         /*Number of times at which the data will be compared.*/
    size_t lat_end;                          /*Number of latitudes at which the data will be compared.*/
    size_t lon_end;                          /*Number of longitudes at which the data will be compared.*/
    size_t pfull_end;                        /*Number of pressure layers at which the data will be compared.*/
    size_t wavenumber_end;                   /*Number of wavenumbers at which the data will be compared.*/
    size_t grt_index = 0;                    /*GRT optical depth index variable for comparison.*/
    size_t rfm_index = 0;                    /*RFM optical depth index variable for comparison.*/
    float diff;                              /*Difference between the optical depth values.*/
    float percent_error;                     /*Precent error for the GRT optical depth compared to RFM.*/
    float max_error = -1.0;                  /*Maximum percent error.*/
    size_t time;                             /*Loop variable.*/
    size_t lat;                              /*Loop variable.*/
    size_t lon;                              /*Loop variable.*/
    size_t pfull;                            /*Loop variable.*/
    size_t wavenumber;                       /*Loop variable.*/

    /*Read in the optical depth from the grtcode output file.*/
    printf("\nReading in the GRT output file %s.\n",
           grt_filename);
    readOpticalDepthFromGrtcodeFile(grt_filename,
                                    &grt_optical_depth);
    printf("Read complete.\n");

    /*Read in the RFM reference data.*/
    printf("\nReading in the RFM reference data from file %s.\n",
           rfm_filename);
    readOpticalDepthFromRfmFile(rfm_filename,
                                &rfm_optical_depth,
                                48,
                                3001);
    printf("Read complete.\n");

    /*Get the array dimension bounds for the comparison.*/
    printf("\nGetting the array dimensions of the RFM and GRT comparison.\n");
    if (rfm_optical_depth.ntime <= grt_optical_depth.ntime)
    {
        time_end = rfm_optical_depth.ntime;
    }
    else
    {
        time_end = grt_optical_depth.ntime;
    }
    if (rfm_optical_depth.nlat <= grt_optical_depth.nlat)
    {
        lat_end = rfm_optical_depth.nlat;
    }
    else
    {
        lat_end = grt_optical_depth.nlat;
    }
    if (rfm_optical_depth.nlon <= grt_optical_depth.nlon)
    {
        lon_end = rfm_optical_depth.nlon;
    }
    else
    {
        lon_end = grt_optical_depth.nlon;
    }
    if (rfm_optical_depth.npfull <= grt_optical_depth.npfull)
    {
        pfull_end = rfm_optical_depth.npfull;
    }
    else
    {
        pfull_end = grt_optical_depth.npfull;
    }
    if (rfm_optical_depth.nwavenumber <= grt_optical_depth.nwavenumber)
    {
        wavenumber_end = rfm_optical_depth.nwavenumber;
    }
    else
    {
        wavenumber_end = grt_optical_depth.nwavenumber;
    }
    printf("The results will be compared at:\n");
    printf("Times:           %d - %lu\n",0,time_end-1);
    printf("Latitudes:       %d - %lu\n",0,lat_end-1);
    printf("Longitudes:      %d - %lu\n",0,lon_end-1);
    printf("Pressure layers: %d - %lu\n",0,pfull_end-1);
    printf("Wavenumbers:     %d - %lu\n",0,wavenumber_end-1);

    /*Loop through the optical depth arrays and calcluate the percent
      error as:
      abs((rfm_optical_depth - grt_optical_depth))/rfm_optical_depth.*/
/*
    grt_index = 0;
    rfm_index = 0;
    max_error = -1.0;
    for (time=0;time<time_end;time++)
    {
        for (lat=0;lat<lat_end;lat++)
        {
            for (lon=0;lon<lon_end;lon++)
            {
                for (pfull=0;pfull<pfull_end;pfull++)
                {
                    for (wavenumber=0;wavenumber<wavenumber_end;wavenumber++)
                    {
                        rfm_index = time*(lat_end*lon_end*pfull_end*wavenumber_end) +
                                    lat*(lon_end*pfull_end*wavenumber_end) +
                                    lon*(pfull_end*wavenumber_end) +
                                    pfull*(wavenumber_end) +
                                    wavenumber;
                        grt_index = time*(lat_end*lon_end*pfull_end*wavenumber_end) +
                                    lat*(lon_end*pfull_end*wavenumber_end) +
                                    lon*(pfull_end*wavenumber_end) +
                                    pfull*(wavenumber_end) +
                                    wavenumber;
                        diff = rfm_optical_depth.OpticalDepth[rfm_index] -
                               grt_optical_depth.OpticalDepth[grt_index];
                        if (diff < 0.0)
                        {
                            diff = diff*(-1.0);
                        }
                        if (rfm_optical_depth.OpticalDepth[rfm_index] != 0.0)
                        {
                            percent_error = diff/(rfm_optical_depth.OpticalDepth[rfm_index]);
                            if (percent_error > max_error)
                            {
                                max_error = percent_error;
                            }
                        }
                    }
                }
            }
        }
    }
    printf("Compared %lu elements.  Max error = %e.\n",
           (time_end*lat_end*lon_end*pfull_end*wavenumber_end),
           max_error);
    printf("Verification complete.\n");
*/

    /*Free all malloc'd data.*/
    printf("\nFreeing remaining memory.\n");
    freeOpticalDepth(&grt_optical_depth);
    freeOpticalDepth(&rfm_optical_depth);
    printf("Free complete.\n\n");

    return EXIT_SUCCESS;
}
