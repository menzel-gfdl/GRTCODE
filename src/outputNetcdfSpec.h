#ifndef SET_OUTPUTNETCDFSPEC_H_
#define SET_OUTPUTNETCDFSPEC_H_

void closeOpticalDepthOutput(int const ncid);

void openOpticalDepthOutput(int * const ncid,
                            int * const varid,
                            char const FNAME[],
                            size_t const nlat,
                            size_t const nlon,
                            size_t const nlayers,
                            size_t const nF);

void writeDimensionData(int const ncid,
                        int const varid,
                        size_t const dim_size,
                        float const * const dim_data);

void writeOpticalDepthOutputByColumn(int const ncid,
                                     int const * const varid,
                                     int const t,
                                     int const lat,
                                     int const lon,
                                     int const nlayers,
                                     int const nF,
                                     float const * const spectra,
                                     float const * const fluxes,
                                     float const * const fluxes_accumulated);

#endif
