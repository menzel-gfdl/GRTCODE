#ifndef CLOUD_OPTICS_H
#define CLOUD_OPTICS_H


#include "debug.h"
#include "device.h"
#include "ice_cloud_optics.h"
#include "incomplete_beta.h"
#include "liquid_cloud_optics.h"
#include "optics.h"
#include "spectral_grid.h"
#include "stochastic_clouds.h"


typedef struct CloudOptics
{
    Device_t device; /**< Device to run on.*/
    int num_layers; /**< Number of layers.*/
    int num_subcolumns; /**< Maximum number of possible subcolumns.*/
    SpectralGrid_t grid; /**<Spectral grid.*/
    IncompleteBeta_t beta; /**< Incomplete beta distribution.*/
    TotalWaterPDF_t pdf; /**< Total water probability distribution function.*/
    LiquidCloudOptics_t liquid; /**< Liquid cloud optics.*/
    IceCloudOptics_t ice; /**< Ice cloud optics.*/
    int * liquid_mapping; /**< Mapping between liquid bands and the spectral grid.*/
    int * ice_mapping; /**< Mapping between ice bands and the spectral grid.*/
    fp_t * tau; /**< Optical depth (layer, band).*/
    fp_t * omega; /**< Single-scatter albedo (layer, band).*/
    fp_t * g; /**< Asymmetry factor (layer, band).*/
    fp_t * random; /**< Random numbers (layer, band).*/
    Optics_t liquid_optics; /**< Liquid cloud optics on the spectral grid.*/
    Optics_t ice_optics; /**< Ice cloud optics on the spectral grid.*/
    fp_t * temperature; /**< Temperature [K] (layer).*/
    fp_t * cc; /**< Cloud fraction (layer).*/
    fp_t * clwc; /**< Liquid cloud condensate (layer).*/
    fp_t * ciwc; /**< Ice cloud condensate (layer).*/
    fp_t * height; /**< Altitude [km] (layer).*/
    fp_t * thickness; /**< Layer thickness (layer).*/
    fp_t * alpha; /**< Overlap parameter (layer-1).*/
    fp_t * random1_d; /**< Random numbers (layer, band).*/
    fp_t * random2_d; /**< Random nubmers (layer-1, band).*/
    fp_t * liquid_condensate; /**< Liquid water condensate (band, layer).*/
    fp_t * ice_condensate; /**< Ice condensate (band, layer).*/
} CloudOptics_t;


int create_cloud_optics(CloudOptics_t * self,
                        char const * beta_path,
                        char const * liquid_path,
                        char const * ice_path,
                        int const num_layers,
                        SpectralGrid_t const grid,
                        Device_t const device);


int destroy_cloud_optics(CloudOptics_t * self);


int calculate_cloud_optics(CloudOptics_t * self,
                           int const num_layers,
                           fp_t const * pressure,
                           fp_t const * temperature,
                           fp_t const * thickness,
                           fp_t const * cloud_fraction,
                           fp_t const * liquid_content,
                           fp_t const * ice_content);


#endif
