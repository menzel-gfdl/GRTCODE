/* GRTCODE is a GPU-able Radiative Transfer Code
 * Copyright (C) 2016  Garrett Wright
 * Modified in 2019 by Raymond Menzel
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; version 2.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "driver.h"
#include "era5.h"
#include "gas_optics.h"
#include "grtcode_utilities.h"
#include "netcdf.h"


#define alloc(p, s, t) {p = (t)malloc(sizeof(*p)*s);}
#define nc_catch(e) { \
    int e_ = e; \
    if (e_ != NC_NOERR) {\
        fprintf(stderr, "[%s: %d] %s\n", __FILE__, __LINE__, nc_strerror(e_)); \
        exit(EXIT_FAILURE); \
    }}
#ifdef SINGLE_PRECISION
#define get_var(a, b, c, d, e) nc_catch(nc_get_vara_float(a, b, c, d, e))
#else
#define get_var(a, b, c, d, e) nc_catch(nc_get_vara_double(a, b, c, d, e))
#endif


enum dimid
{
    LAT = 0,
    LON,
    LEVEL,
    TIME,
    NUM_DIMS
};


struct Output
{
    int *dimid;
    int ncid;
    int *varid;
};


static size_t nlon;
static size_t nlat;


/*Reorder from (t,z,y,x) to (t,y,x,z).*/
static void tzyx_to_tyxz(fp_t *dest, fp_t *src, int nx, int ny, int nz, int nt)
{
    int i;
    for (i=0; i<nt; ++i)
    {
        int j;
        for (j=0; j<nz; ++j)
        {
            int k;
            for (k=0; k<ny; ++k)
            {
                int m;
                for (m=0; m<nx; ++m)
                {
                    int offset1 = i*nz*ny*nx + j*ny*nx + k*nx + m;
                    int offset2 = i*ny*nx*nz + k*nx*nz + m*nz + j;
                    dest[offset2] = src[offset1];
                }
            }
        }
    }
    return;
}


/*Reorder from (t,y,x,z) to (t,z,y,x).*/
static void tyxz_to_tzyx(fp_t *dest, fp_t *src, int nx, int ny, int nz, int nt)
{
    int i;
    for (i=0; i<nt; ++i)
    {
        int j;
        for (j=0; j<ny; ++j)
        {
            int k;
            for (k=0; k<nx; ++k)
            {
                int m;
                for (m=0; m<nz; ++m)
                {
                    int offset1 = i*ny*nx*nz + j*nx*nz + k*nz + m;
                    int offset2 = i*nz*ny*nx + m*ny*nx + j*nx + k;
                    dest[offset2] = src[offset1];
                }
            }
        }
    }
    return;
}




/*Reserve memory and read in atmospheric data.*/
Atmosphere_t create_atmosphere(Parser_t *parser)
{
    /*Add/parse command line arguments.*/
/*
    parser->description = "Calculates the radiative fluxes for the NASA CIRC test cases.";
*/
    add_argument(parser, "level_file", NULL, "Input data file.", NULL);
    add_argument(parser, "single_file", NULL, "Input data file.", NULL);
    int one = 1;
    add_argument(parser, "-clean", NULL, "Run without aerosols.", NULL);
    add_argument(parser, "-clear", NULL, "Run without clouds.", NULL);
    add_argument(parser, "-H2O", NULL, "Include H2O.", NULL);
    add_argument(parser, "-h2o-ctm", NULL, "Directory containing H2O continuum files", &one);
    add_argument(parser, "-O3", NULL, "Include O3.", NULL);
    add_argument(parser, "-o3-ctm", NULL, "Directory containing O3 continuum files", &one);
    add_argument(parser, "-t", "--time-lower-bound", "Starting time index.", &one);
    add_argument(parser, "-T", "--Time-upper-bound", "Ending time index.", &one);
    add_argument(parser, "-x", "--lon-lower-bound", "Starting longitude index.", &one);
    add_argument(parser, "-X", "--lon-upper-bound", "Ending longitude index.", &one);
    add_argument(parser, "-y", "--lat-lower-bound", "Starting latitude index.", &one);
    add_argument(parser, "-Y", "--lat-upper-bound", "Ending latitude index.", &one);
    add_argument(parser, "-z", "--level-lower-bound", "Starting level index.", &one);
    add_argument(parser, "-Z", "--level-upper-bound", "Ending level index.", &one);
    parse_args(*parser);

    /*Open the level input file.*/
    char buffer[valuelen];
    get_argument(*parser, "level_file", buffer);
    int ncid;
    nc_catch(nc_open(buffer, NC_NOWRITE, &ncid));

    /*Determine the number of times.*/
    Atmosphere_t atm;
    int t = get_argument(*parser, "-t", buffer) ? atoi(buffer) : 0;
    int T;
    if (get_argument(*parser, "-T", buffer))
    {
        T = atoi(buffer);
    }
    else
    {
        int dimid;
        nc_catch(nc_inq_dimid(ncid, "time", &dimid));
        size_t num_times;
        nc_catch(nc_inq_dimlen(ncid, dimid, &num_times));
        T = (int)num_times - 1;
    }
    atm.num_times = T - t + 1;

    /*Determine the number of columns.*/
    int x = get_argument(*parser, "-x", buffer) ? atoi(buffer) : 0;
    int X;
    if (get_argument(*parser, "-X", buffer))
    {
        X = atoi(buffer);
    }
    else
    {
        int dimid;
        nc_catch(nc_inq_dimid(ncid, "lon", &dimid));
        size_t num_lon;
        nc_catch(nc_inq_dimlen(ncid, dimid, &num_lon));
        X = (int)num_lon - 1;
    }
    int y = get_argument(*parser, "-y", buffer) ? atoi(buffer) : 0;
    int Y;
    if (get_argument(*parser, "-Y", buffer))
    {
        Y = atoi(buffer);
    }
    else
    {
        int dimid;
        nc_catch(nc_inq_dimid(ncid, "lat", &dimid));
        size_t num_lat;
        nc_catch(nc_inq_dimlen(ncid, dimid, &num_lat));
        Y = (int)num_lat - 1;
    }
    nlon = X - x + 1;
    nlat = Y - y + 1;
    atm.num_columns = nlon*nlat;

    /*Determine the number of levels.*/
    int z = get_argument(*parser, "-z", buffer) ? atoi(buffer) : 0;
    int Z;
    if (get_argument(*parser, "-Z", buffer))
    {
        Z = atoi(buffer);
    }
    else
    {
        int dimid;
        nc_catch(nc_inq_dimid(ncid, "level", &dimid));
        size_t num_levels;
        nc_catch(nc_inq_dimlen(ncid, dimid, &num_levels));
        Z = (int)num_levels - 1;
    }
    atm.num_levels = Z - z + 1;
    atm.num_layers = atm.num_levels - 1;

    /*Pressure.*/
    alloc(atm.level_pressure, atm.num_times*atm.num_columns*atm.num_levels, fp_t *);
    fp_t *pressure;
    alloc(pressure, atm.num_times*atm.num_levels*nlat*nlon, fp_t *);
    int varid;
    nc_catch(nc_inq_varid(ncid, "p", &varid));
    size_t start[4] = {t, z, y, x};
    size_t count[4] = {atm.num_times, atm.num_levels, nlat, nlon};
    get_var(ncid, varid, start, count, pressure);
    tzyx_to_tyxz(atm.level_pressure, pressure, nlon, nlat, atm.num_levels, atm.num_times);
    free(pressure);
    alloc(atm.layer_pressure, atm.num_times*atm.num_columns*atm.num_layers, fp_t *);
    int i;
    for (i=0; i<atm.num_times*atm.num_columns; ++i)
    {
        int j;
        for (j=0; j<atm.num_layers; ++j)
        {
            int offset1 = i*atm.num_layers + j;
            int offset2 = i*atm.num_levels + j;
            atm.layer_pressure[offset1] = 0.5*(atm.level_pressure[offset2] +
                                               atm.level_pressure[offset2+1]);
        }
    }

    /*Temperature.*/
    alloc(atm.level_temperature, atm.num_times*atm.num_columns*atm.num_levels, fp_t *);
    fp_t *temperature;
    alloc(temperature, atm.num_times*atm.num_levels*nlat*nlon, fp_t *);
    nc_catch(nc_inq_varid(ncid, "t", &varid));
    start[0] = t; start[1] = z; start[2] = y; start[3] = x;
    count[0] = atm.num_times; count[1] = atm.num_levels; count[2] = nlat; count[3] = nlon;
    get_var(ncid, varid, start, count, temperature);
    double scale;
    double add;
    nc_catch(nc_get_att_double(ncid, varid, "scale_factor", &scale));
    nc_catch(nc_get_att_double(ncid, varid, "add_offset", &add));
    for (i=0; i<atm.num_times*atm.num_levels*nlat*nlon; ++i)
    {
        temperature[i] = temperature[i]*scale + add;
    }
    tzyx_to_tyxz(atm.level_temperature, temperature, nlon, nlat, atm.num_levels, atm.num_times);
    free(temperature);
    alloc(atm.layer_temperature, atm.num_times*atm.num_columns*atm.num_layers, fp_t *);
    for (i=0; i<atm.num_times*atm.num_columns; ++i)
    {
        int j;
        for (j=0; j<atm.num_layers; ++j)
        {
            int offset_lay = i*atm.num_layers + j;
            int offset_lev = i*atm.num_levels + j;
            fp_t *tlay = &(atm.layer_temperature[offset_lay]);
            fp_t const *tlev = &(atm.level_temperature[offset_lev]);
            fp_t const *play = &(atm.layer_pressure[offset_lay]);
            fp_t const *plev = &(atm.level_pressure[offset_lev]);
            tlay[j] = tlev[j] + (tlev[j+1] - tlev[j])*(play[j] - plev[j])/(plev[j+1] - plev[j]);
        }
    }

    /*Molecular abundances.*/
    fp_t const to_ppmv = 1.e6;
    struct MoleculeMeta
    {
        int id;
        char *flag;
        char *name;
    };
    int const num_molecules = 2;
    struct MoleculeMeta molecules[num_molecules] = {{H2O, "-H2O", "q"}, {O3, "-O3", "o3"}};
    alloc(atm.molecules, num_molecules, int *);
    atm.num_molecules = 0;
    alloc(atm.ppmv, num_molecules, fp_t **);
    fp_t *abundance;
    alloc(abundance, atm.num_times*atm.num_levels*nlat*nlon, fp_t *);
    for (i=0; i<num_molecules; ++i)
    {
        if (get_argument(*parser, molecules[i].flag, NULL))
        {
            atm.molecules[atm.num_molecules] = molecules[i].id;
            alloc(atm.ppmv[atm.num_molecules], atm.num_times*atm.num_columns*atm.num_levels, fp_t *);
            fp_t *ppmv = atm.ppmv[atm.num_molecules];
            nc_catch(nc_inq_varid(ncid, molecules[i].name, &varid));
            start[0] = t; start[1] = z; start[2] = y; start[3] = x;
            count[0] = atm.num_times; count[1] = atm.num_levels; count[2] = nlat; count[3] = nlon;
            get_var(ncid, varid, start, count, abundance);
            nc_catch(nc_get_att_double(ncid, varid, "scale_factor", &scale));
            nc_catch(nc_get_att_double(ncid, varid, "add_offset", &add));
            int j;
            for (j=0; j<atm.num_times*atm.num_levels*nlat*nlon; ++j)
            {
                abundance[i] = (abundance[i]*scale + add)*to_ppmv;
            }
            tzyx_to_tyxz(ppmv, abundance, nlon, nlat, atm.num_levels, atm.num_times);
            atm.num_molecules++;
        }
    }
    free(abundance);

    /*Molecular continua.*/
    if (!get_argument(*parser, "-h2o-ctm", atm.h2o_ctm))
    {
        snprintf(atm.h2o_ctm, valuelen, "%s", "none");
    }
    if (!get_argument(*parser, "-o3-ctm", atm.o3_ctm))
    {
        snprintf(atm.o3_ctm, valuelen, "%s", "none");
    }

    /*Close level file and open single file.*/
    nc_catch(nc_close(ncid));
    get_argument(*parser, "single_file", buffer);
    nc_catch(nc_open(buffer, NC_NOWRITE, &ncid));

    /*Surface temperature.*/
    alloc(atm.surface_temperature, atm.num_times*atm.num_columns, fp_t *);
    nc_catch(nc_inq_varid(ncid, "skt", &varid));
    start[0] = t; start[1] = y; start[2] = x; start[3] = 0;
    count[0] = atm.num_times; count[1] = nlat; count[2] = nlon; count[3] = 1;
    get_var(ncid, varid, start, count, atm.surface_temperature);
    nc_catch(nc_get_att_double(ncid, varid, "scale_factor", &scale));
    nc_catch(nc_get_att_double(ncid, varid, "add_offset", &add));
    for (i=0; i<atm.num_times*atm.num_columns; ++i)
    {
        atm.surface_temperature[i] = atm.surface_temperature[i]*scale + add;
    }

    /*Solar zenith angle.*/
    alloc(atm.solar_zenith_angle, atm.num_times*atm.num_columns, fp_t *);
    for (i=0; i<atm.num_times*atm.num_columns; ++i)
    {
        atm.solar_zenith_angle[i] = 0.5;
    }

    /*Solar irradiance.*/
    alloc(atm.total_solar_irradiance, atm.num_times*atm.num_columns, fp_t *);
    nc_catch(nc_inq_varid(ncid, "tisr", &varid));
    start[0] = t; start[1] = y; start[2] = x; start[3] = 0;
    count[0] = atm.num_times; count[1] = nlat; count[2] = nlon; count[3] = 1;
    get_var(ncid, varid, start, count, atm.total_solar_irradiance);
    nc_catch(nc_get_att_double(ncid, varid, "scale_factor", &scale));
    nc_catch(nc_get_att_double(ncid, varid, "add_offset", &add));
    for (i=0; i<atm.num_columns; ++i)
    {
        atm.total_solar_irradiance[i] = (atm.total_solar_irradiance[i]*scale + add)/
                                        atm.solar_zenith_angle[i];
    }

    /*Surface albedo.*/
    atm.albedo_grid_size = 2;
    alloc(atm.albedo_grid, atm.albedo_grid_size, fp_t *);
    fp_t const ir_uv_boundary = 10000.;
    fp_t const ir_uv_offset = 1.e-5;
    atm.albedo_grid[0] = ir_uv_boundary - ir_uv_offset;
    atm.albedo_grid[1] = ir_uv_boundary + ir_uv_offset;
    alloc(atm.surface_albedo, atm.num_times*atm.num_columns*atm.albedo_grid_size, fp_t *);
    fp_t *albedo;
    alloc(albedo, atm.num_times*atm.num_columns, fp_t *);
    nc_catch(nc_inq_varid(ncid, "alnip", &varid));
    start[0] = t; start[1] = y; start[2] = x; start[3] = 0;
    count[0] = atm.num_times; count[1] = nlat; count[2] = nlon; count[3] = 1;
    get_var(ncid, varid, start, count, albedo);
    nc_catch(nc_get_att_double(ncid, varid, "scale_factor", &scale));
    nc_catch(nc_get_att_double(ncid, varid, "add_offset", &add));
    for (i=0; i<atm.num_times*atm.num_columns; ++i)
    {
        atm.surface_albedo[i*atm.albedo_grid_size] = albedo[i]*scale + add;
    }
    nc_catch(nc_inq_varid(ncid, "aluvp", &varid));
    start[0] = t; start[1] = y; start[2] = x; start[3] = 0;
    count[0] = atm.num_times; count[1] = nlat; count[2] = nlon; count[3] = 1;
    get_var(ncid, varid, start, count, albedo);
    nc_catch(nc_get_att_double(ncid, varid, "scale_factor", &scale));
    nc_catch(nc_get_att_double(ncid, varid, "add_offset", &add));
    for (i=0; i<atm.num_times*atm.num_columns; ++i)
    {
        atm.surface_albedo[i*atm.albedo_grid_size+1] = albedo[i]*scale + add;
    }
    free(albedo);

    /*Close single file.*/
    nc_catch(nc_close(ncid));

    /*Surface emissivity. CIRC cases assume this is 1.*/
    atm.emissivity_grid_size = 2;
    alloc(atm.emissivity_grid, atm.emissivity_grid_size, fp_t *);
    atm.emissivity_grid[0] = -1.;
    atm.emissivity_grid[1] = 0.;
    alloc(atm.surface_emissivity, atm.emissivity_grid_size, fp_t *);
    atm.surface_emissivity[0] = 0.98;
    atm.surface_emissivity[1] = atm.surface_emissivity[0];

    /*Aerosols.*/
    atm.clean = 1;

    /*Clouds.*/
    atm.clear = 1;
    return atm;
}


/*Free memory for atmosphere.*/
void destroy_atmosphere(Atmosphere_t * const atm)
{
    free(atm->layer_pressure);
    free(atm->level_pressure);
    free(atm->layer_temperature);
    free(atm->level_temperature);
    free(atm->solar_zenith_angle);
    free(atm->surface_temperature);
    free(atm->total_solar_irradiance);
    free(atm->albedo_grid);
    free(atm->surface_albedo);
    free(atm->emissivity_grid);
    free(atm->surface_emissivity);
    if (!atm->clean)
    {
        free(atm->aerosol_grid);
        free(atm->aerosol_optical_depth);
        free(atm->aerosol_single_scatter_albedo);
        free(atm->aerosol_asymmetry_factor);
    }
    if (!atm->clear)
    {
        free(atm->liquid_water_droplet_radius);
        free(atm->liquid_water_path);
    }
    int i;
    for (i=0; i<atm->num_molecules; ++i)
    {
        free(atm->ppmv[i]);
    }
    free(atm->ppmv);
    free(atm->molecules);
    for (i=0; i<atm->num_cfcs; ++i)
    {
        free(atm->cfc_ppmv[i]);
    }
    free(atm->cfc_ppmv);
    free(atm->cfc);
    for (i=0; i<atm->num_cia_species; ++i)
    {
        free(atm->cia_ppmv[i]);
    }
    free(atm->cia_ppmv);
    free(atm->cia);
    free(atm->cia_species);
    return;
}


/*Add a variable to the output file.*/
static void add_flux_variable(Output_t * const o, /*Output object.*/
                              VarId_t const index, /*Variable index.*/
                              char const * const name, /*Variable name.*/
                              char const * const standard_name, /*Variable standard name.*/
                              fp_t const * const fill_value /*Fill value.*/
                             )
{
#ifdef SINGLE_PRESCISION
    nc_type type = NC_FLOAT;
#else
    nc_type type = NC_DOUBLE;
#endif
    int varid;
    nc_catch(nc_def_var(o->ncid, name, type, NUM_DIMS, o->dimid, &varid));
    char *unit = "W m-2";
    nc_catch(nc_put_att_text(o->ncid, varid, "units", strlen(unit), unit));
    nc_catch(nc_put_att_text(o->ncid, varid, "standard_name", strlen(standard_name),
                             standard_name));
    if (fill_value != NULL)
    {
#ifdef SINGLE_PRESCISION
        nc_catch(nc_put_att_float(o->ncid, varid, "_FillValue", type, 1, fill_value));
#else
        nc_catch(nc_put_att_double(o->ncid, varid, "_FillValue", type, 1, fill_value));
#endif
    }
    o->varid[index] = varid;
    return;
}


/*Create an output file and write metadata.*/
void create_flux_file(Output_t **output, char const * const filepath,
                      Atmosphere_t const * const atm)
{
    Output_t *file = (Output_t *)malloc(sizeof(*file));
    file->dimid = (int *)malloc(sizeof(*(file->dimid))*NUM_DIMS);
    file->varid = (int *)malloc(sizeof(*(file->varid))*NUM_VARS);
    nc_catch(nc_create(filepath, NC_NETCDF4, &(file->ncid)));
    nc_catch(nc_def_dim(file->ncid, "level", atm->num_levels, &(file->dimid[LEVEL])));
    add_flux_variable(file, RLU, "rlu", "upwelling_longwave_flux_in_air", NULL);
    add_flux_variable(file, RLD, "rld", "downwelling_longwave_flux_in_air", NULL);
    fp_t const zero = 0;
    add_flux_variable(file, RSU, "rsu", "upwelling_shortwave_flux_in_air", &zero);
    add_flux_variable(file, RSD, "rsd", "downwelling_shortwave_flux_in_air", &zero);
    *output = file;
    return;
}


/*Close output file.*/
void close_flux_file(Output_t * const o)
{
    nc_catch(nc_close(o->ncid));
    free(o->varid);
    free(o->dimid);
    return;
}


/*Write fluxes to the output file.*/
void write_output(Output_t *output, VarId_t id, fp_t *data, int time, int column)
{
    size_t num_levels;
    nc_catch(nc_inq_dimlen(output->ncid, output->dimid[LEVEL], &num_levels));
    size_t start = 0;
    size_t count = num_levels;
#ifdef SINGLE_PRECISION
    nc_catch(nc_put_vara_float(output->ncid, output->varid[id], &start, &count, data))
#else
    nc_catch(nc_put_vara_double(output->ncid, output->varid[id], &start, &count, data))
#endif
    (void)time;
    (void)column;
    return;
}
