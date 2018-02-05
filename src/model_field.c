#include <stdio.h>
#include <stdlib.h>
#include "debug.h"
#include "model_field.h"

int init_model_field(model_field_t **field,
                     char const * const name,
                     char const * const units,
                     int const ordering)
{
    not_null(field);
    is_null(*field);
    not_null(name);
    not_null(units);
    not_null(dims);

    model_field_t *f = malloc(sizeof(*f));
    not_null(f);
    snprintf(f->name,
             FSL,
             "%s",
             name);
    snprintf(f->units,
             FSL,
             "%s",
             units);
    f->ordering = ordering;
    f->data = NULL;
    *field = f;

    return SUCCESS;
}

int init_req_model_fields(req_model_fields_t *fields)
{
    not_null(fields);

    check(init_model_field(&(fields->P),
                           "level_pressure",
                           "atm",
                           T_LAT_LON_LEV));

    check(init_model_field(&(fields->T),
                           "level_temperature",
                           "K",
                           T_LAT_LON_LEV));

    check(init_model_field(&(fields->TSURF),
                           "surface_temperature",
                           "K",
                           T_LAT_LON));

    check(init_model_field(&(fields->EMIS),
                           "surface_emissivity",
                           "",
                           T_LAT_LON));

    check(init_model_field(&(fields->xh2o),
                           "water_vapor_abundance",
                           "",
                           T_LAT_LON_LEV));

    check(init_model_field(&(fields->xco2),
                           "carbon_dioxide_abundance",
                           "",
                           T_LAT_LON_LEV));

    check(init_model_field(&(fields->xo3),
                           "ozone_abundance",
                           "",
                           T_LAT_LON_LEV));

    check(init_model_field(&(fields->xn2o),
                           "nitrous_oxide_abundance",
                           "",
                           T_LAT_LON_LEV));

    check(init_model_field(&(fields->xco),
                           "carbon_monoxide_abundance",
                           "",
                           T_LAT_LON_LEV));

    check(init_model_field(&(fields->xch4),
                           "methane_abundance",
                           "",
                           T_LAT_LON_LEV));

    check(init_model_field(&(fields->xo2),
                           "oxygen_abundance",
                           "",
                           T_LAT_LON_LEV));

    return SUCCESS;
}
