#ifndef MODEL_FIELD_H_
#define MODEL_FIELD_H_

#define FSL 64
typedef float fp_t;

enum orderings
{
    T_LAT_LON,
    T_LAT_LON_LEV,
};

struct model_field
{
    char name[FSL];
    char units[FSL];
    int ordering;
    fp_t *data;
};
typedef struct model_field model_field_t;

struct req_model_fields
{
    model_field *P;
    model_field *T;
    model_field *TSURF;
    model_field *EMIS;
    model_field *xh2o;
    model_field *xco2;
    model_field *xo3;
    model_field *xn2o;
    model_field *xco;
    model_field *xch4;
    model_field *xo2;
};
typedef struct req_model_fields req_model_fields_t;

int init_model_field(model_field_t **field,
                     char const * const name,
                     char const * const units,
                     int const ordering);

int init_req_model_fields(req_model_fields_t *fields);

#endif
