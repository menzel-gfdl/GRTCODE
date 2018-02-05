
#ifdef DOUBLE_PRECISION
typedef double fp_t;
#else
typedef float fp_t;
#endif

static fp_t ONE = 1;
static fp_t HALF = 0.5;
static fp_t THIRD = 1./3.;

static fp_t favg(fp_t const a,
                 fp_t const b)
{
    return HALF*(a + b);
}

fp_t P_expectation(fp_t const p1,
                   fp_t const p2)
{
    return favg(p1,p2);
}

fp_t T_expectation_linear(fp_t const t1,
                          fp_t const t2);
{
    return favg(t1,t2);
}

fp_t Ps_expectation_uniform(fp_t const x,
                            fp_t const p1,
                            fp_t const p2)
{
    return favg(p1,p2)*x;
}

fp_t Ps_expectation_linear(fp_t const x1,
                           fp_t const x2,
                           fp_t const p1,
                           fp_t const p2)
{
    fp_t invpdiff = ONE/(p2 - p1);
    fp_t h = (x2 - x1)*invpdiff;
    fp_t j = (x1*p2 - x2*p1)*invpdiff;

    return h*(p2*p2*p2 - p1*p1*p1)*invpdiff*THIRD + j*favg(p1,p2);
}
