#include <math.h>
#include "floating_point_type.h"
#include "line_shape.h"
#include "RfmVoigtFuncs.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*pi^(-1/2).*/
#ifdef RSQRPI
#error
#else
#define RSQRPI 0.56418958
#endif

/*(ln(2))^(1/2).*/
#ifdef SQRLN2
#error
#else
#define SQRLN2 0.832554611
#endif


/*Calculate the Voigt line shape function using the Humlicek algorithm, as
  implemented in the Reference Forward Model (RFM).

  Arguments:
      vals [in]  Structure containing line shape input values.

  Returns:
      Voigt line shape value (cm).
*/
#ifdef __NVCC__
__host__ __device__
#endif
fp_t rfm_voigt_line_shape(LineShapeInputs_t const vals)
{
    fp_t const DWNO = vals.freq;
    fp_t const WNOADJ = vals.lineCenter;
    fp_t const WIDADJ = vals.lorHWHM;
    fp_t const DOPADJ = vals.gauHWHM;
    fp_t K;
    const float Y0 = 1.5;
    const float Y0PY0 = Y0 + Y0; 
    const float Y0Q = Y0*Y0;
    int J;                       /*Loop variable.*/
    int RG1;                     /*y polynomial flag.*/
    int RG2;                     /*y polynomial flag.*/
    int RG3;                     /*y polynomial flag.*/
    float ABX;                   /*|x|.*/
    float XQ;                    /*x^2.*/
    float YQ;                    /*y^2.*/
    float YRRTPI;                /*y/sqrt(pi).*/
    float XLIM0;                 /*|x| on region boundary.*/
    float XLIM1;                 /*|x| on region boundary.*/
    float XLIM2;                 /*|x| on region boundary.*/
    float XLIM3;                 /*|x| on region boundary.*/
    float XLIM4;                 /*|x| on region boundary.*/
    float A0;
    float D0;
    float D2;
    float E0;
    float E2;
    float E4;
    float H0;
    float H2;
    float H4;
    float H6;
    float C[6];
    float S[6];
    float T[6];
    float P0;
    float P2;
    float P4;
    float P6;
    float P8;
    float Z0;
    float Z2;
    float Z4;
    float Z6;
    float Z8;
    float XP[6];
    float XM[6];
    float YP[6];
    float YM[6];
    float MQ[6];
    float PQ[6];
    float MF[6];
    float PF[6];
    float D;
    float YF;
    float YPY0;
    float YPY0Q;
    float REPWID;
    float XI;
    float Y;

    C[0] = 1.0117281;
    C[1] = -0.75197147;
    C[2] = 0.012557727;
    C[3] = 0.010022008;
    C[4] = -0.00024206814;
    C[5] = 0.00000050084806;
    S[0] = 1.393237;
    S[1] = 0.23115241;
    S[2] = -0.15535147;
    S[3] = 0.0062183662;
    S[4] = 0.000091908299;
    S[5] = -0.00000062752596;
    T[0] = 0.31424038;
    T[1] = 0.94778839;
    T[2] = 1.5976826;
    T[3] = 2.2795071;
    T[4] = 3.0206370;
    T[5] = 3.8897249;

    REPWID = SQRLN2/DOPADJ;
    Y = REPWID*WIDADJ;
    YQ = Y*Y;

    if (Y >= 70.55)
    {
        XI = ((float)(DWNO-WNOADJ))*REPWID;
        K = REPWID*Y/(M_PI*(XI*XI+YQ));
        return K;
    }

    RG1 = 1;
    RG2 = 1;
    RG3 = 1;

    YRRTPI = Y*RSQRPI;
    XLIM0 = sqrt(15100.0 + Y*(40.0 - Y*3.6));
    if (Y >= 8.425)
    {
        XLIM1 = 0.0;
    }
    else
    {
        XLIM1 = sqrt(164.0 - Y*(4.3 + Y*1.8));
    }
    XLIM2 = 6.8 - Y;
    XLIM3 = 2.4*Y;
    XLIM4 = 18.1*Y + 1.65;

    if (Y <= 0.000001)
    {
        XLIM1 = XLIM0;
        XLIM2 = XLIM0;
    }

    XI = ((float)(DWNO-WNOADJ))*REPWID;
    ABX = fabs(XI);
    XQ = ABX*ABX;
    if (ABX >= XLIM0)
    {
        K = YRRTPI/(XQ + YQ);
    }
    else if (ABX >= XLIM1)
    {
        if (RG1 != 0)
        {
            RG1 = 0;
            A0 = YQ + 0.5;
            D0 = A0*A0;
            D2 = YQ + YQ - 1.0;
        }
        D = RSQRPI/(D0 + XQ*(D2 + XQ));
        K = D*Y*(A0 + XQ);
    }
    else if (ABX >= XLIM2)
    {
        if (RG2 != 0)
        {
            RG2 = 0;
            H0 = 0.5625 + YQ*(4.5 + YQ*(10.5 + YQ*(6.0 + YQ)));
            H2 = -4.5 + YQ*(9.0 + YQ*(6.0 + YQ*4.0));
            H4 = 10.5 - YQ*(6.0 - YQ*6.0);
            H6 = -6.0 + YQ* 4.0;
            E0 = 1.875 + YQ*(8.25 + YQ*(5.5 + YQ));
            E2 = 5.25 + YQ*(1.0 + YQ*3.0);
            E4 = 0.75*H6;
        }
        D = RSQRPI/(H0 + XQ*(H2 + XQ*(H4 + XQ*(H6 + XQ))));
        K = D*Y*(E0 + XQ*(E2 + XQ*(E4 + XQ)));
    }
    else if (ABX < XLIM3)
    {
        if (RG3 != 0)
        {
            RG3 = 0;
            Z0 = 272.1014 + Y*(1280.829 + Y*(2802.870 + Y*(3764.966
                 + Y*(3447.629 + Y*(2256.981 + Y*(1074.409 + Y*(369.1989
                 + Y*(88.26741 + Y*(13.39880 + Y)))))))));
            Z2 = 211.678 + Y*(902.3066 + Y*(1758.336 + Y*(2037.310
                 + Y*(1549.675 + Y*(793.4273 + Y*(266.2987
                 + Y*(53.59518 + Y*5.0)))))));
            Z4 = 78.86585 + Y*(308.1852 + Y*(497.3014 + Y*(479.2576
                 + Y*(269.2916 + Y*(80.39278 + Y*10.0)))));
            Z6 = 22.03523 + Y*(55.02933 + Y*(92.75679 + Y*(53.59518
                 + Y*10.0)));
            Z8 = 1.496460 + Y*(13.39880 + Y*5.0);
            P0 = 153.5168 + Y*(549.3954 + Y*(919.4955 + Y*(946.8970
                 + Y*(662.8097 + Y*(328.2151 + Y*(115.3772 + Y*(27.93941
                 + Y*(4.264678 + Y*0.3183291))))))));
            P2 = -34.16955 + Y*(-1.322256+ Y*(124.5975 + Y*(189.7730
                 + Y*(139.4665 + Y*(56.81652 + Y*(12.79458
                 + Y*1.2733163))))));
            P4 = 2.584042 + Y*(10.46332 + Y*(24.01655 + Y*(29.81482
                 + Y*(12.79568 + Y*1.9099744))));
            P6 = -0.07272979 + Y*(0.9377051 + Y*(4.266322 + Y*1.273316));
            P8 = 0.0005480304 + Y*0.3183291;
        }
        D = 1.7724538/(Z0 + XQ*(Z2 + XQ*(Z4 + XQ*(Z6 + XQ*(Z8+XQ)))));
        K = D*(P0 + XQ*(P2 + XQ*(P4 + XQ*(P6 + XQ*P8))));
    }
    else
    {
        YPY0 = Y + Y0;
        YPY0Q = YPY0*YPY0;
        K = 0.0;
        for (J=0;J<=5;J++)
        {
            D = XI - T[J];
            MQ[J] = D*D;
            MF[J] = 1.0/(MQ[J] + YPY0Q);
            XM[J] = MF[J]*D;
            YM[J] = MF[J]*YPY0;
            D = XI + T[J];
            PQ[J] = D*D;
            PF[J] = 1.0/(PQ[J] + YPY0Q);
            XP[J] = PF[J]*D;
            YP[J] = PF[J]*YPY0;
        }

        if (ABX <= XLIM4)
        {
            for (J=0;J<=5;J++)
            {
                K = K + C[J]*(YM[J]+YP[J]) - S[J]*(XM[J]-XP[J]);
            }
        }
        else
        {
            YF = Y + Y0PY0;
            for (J=0;J<=5;J++)
            {
                K = K
                    + (C[J]*(MQ[J]*MF[J]-Y0*YM[J]) + S[J]*YF*XM[J])/
                    (MQ[J]+Y0Q)
                    + (C[J]*(PQ[J]*PF[J]-Y0*YP[J]) - S[J]*YF*XP[J])/
                    (PQ[J]+Y0Q);
            }
            K = Y*K + exp(-XQ);
        }
    }
    K = RSQRPI*REPWID*K;

    return K;
}
