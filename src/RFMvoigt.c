#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/*pi^(-1/2).*/
#define RSQRPI 0.56418958

/*(ln(2))^(1/2).*/
#define SQRLN2 0.832554611

/*---------------------------------------------------------------------------*/
/*Calculate Voigt Line shape function. The Voigt lineshape formulation:

  g(X,Y) = S * g0 * K(X,Y)
  g0 = 1/Ad * SQRT(ln2/pi)
  X = (nu - nu0)/Ad *SQRT(ln2)
  Y = Al/Ad *SQRT(ln2)
  K(X,Y) = Y/pi * 
  INT^(+infty)_(-infty){exp(-t**2)/[Y**2 + (X-t)**2]}dt

  This routine calculates the complex probability function using a
  modified version of the Humlicek algorithm (JQSRT V27 437 1982)
  accepted for publication in JQSRT 1999.
  The calculation is performed for the array of x,y pairs for a given line
  over the fine mesh points of the current wide mesh.

  Arguments:
      N      [in]      Size of the frequency and voigt line shape arrays.
      DWNO   [in]      Array of frequencies where the voigt line shape will
                           be calculated.
      WNOADJ [in]      Frequency of the line center.
      DOPADJ [in]      Doppler (gaussian) half-width for the line.
      WIDADJ [in]      Pressure (lorentzian) half-width for the line.
      STRADJ [in]      Line strength.
      K      [in,out]  Array of voigt line shape values.
*/
#ifdef __NVCC__
__host__ __device__
#endif
void voigt_shape_function(int const N,
                          double const * const DWNO,
                          double const WNOADJ,
                          float const DOPADJ,
                          float const WIDADJ,
                          float const STRADJ,
                          float *K)
{
    /*Local variables*/
    const float Y0 = 1.5;
    const float Y0PY0 = Y0 + Y0; 
    const float Y0Q = Y0*Y0;
    int I;                       /*Loop variable.*/
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
        for (I=0;I<N;I++)
        {
            XI = ((float)(DWNO[I]-WNOADJ))*REPWID;
            K[I] = STRADJ*REPWID*Y/(M_PI*(XI*XI+YQ));
        }
        return;
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

    for (I=0;I<N;I++)
    {
        XI = ((float)(DWNO[I]-WNOADJ))*REPWID;
        ABX = fabs(XI);
        XQ = ABX*ABX;
        if (ABX >= XLIM0)
        {
            K[I] = YRRTPI/(XQ + YQ);
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
            K[I] = D*Y*(A0 + XQ);
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
            K[I] = D*Y*(E0 + XQ*(E2 + XQ*(E4 + XQ)));
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
            K[I] = D*(P0 + XQ*(P2 + XQ*(P4 + XQ*(P6 + XQ*P8))));
        }
        else
        {
            YPY0 = Y + Y0;
            YPY0Q = YPY0*YPY0;
            K[I] = 0.0;
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
                    K[I] = K[I] + C[J]*(YM[J]+YP[J]) - S[J]*(XM[J]-XP[J]);
                }
            }
            else
            {
                YF = Y + Y0PY0;
                for (J=0;J<=5;J++)
                {
                    K[I] = K[I]
                           + (C[J]*(MQ[J]*MF[J]-Y0*YM[J]) + S[J]*YF*XM[J])/
                           (MQ[J]+Y0Q)
                           + (C[J]*(PQ[J]*PF[J]-Y0*YP[J]) - S[J]*YF*XP[J])/
                           (PQ[J]+Y0Q);
                }
                K[I] = Y*K[I] + exp(-XQ);
            }
        }
        K[I] = STRADJ*RSQRPI*REPWID*K[I];
    }

    return;
}

#ifdef UNIT_TEST
/*---------------------------------------------------------------------------*/
/*Unit test.*/
int main(int argc,char **argv)
{
    /*Local variables*/
    int numFreqs;             /*Number of frequencies that the line shape
                                function will be calculated at.*/
    double *freqs = NULL;     /*Array of frequencies that the line shape
                                function will be calculated at.*/
    float *Knn = NULL;        /*Array of line shape function values.*/
    double lineCenter;        /*Freqency of the line center.*/
    int wingCutOff;           /*Number of frequencies from the line center to
                                the cutoff for the line shape.*/
    float resolution;         /*Freqency resolution.*/
    float dopplerHWHM;        /*Doppler HWHM.*/
    float lorentzHWHM;        /*Lorentz HWHM.*/
    float lineStrength;       /*Line strength.*/
    int i;                    /*Loop variable.*/
    FILE *outFile = NULL;     /*Output file pointer.*/
    char *outFileName = NULL; /*Name of the output file.*/
    int ioerr;                /*I/O error code.*/

    /*Make sure that the command line arguments are correct.*/
    if (argc != 4)
    {
        fprintf(stderr,
                "\nUsage: ./<executable> <line center frequency>"
                    " <wing cutoff> <resolution>\n\n");
        exit(EXIT_FAILURE);
    }

    /*Set the line center, wing cutoff, and resolution.*/
    lineCenter = atof(argv[1]);
    wingCutOff = (int)(ceil(atof(argv[2])));
    resolution = atof(argv[3]);
    fprintf(stdout,
            "\nComputing line shape for line with:\n"
                "Line center frequency: %e\n"
                "Frequeny range:        %e to %e\n",
            lineCenter,
            lineCenter - wingCutOff*resolution,
            lineCenter + wingCutOff*resolution);

    /*Set the half-widths and line strength.*/
    dopplerHWHM = 0.0004600089;
    lorentzHWHM = 0.0001201492;
    lineStrength = 0.0000001;

    /*Malloc and set the array of frequencies.*/
    numFreqs = 2*wingCutOff + 1;
    if (freqs == NULL)
    {
        freqs = (double *)malloc(sizeof(double)*numFreqs);
        if (freqs == NULL)
        {
            fprintf(stderr,
                    "Error: Malloc failed for freqs array.\n");
            exit(EXIT_FAILURE);
        }
    }
    else
    {
        fprintf(stderr,
                "Error: Cannot malloc a non-null pointer (%p, freqs) at"
                    " %p.\n",
                (void *)freqs,
                (void *)(&freqs));
        exit(EXIT_FAILURE);
    }
    for (i=0;i<numFreqs;i++)
    {
        freqs[i] = lineCenter - resolution*(wingCutOff - i);
    }

    /*Malloc the array of line shape function values.*/
    if (Knn == NULL)
    {
        Knn = (float *)malloc(sizeof(float)*numFreqs);
        if (Knn == NULL)
        {
            fprintf(stderr,
                    "Error: Malloc failed for Knn array.\n");
            exit(EXIT_FAILURE);
        }
    }
    else
    {
        fprintf(stderr,
                "Error: Cannot malloc a non-null pointer (%p, Knn) at"
                    " %p.\n",
                (void *)Knn,
                (void *)(&Knn));
        exit(EXIT_FAILURE);
    }

    /*Calculate the line shape function values.*/
    voigt_shape_function(numFreqs,
                         freqs,
                         lineCenter,
                         dopplerHWHM,
                         lorentzHWHM,
                         lineStrength,
                         Knn);

    /*Print the resulting line shape function to a file.*/
    outFileName = "rfmCPort.test";
    outFile = fopen(outFileName,
                    "w");
    if (outFile == NULL)
    {
        fprintf(stderr,
                "Error: Failed to open file %s.\n",
                outFileName);
        exit(EXIT_FAILURE);
    }
    fprintf(outFile,
            "freq,Knn\n");
    for (i=0;i<numFreqs;i++)
    {
        fprintf(outFile,
                "%e    %e\n",
                freqs[i],
                Knn[i]);
    }
    ioerr = fclose(outFile);
    if (ioerr != 0)
    {
        fprintf(stderr,
                "Error: Failed to close file %s.\n",
                outFileName);
    }
    outFile = NULL;
    fprintf(stdout,
            "\nResults placed in file %s.\n\n",
            outFileName);
    outFileName = NULL;

    /*Free arrays.*/
    if (freqs != NULL)
    {
        free(freqs);
        freqs = NULL;
    }
    else
    {
        fprintf(stderr,
                "Error: Cannot free a null pointer (freqs) at"
                    " %p.\n",
                (void *)(&freqs));
        exit(EXIT_FAILURE);
    }
    if (Knn != NULL)
    {
        free(Knn);
        Knn = NULL;
    }
    else
    {
        fprintf(stderr,
                "Error: Cannot free a null pointer (Knn) at"
                    " %p.\n",
                (void *)(&Knn));
        exit(EXIT_FAILURE);
    }

    return EXIT_SUCCESS;
}
#endif
