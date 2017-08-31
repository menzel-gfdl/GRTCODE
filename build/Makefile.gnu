CC=gcc

CFLAGS=-Wall -Wextra -Wno-missing-field-initializers -g -O3 -fPIC -std=gnu99
CLIBS= -lnetcdf -lhdf5 -lhdf5_hl -lm
CFLAGS+= $(CLIBS)
#NETCDF_INCLUDES = -I${CPATH}
#INCLUDES = $(NETCDF_INCLUDES)
#CFLAGS += $(INCLUDES)
LDFLAGS= --shared

GRTCODE_REQ_BINS = grtcode.o TIPS_2011.o parseHITRANfile.o parseNetcdfRadiation.o continuum.o outputNetcdfSpec.o eval_gamma.o eval_pShift.o eval_Snn_correction.o pre_eval_Snn.o GasProps.o GaussianFuncs.o IdaVoigtFuncs.o LineShapeUtils.o LorentzFuncs.o RfmVoigtFuncs.o flux.o

OPENMP =
ifneq ($(OPENMP),)
    CFLAGS += -fopenmp
endif

.PHONY: all_grtcode all clean

all_grtcode: grtcode.x grtcodeGauss.x grtcodeLorentz.x grtcodeIdaVoigt.x

grtcode.x: $(GRTCODE_REQ_BINS) eval_profile.o
	$(CC) $(CFLAGS)  $^ -o $@

grtcodeGauss.x: $(GRTCODE_REQ_BINS) eval_profile_gauss.o
	$(CC) $(CFLAGS)  $^ -o $@

grtcodeLorentz.x: $(GRTCODE_REQ_BINS) eval_profile_lorentz.o
	$(CC) $(CFLAGS)  $^ -o $@

grtcodeIdaVoigt.x: $(GRTCODE_REQ_BINS) eval_profile_ida_voigt.o
	$(CC) $(CFLAGS)  $^ -o $@

eval_profile_gauss.o: ../src/eval_profile.c
	$(CC) $(CFLAGS) -DSKIPMAIN -DGAUSSIAN -c $< -o $@

eval_profile_lorentz.o: ../src/eval_profile.c
	$(CC) $(CFLAGS) -DSKIPMAIN -DLORENTZ -c $< -o $@

eval_profile_ida_voigt.o: ../src/eval_profile.c
	$(CC) $(CFLAGS) -DSKIPMAIN -DIDA_VOIGT -c $< -o $@

%.o: ../src/%.c
	$(CC) $(CFLAGS) -DSKIPMAIN -c $< -o $@

clean:
	-rm -f *.o
	-rm -f *.x
