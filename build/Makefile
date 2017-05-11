CC=cc
NVCC=nvcc

CFLAGS=' -Wall -Wextra -Wno-missing-field-initializers -g -O3 '
CLIBS= -lnetcdf -lhdf5 -lhdf5_hl
CFLAGS+= $(CLIBS)
#NETCDF_INCLUDES = -I${CPATH}
#INCLUDES = $(NETCDF_INCLUDES)
#CFLAGS += $(INCLUDES)
#-Wsign-conversion (some nvidia libs can make this a noisy warning, might be fixed now)
#LDFLAGS= --shared
NVCCFLAGS= --compiler-options $(CFLAGS) -g -O3 --use_fast_math 
NVCCLDFLAGS= --compiler-options '-shared -fPIC '
NVOPT= --maxrregcount=32
NVCCFLAGS += $(NVOPT)

GRTCODE_REQ_BINS = grtcode.o TIPS_2011.o parseHITRANfile.o cudaHelpers.o parseNetcdfRadiation.o continuum.o outputNetcdfSpec.o eval_gamma.o eval_pShift.o eval_Snn_correction.o pre_eval_Snn.o GasProps.o GaussianFuncs.o IdaVoigtFuncs.o LineShapeUtils.o LorentzFuncs.o RfmVoigtFuncs.o

.PHONY: all_grtcode all clean

all_grtcode: grtcode.x grtcodeGauss.x grtcodeLorentz.x grtcodeIdaVoigt.x

grtcode.x: $(GRTCODE_REQ_BINS) eval_profile.o
	$(NVCC) $(NVCCFLAGS)  $^ -o $@

grtcodeGauss.x: $(GRTCODE_REQ_BINS) eval_profile_gauss.o
	$(NVCC) $(NVCCFLAGS)  $^ -o $@

grtcodeLorentz.x: $(GRTCODE_REQ_BINS) eval_profile_lorentz.o
	$(NVCC) $(NVCCFLAGS)  $^ -o $@

grtcodeIdaVoigt.x: $(GRTCODE_REQ_BINS) eval_profile_ida_voigt.o
	$(NVCC) $(NVCCFLAGS)  $^ -o $@

#eval_profile.o: ../src/eval_profile.c
#	$(NVCC) -x cu -dc $(NVCCLDFLAGS) $(NVCCFLAGS) -DSKIPMAIN -c $< -o $@

eval_profile_gauss.o: ../src/eval_profile.c
	$(NVCC) -x cu -dc $(NVCCLDFLAGS) $(NVCCFLAGS) -DSKIPMAIN -DGAUSSIAN -c $< -o $@

eval_profile_lorentz.o: ../src/eval_profile.c
	$(NVCC) -x cu -dc $(NVCCLDFLAGS) $(NVCCFLAGS) -DSKIPMAIN -DLORENTZ -c $< -o $@

eval_profile_ida_voigt.o: ../src/eval_profile.c
	$(NVCC) -x cu -dc $(NVCCLDFLAGS) $(NVCCFLAGS) -DSKIPMAIN -DIDA_VOIGT -c $< -o $@

%.o: ../src/%.c
	$(NVCC) -x cu -dc $(NVCCLDFLAGS) $(NVCCFLAGS) -DSKIPMAIN -c $< -o $@

%.o: ../src/%.cu
	$(NVCC) -x cu -dc $(NVCCLDFLAGS) $(NVCCFLAGS) -DSKIPMAIN -c $< -o $@





all: libparseNetcdfRadiation.so libgrtcode.so grtcode.x parseNetcdfRadiation.x TIPS_2011.unittest.x

libgrtcode.so: grtcode.o TIPS_2011.o parseHITRANfile.o cudaHelpers.o parseNetcdfRadiation.o voigt.o outputNetcdfSpec.o continuum.o RFMvoigt.o
	$(NVCC) $(NVCCLDFLAGS) $(NVCCFLAGS) $^ -o $@

libparseNetcdfRadiation.so: parseNetcdfRadiation.o
	$(NVCC) $(NVCCLDFLAGS) $(NVCCFLAGS) $^ -o $@

TIPS_2011.unittest.x: 
	$(NVCC) -x cu $(NVCCFLAGS) ../src/TIPS_2011.c ../src/cudaHelpers.cu -o $@

parseNetcdfRadiation.x: ../src/parseNetcdfRadiation.c
	$(NVCC) $(NVCCFLAGS) $^ -o $@
clean:
	-rm -f *.o
	-rm -f *.oo
	-rm -f *.so
	-rm -f *.x
