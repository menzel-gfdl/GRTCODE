#Makefile for verification of Grtcode output.

CC = cc
NVCC = nvcc
CFLAGS = ' -Wall -Wextra -Wno-missing-field-initializers -g -O0 '
CLIBS = -lnetcdf -lhdf5 -lhdf5_hl
CFLAGS += $(CLIBS)
NVCCFLAGS = --compiler-options $(CFLAGS) -g -O0
NVCCLDFLAGS = --compiler-options '-shared -fPIC '
NVOPT = --maxrregcount=32
NVCCFLAGS += $(NVOPT)

build: verification.x

verification.x: verify_results.o parseOpticalDepthOutput.o
	$(NVCC) $(NVCCFLAGS) $^ -o $@

verify_results.o: verify_results.c
	$(NVCC) -x cu -dc $(NVCCLDFLAGS) $(NVCCFLAGS) -DSKIPMAIN -c $< -o $@

parseOpticalDepthOutput.o: parseOpticalDepthOutput.c
	$(NVCC) -x cu -dc $(NVCCLDFLAGS) $(NVCCFLAGS) -DSKIPMAIN -c $< -o $@

clean:
	rm -f verification.x
	rm -f verify_results.o parseOpticalDepthOutput.o
