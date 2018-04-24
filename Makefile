CC = gcc
CFLAGS = -g -Wall -Wextra -pedantic -std=c99 -fopenmp

OBJ = grtcode.o write_output.o launch.o sw_flux.o lw_flux.o output_fields.o \
      eval_profile.o RfmVoigtFuncs.o IdaVoigtFuncs.o GaussianFuncs.o \
      eval_Snn_correction.o eval_pShift.o eval_gamma.o LorentzFuncs.o \
      pre_eval_Snn.o line_shape_utils.o GasProps.o integrate_layer.o \
      TIPS_2011.o continuum.o o3_continuum.o solar_flux.o \
      continuum_helpers.o parseHITRANfile.o input_fields.o model_fields.o \
      arguments.o molecules.o parse_csv.o utils.o constants.o
EXECUTABLE = grtcode.x

INCLUDE = -Isrc \
          -Ipackages/radiation_solvers/src \
          -Ipackages/netcdf-4.6.1/include
LDFLAGS = -Lpackages/radiation_solvers \
          -Lpackages/netcdf-4.6.1/lib \
          -Lpackages/hdf5-1.10.1/lib
LIBS = -lradiation_solvers -lnetcdf -lhdf5 -lm

RPATH = -Wl,-rpath=`readlink -f packages/netcdf-4.6.1/lib` \
        -Wl,-rpath=`readlink -f packages/hdf5-1.10.1/lib` \
        -Wl,-rpath=`readlink -f packages/radiation_solvers`

all: ${EXECUTABLE}

${EXECUTABLE}: ${OBJ}
	${CC} ${CFLAGS} -o $@ $^ ${LDFLAGS} ${RPATH} ${LIBS}

%.o:src/%.c
	${CC} ${CFLAGS} ${INCLUDE} -o $@ -c $<

clean:
	rm -f *.o *.x
