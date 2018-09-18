CC = gcc
CFLAGS = -g -Wall -Wextra -pedantic -std=c99 -fopenmp -O2
CPPFLAGS =
FC = gfortran
FFLAGS = -g -Wall -Wextra -pedantic -std=f2008ts -fopenmp -O2

LIB_NAME = molecular_lines
OBJ = ${LIB_NAME}_fhl.o \
      ${LIB_NAME}_f.o \
      ${LIB_NAME}.o \
      host_launch.o \
      ozone_continuum.o \
      water_vapor_continuum.o \
      eval_profile.o \
      eval_snn_correction.o \
      eval_pshift.o \
      eval_gamma.o \
      pre_eval_snn.o \
      integrate_layer.o \
      line_shape_utils.o \
      gas_properties.o \
      RFM_voigt.o \
      ida_voigt.o \
      doppler.o \
      lorentz.o \
      TIPS_2011.o \
      parse_HITRAN_file.o \
      parse_csv.o \
      molecules.o \
      utils.o \
      verbosity.o
STATIC_LIB = lib${LIB_NAME}.a
SHARED_LIB = lib${LIB_NAME}.so
LIB_HEADERS = src/${LIB_NAME}.h src/floating_point_type.h
LIB_MODS = ${LIB_NAME}_f.mod ${LIB_NAME}_fhl.mod
LIBS = ${STATIC_LIB} ${SHARED_LIB}
LIB_PC = ${LIB_NAME}.pc
PREFIX = $(shell pwd)
TESTX = example.x
TESTFX = example_f.x
TESTFX_HL = example_fhl.x

all: ${LIBS}

test: ${TESTX} ${TESTFX} ${TESTFX_HL}
	@./examples/run_tests --host

${STATIC_LIB}: ${OBJ}
	ar rcs $@ $^
	ranlib $@

${SHARED_LIB}: ${OBJ}
	${FC} -shared -fPIC -o $@ $^

%.o: src/%.c
	${CC} ${CFLAGS} ${CPPFLAGS} -fPIC -o $@ -c $<

${LIB_NAME}_f.o: src/${LIB_NAME}_f.F90
	${FC} ${FFLAGS} ${CPPFLAGS} -fPIC -o $@ -c $<

${LIB_NAME}_fhl.o: src/${LIB_NAME}_fhl.F90 ${LIB_NAME}_f.o
	${FC} ${FFLAGS} ${CPPFLAGS} -fPIC -o $@ -c $<

%.o: examples/%.F90 install
	${FC} ${FFLAGS} ${CPPFLAGS} -o $@ -c $< \
    $(shell pkg-config --cflags ${PREFIX}/pkg-config/${LIB_PC})

${TESTX}: examples/example.c install
	${CC} ${CFLAGS} ${CPPFLAGS} -o $@ $< \
    $(shell pkg-config --cflags ${PREFIX}/pkg-config/${LIB_PC}) \
    $(shell pkg-config --libs ${PREFIX}/pkg-config/${LIB_PC}) \
    -Wl,-rpath=$(shell pkg-config --variable=libdir ${PREFIX}/pkg-config/${LIB_PC})

${TESTFX}: examples/example_f.F90 install
	${FC} ${FFLAGS} ${CPPFLAGS} -o $@ $< \
    $(shell pkg-config --cflags ${PREFIX}/pkg-config/${LIB_PC}) \
    $(shell pkg-config --libs ${PREFIX}/pkg-config/${LIB_PC}) \
    -Wl,-rpath=$(shell pkg-config --variable=libdir ${PREFIX}/pkg-config/${LIB_PC})

${TESTFX_HL}: examples/example_fhl.F90 install
	${FC} ${FFLAGS} ${CPPFLAGS} -o $@ $< \
    $(shell pkg-config --cflags ${PREFIX}/pkg-config/${LIB_PC}) \
    $(shell pkg-config --libs ${PREFIX}/pkg-config/${LIB_PC}) \
    -Wl,-rpath=$(shell pkg-config --variable=libdir ${PREFIX}/pkg-config/${LIB_PC})

${LIB_PC}: ${LIBS}
	@echo "prefix=${PREFIX}" > $@
	@echo 'includedir=$${prefix}/include' >> $@
	@echo 'libdir=$${prefix}/lib' >> $@
	@echo "" >> $@
	@echo "Name: ${SHARED_LIB}" >> $@
	@echo "Description: The ${LIB_NAME} library." >> $@
	@echo "Version: 0.0" >> $@
	@echo 'Cflags: -I$${includedir}' >> $@
	@echo 'Libs: -L$${libdir} -l${LIB_NAME}' >> $@

install: ${LIB_PC}
	install -d ${PREFIX}/include
	install -t ${PREFIX}/include ${LIB_MODS} ${LIB_HEADERS}
	install -d ${PREFIX}/lib
	install -t ${PREFIX}/lib ${LIBS}
	install -d ${PREFIX}/pkg-config
	install -t ${PREFIX}/pkg-config ${LIB_PC}

clean:
	rm -f ${LIBS} ${LIB_PC} ${LIB_MODS} ${OBJ} *.mod
	rm -f ${TESTX} ${TESTFX} ${TESTFX_HL}
