CC = gcc
CFLAGS = -g -Wall -Wextra -pedantic -std=c99 -fopenmp
CPPFLAGS =
FC = gfortran
FFLAGS = -g -Wall -Wextra -pedantic -std=f2008ts -fopenmp

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
      utils.o
STATIC_LIB = lib${LIB_NAME}.a
SHARED_LIB = lib${LIB_NAME}.so
LIB_MODS = ${LIB_NAME}_f.mod ${LIB_NAME}_fhl.mod
LIBS = ${STATIC_LIB} ${SHARED_LIB}
PREFIX = .

TESTOBJ = example.o
TESTX = example.x
TESTFOBJ = example_f.o
TESTFX = example_f.x
TESTFOBJ_HL = example_fhl.o
TESTFX_HL = example_fhl.x

all: ${LIBS}

test: ${TESTX} ${TESTFX} ${TESTFX_HL}
	@printf "\n\nRunning ${TESTX}.\n\n"
	./${TESTX}
	@printf "\n\nCompleted.\n\nRunning ${TESTFX}\n\n"
	./${TESTFX}
	@printf "\n\nCompleted.\n\nRunning ${TESTFX_HL}\n\n"
	./${TESTFX_HL}
	@printf "\n\nCompleted.\n\nAll tests complete.\n\n"

${STATIC_LIB}: ${OBJ}
	ar rcs $@ $^
	ranlib $@

${SHARED_LIB}: ${OBJ}
	${FC} -shared -fPIC -o $@ $^

%.o: src/%.c
	${CC} ${CFLAGS} ${CPPFLAGS} -fPIC -o $@ -c $<

%.o: src/%.F90
	${FC} ${FFLAGS} -fPIC -o $@ -c $<

%.o: examples/%.c
	${CC} ${CFLAGS} ${CPPFLAGS} -Isrc -o $@ -c $<

%.o: examples/%.F90
	${FC} ${FFLAGS} -o $@ -c $<

${LIB_NAME}_fhl.o: src/${LIB_NAME}_fhl.F90 ${LIB_NAME}_f.o
	${FC} ${FFLAGS} -fPIC -o $@ -c $<

${TESTOBJ}: ${LIBS}
${TESTFOBJ}: ${LIBS}
${TESTFOBJ_HL}: ${LIBS}

${TESTX}: ${TESTOBJ} ${LIBS}
	${CC} ${CFLAGS} -o $@ $< -L. -l${LIB_NAME} -Wl,-rpath='$$ORIGIN'

${TESTFX}: ${TESTFOBJ} ${LIBS}
	${FC} ${FFLAGS} -o $@ $< -L. -l${LIB_NAME} -Wl,-rpath='$$ORIGIN'

${TESTFX_HL}: ${TESTFOBJ_HL} ${LIBS}
	${FC} ${FFLAGS} -o $@ $< -L. -l${LIB_NAME} -Wl,-rpath='$$ORIGIN'

install:
	install -d ${PREFIX}/include
	install -t ${PREFIX}/include ${LIB_MODS} src/${LIB_NAME}.h
	install -d ${PREFIX}/lib
	install -t ${PREFIX}/lib ${LIBS}

clean:
	rm -f ${LIBS} ${LIB_MODS} ${OBJ} *.mod
	rm -f ${TESTX} ${TESTOBJ} ${TESTFX} ${TESTFOBJ} ${TESTFX_HL} ${TESTFOBJ_HL}
