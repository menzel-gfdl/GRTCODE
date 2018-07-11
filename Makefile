CC = gcc
CFLAGS = -g -Wall -Wextra -pedantic -std=c99 -fopenmp
CPPFLAGS =
FC = gfortran
FFLAGS = -g -Wall -Wextra -pedantic -std=f2008ts -fopenmp

OBJ = molecular_lines_hl.o \
      molecular_lines_ll.o \
      new.o \
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
LIB_NAME = molecular_lines
STATIC_LIB = lib${LIB_NAME}.a
SHARED_LIB = lib${LIB_NAME}.so
LIB_MOD = ${LIB_NAME}.mod
LIBS = ${STATIC_LIB} ${SHARED_LIB}
PREFIX = .

TESTOBJ = example.o
TESTX = example.x
TESTFOBJ = examplef.o
TESTFX = examplef.x
TESTFOBJ_HL = examplef_hl.o
TESTFX_HL = examplef_hl.x

all: ${LIBS}

test: ${TESTX} ${TESTFX} ${TESTFX_HL}
	./${TESTX}
	./${TESTFX}
	./${TESTFX_HL}

${STATIC_LIB}: ${OBJ}
	ar rcs $@ $^
	ranlib $@

${SHARED_LIB}: ${OBJ}
	${FC} -shared -fPIC -o $@ $^

%.o: src/%.c
	${CC} ${CFLAGS} ${CPPFLAGS} -fPIC -o $@ -c $<

%.o: src/%.F90
	${FC} ${FFLAGS} -fPIC -o $@ -c $<

%.o: tests/%.c
	${CC} ${CFLAGS} ${CPPFLAGS} -Isrc -o $@ -c $<

%.o: tests/%.F90
	${FC} ${FFLAGS} -o $@ -c $<

molecular_lines_hl.o: src/molecular_lines_hl.F90 molecular_lines_ll.o
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
	install -t ${PREFIX}/include ${LIB_MOD} src/${LIB_NAME}.h
	install -d ${PREFIX}/lib
	install -t ${PREFIX}/lib ${LIBS}

clean:
	rm -f ${LIBS} ${LIB_MOD} ${OBJ} *.mod
	rm -f ${TESTX} ${TESTOBJ} ${TESTFX} ${TESTFOBJ} ${TESTFX_HL} ${TESTFOBJ_HL}
