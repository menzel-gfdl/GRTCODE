#!/bin/bash

source $MODULESHOME/init/bash
module load cuda/9.1
module load intel/17.0.5.239
module load impi/5.1.2.150

GRT="grtcode.x"
RUN="run"
MAKEFILE="Makefile.nvcc"
CC="`which mpiicc`"

cd packages/radiation_solvers
make -f $MAKEFILE clean
make -f $MAKEFILE -j4 CC=$CC CFLAGS='-O3'
if [ $? -ne 0 ]; then
    printf "Radiation_solvers library build failed.\n"
    exit 1
fi
cd -

make -f $MAKEFILE clean
make -f $MAKEFILE -j6 CC=$CC CFLAGS='-O3 -qopenmp -Duse_MPI'
if [ $? -ne 0 ]; then
    printf "GRTcode build failed.\n"
    exit 1
fi

export OMP_STACKSIZE="1M"

mv $GRT $RUN
cd $RUN
time mpirun -genv I_MPI_FABRICS=dapl -np 18 ./$GRT -aINPUT/new.multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-0-3.0_none.nc -ofoo --h2octm --o3ctm \
 HITFILES/01_hit12.par HITFILES/02_hit12.par HITFILES/03_hit12.par HITFILES/04_hit08.par HITFILES/06_hit12.par HITFILES/07_hit12.par \
 -r1 -w1 -W50000 --workers=8
