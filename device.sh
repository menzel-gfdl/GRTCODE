#!/bin/bash

source $MODULESHOME/init/bash
module load cuda/9.1

GRT="grtcode.x"
RUN="run"
MAKEFILE="Makefile.nvcc"

cd packages/radiation_solvers
make -f $MAKEFILE clean
make -f $MAKEFILE -j4 CFLAGS="-O3"
if [ $? -ne 0 ]; then
    printf "Radiation_solvers library build failed.\n"
    exit 1
fi
cd -

make -f $MAKEFILE clean
make -f $MAKEFILE -j6
if [ $? -ne 0 ]; then
    printf "GRTcode build failed.\n"
    exit 1
fi

#exit 0

export OMP_STACKSIZE="1M"

mv $GRT $RUN
cd $RUN
time ./$GRT -aINPUT/new.multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-0-3.0_none.nc -ofoo --h2octm --o3ctm \
 HITFILES/01_hit12.par HITFILES/02_hit12.par HITFILES/03_hit12.par HITFILES/04_hit08.par HITFILES/06_hit12.par HITFILES/07_hit12.par \
 -r0.1 -w1 -W50000 -t0 -T1 -y0 -Y0
