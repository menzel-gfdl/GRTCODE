#!/bin/bash

GRT="./grtcode.x"
RUN="../run"

MAKEFILE="Makefile.nvcc"

make -f $MAKEFILE clean
make -f $MAKEFILE
if [ $? -ne 0 ]; then
    printf "Make failed.\n"
    exit 1
fi

mv grtcode.x $RUN
cd $RUN
time $GRT -aINPUT/new.multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-0-3.0_none.nc -r0.1 -ofoo -C \
    HITFILES/01_hit12.par HITFILES/02_hit12.par HITFILES/03_hit12.par HITFILES/04_hit08.par HITFILES/06_hit12.par HITFILES/07_hit12.par
