#!/bin/bash

#Set platform and architecture
declare platform="gaea.c4"
declare architecture="cpu_openmp"
export OMP_NUM_THREADS=36

#Create a log file.
declare timestamp=$(date +"%Y.%m.%d-%H.%M.%S")
declare logfile="logs/grtcode.out.$timestamp"
touch "$logfile"

#Check return code.
function check {
    if [ $1 -ne 0 ]; then
        printf "Error.  See $logfile\n"
        exit 1
    fi
}

#Run grtcode.
function run_tests {
    declare config="tests.config.${platform}.r$1"
    declare -a tests=("rfmip_${architecture}_10col_r$1" "rfmip_${architecture}_25col_r$1" \
                      "rfmip_${architecture}_50col_r$1" "rfmip_${architecture}_75col_r$1" \
                      "rfmip_${architecture}_100col_r$1")

    for i in "${tests[@]}"; do
        printf "\n./run_tests.py -c $config -t $i &>> $logfile ..."
        if [ $build -eq 0 ]; then
            ./run_tests.py -c "$config" -t "$i" &>> "$logfile"
            check "$?"
            build="1"
        else
            ./run_tests.py -c "$config" -t "$i" -s &>> "$logfile"
            check "$?"
        fi
        printf "done.\n"
    done
}

#Set environment.
source /opt/cray/pe/modules/default/init/bash
source $MODULESHOME/init/bash
module use -a /ncrc/home2/fms/local/modulefiles
module unload PrgEnv-pgi PrgEnv-intel PrgEnv-gnu PrgEnv-cray
module unload cray-netcdf cray-hdf5 fre
module load PrgEnv-intel/6.0.3
module swap intel intel/16.0.3.210
module load fre/bronx-12
module load cray-hdf5/1.8.16

#Run tests.
declare build="0"

printf "Running grtcode at $timestamp on $platform ($HOSTNAME) using $architecture\n\n" &>> $logfile
printf "Using $OMP_NUM_THREADS OpenMP threads.\n" &>> $logfile

run_tests "1"
run_tests "0.1"
run_tests "0.01"
#run_tests "0.001"

printf "\nFinished.\n"

exit 0
