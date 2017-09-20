#!/bin/bash

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
    declare config="tests.config.gpu.r$1"
    declare -a tests=("rfmip_gpu_10col_r$1" "rfmip_gpu_25col_r$1" \
                      "rfmip_gpu_50col_r$1" "rfmip_gpu_75col_r$1" \
                      "rfmip_gpu_100col_r$1")

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

#Run tests.
declare build="0"

printf "Running grtcode at $timestamp\n\n" &>> $logfile

run_tests "1"
run_tests "0.1"
run_tests "0.01"
#run_tests "0.001"

exit 0
