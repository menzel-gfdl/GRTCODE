#!/bin/tcsh -f

#Generate test data files for the Smallsubset input atmosphere.

#Set the name of the HITRAN files.
set h2o_hitran_file = "01_hit12.par"
set co2_hitran_file = "02_hit12.par"
set o3_hitran_file = "03_hit12.par"
set n2o_hitran_file = "04_hit08.par"
set co_hitran_file = "05_hit12.par"
set ch4_hitran_file = "06_hit12.par"
set o2_hitran_file = "07_hit12.par"

#Make sure that the HITRAN files exist in the correct directory.
if (! -f "HITFILES/${h2o_hitran_file}") then
    echo "Error: the file ${h2o_hitran_file} does not exist in the HITRAN directory."
    exit
endif
if (! -f "HITFILES/${co2_hitran_file}") then
    echo "Error: the file ${co2_hitran_file} does not exist in the HITRAN directory."
    exit
endif
if (! -f "HITFILES/${o3_hitran_file}") then
    echo "Error: the file ${o3_hitran_file} does not exist in the HITRAN directory."
    exit
endif
if (! -f "HITFILES/${n2o_hitran_file}") then
    echo "Error: the file ${n2o_hitran_file} does not exist in the HITRAN directory."
    exit
endif
if (! -f "HITFILES/${co_hitran_file}") then
    echo "Error: the file ${co_hitran_file} does not exist in the HITRAN directory."
    exit
endif
if (! -f "HITFILES/${ch4_hitran_file}") then
    echo "Error: the file ${ch4_hitran_file} does not exist in the HITRAN directory."
    exit
endif
if (! -f "HITFILES/${o2_hitran_file}") then
    echo "Error: the file ${o2_hitran_file} does not exist in the HITRAN directory."
    exit
endif

#Set default values
set test_type = "smallsubset"
set atmos_data_file = "smallSubset_2t.nc"
set min_wavenumber = "1"
set max_wavenumber = "3000"
set wavenumber_res = "1"
set h2o_ppmv = "from_atmos_file"
set co2_ppmv = "400"
set o3_ppmv = "from_atmos_file"
set n2o_ppmv = "0.32"
set co_ppmv = "0"
set ch4_ppmv = "1.7"
set o2_ppmv = "200000"

#Make sure that the default atmosphere input data file exists in the correct
#directory.
if (! -f "INPUT/${atmos_data_file}") then
    echo "Error: the file ${atmos_data_file} does not exist in the INPUT directory."
    exit
endif

#Handle command line arguments.
if ($#argv == 0) then
    echo "No arguments given, using default values."
else
    echo "Error: No command line arguments supported."
    exit
endif

#Set the test_type and atmosphere input file.
#if ($argv[1] == "smallsubset") then
#    set test_type = "smallsubset"
#    set atmos_data_file = "smallSubset_2t.nc"
#else
#    echo "Error: test type must be smallsubset."
#    exit
#endif

#Print out starting test suite message.
echo "Running ${test_type} test suite ..."

#Set the output file names.
set h2o_output_file = "h2o_${test_type}_test_W_${max_wavenumber}_ppmv_${h2o_ppmv}.nc"
set co2_output_file = "co2_${test_type}_test_W_${max_wavenumber}_ppmv_${co2_ppmv}.nc"
set o3_output_file = "o3_${test_type}_test_W_${max_wavenumber}_ppmv_${o3_ppmv}.nc"
#set n2o_output_file = "n2o_${test_type}_test_W_${max_wavenumber}_ppmv_${n2o_ppmv}.nc"
#set co_output_file = "co_${test_type}_test_W_${max_wavenumber}_ppmv_${co_ppmv}.nc"
#set ch4_output_file = "ch4_${test_type}_test_W_${max_wavenumber}_ppmv_${ch4_ppmv}.nc"
#set o2_output_file = "o2_${test_type}_test_W_${max_wavenumber}_ppmv_${o2_ppmv}.nc"
set gas5_output_file = "gas5_${test_type}_test_W_${max_wavenumber}.nc"

#Get rid of any old GRTcode binaries and executables.
echo "Removing old GRTcode binary files and executables ..."
make clean

#Build the GRTcode executable.
echo "Building the GRTcode executable ..."
make -j 12

#Perform the runs.

#Water
echo "Calculating the water spectra ..."
./grtcode.x -aINPUT/$atmos_data_file -o$h2o_output_file -w$min_wavenumber -W$max_wavenumber -r$wavenumber_res -1a HITFILES/$h2o_hitran_file

#Carbon dioxide
echo "Calculating the carbon dioxide spectra ..."
./grtcode.x -aINPUT/$atmos_data_file -o$co2_output_file -w$min_wavenumber -W$max_wavenumber -r$wavenumber_res -2$co2_ppmv HITFILES/$co2_hitran_file

#Ozone
echo "Calculating the ozone spectra ..."
./grtcode.x -aINPUT/$atmos_data_file -o$o3_output_file -w$min_wavenumber -W$max_wavenumber -r$wavenumber_res -3a HITFILES/$o3_hitran_file

#Nitrous oxide
#echo "Calculating the nitrous oxide spectra ..."
#./grtcode.x -aINPUT/$atmos_data_file -o$n2o_output_file -w$min_wavenumber -W$max_wavenumber -r$wavenumber_res -4$n2o_ppmv HITFILES/$n2o_hitran_file

#Carbon monoxide
#echo "Calculating the carbon monoxide spectra ..."
#./grtcode.x -aINPUT/$atmos_data_file -o$co_output_file -w$min_wavenumber -W$max_wavenumber -r$wavenumber_res -5$co_ppmv HITFILES/$co_hitran_file

#Methane
#echo "Calculating the methane spectra ..."
#./grtcode.x -aINPUT/$atmos_data_file -o$ch4_output_file -w$min_wavenumber -W$max_wavenumber -r$wavenumber_res -6$ch4_ppmv HITFILES/$ch4_hitran_file

#Oxygen
#echo "Calculating the oxygen spectra ..."
#./grtcode.x -aINPUT/$atmos_data_file -o$o2_output_file -w$min_wavenumber -W$max_wavenumber -r$wavenumber_res -7$o2_ppmv HITFILES/$o2_hitran_file

#5 Gases
echo "Calculating the 5 gas spectra ..."
./grtcode.x -aINPUT/$atmos_data_file -o$gas5_output_file -w$min_wavenumber -W$max_wavenumber -r$wavenumber_res -1a HITFILES/$h2o_hitran_file -2$co2_ppmv HITFILES/$co2_hitran_file -3a HITFILES/$o3_hitran_file -4$n2o_ppmv HITFILES/$n2o_hitran_file -6$ch4_ppmv HITFILES/$ch4_hitran_file

#Move output files to the RESULTS directory.
mv $h2o_output_file ./RESULTS/
mv $co2_output_file ./RESULTS/
mv $o3_output_file ./RESULTS/
#mv $n2o_output_file ./RESULTS/
#mv $co_output_file ./RESULTS/
#mv $ch4_output_file ./RESULTS/
#mv $o2_output_file ./RESULTS/
mv $gas5_output_file ./RESULTS/

#Write out that the runs have finished.
echo "Runs for ${test_type} test suite complete ..."

#Test the output from the runs against the reference results.

#Change to the verification directory.
cd ./verification

#Set the name of the RFM reference results.
set h2o_rfm_reference_file = h2o/0-0.spc
set co2_rfm_reference_file = co2/0-0.spc
set o3_rfm_reference_file = o3/0-0.spc
#set n2o_rfm_reference_file = 
#set co_rfm_reference_file = 
#set ch4_rfm_reference_file = 
#set o2_rfm_reference_file = 
set gas5_rfm_reference_file = 5gas/0-0.spc

#Set the output file names.
set h2o_verification_results = "${h2o_output_file}.verification_results"
set co2_verification_results = "${co2_output_file}.verification_results"
set o3_verification_results = "${o3_output_file}.verification_results"
#set n2o_verification_results = "${n2o_output_file}.verification_results"
#set co_verification_results = "${co_output_file}.verification_results"
#set ch4_verification_results = "${ch4_output_file}.verification_results"
#set o2_verification_results = "${o2_output_file}.verification_results"
set gas5_verification_results = "${gas5_output_file}.verification_results"

#Get rid of any old verification binaries and executables.
echo "Removing old verification binary files and executables ..."
make clean

#Build the verification executable.
echo "Building the verification executable ..."
make

#Perform the verification if the reference file exists.

#Water
if ( -f "RFM_SMALLSUBSET_RESULTS/${h2o_rfm_reference_file}" ) then
    echo "Verifiying h2o results against the RFM file."
    ./verification.x -rRFM_SMALLSUBSET_RESULTS/$h2o_rfm_reference_file -o$h2o_verification_results ../RESULTS/$h2o_output_file
endif

#Carbon dioxide
if ( -f "RFM_SMALLSUBSET_RESULTS/${co2_rfm_reference_file}" ) then
    echo "Verifiying co2 results against the RFM file."
    ./verification.x -rRFM_SMALLSUBSET_RESULTS/$co2_rfm_reference_file -o$co2_verification_results ../RESULTS/$co2_output_file
endif

#Ozone
if ( -f "RFM_SMALLSUBSET_RESULTS/${o3_rfm_reference_file}" ) then
    echo "Verifiying o3 results against the RFM file."
    ./verification.x -rRFM_SMALLSUBSET_RESULTS/$o3_rfm_reference_file -o$o3_verification_results ../RESULTS/$o3_output_file
endif

#Nitrous oxide.
#if ( -f "RFM_SMALLSUBSET_RESULTS/${n2o_rfm_reference_file}" ) then
#    echo "Verifiying n2o results against the RFM file."
#    ./verification.x -rRFM_SMALLSUBSET_RESULTS/$n2o_rfm_reference_file -o$n2o_verification_results ../RESULTS/$n2o_output_file
#endif

#Carbon monoxide
#if ( -f "RFM_SMALLSUBSET_RESULTS/${co_rfm_reference_file}" ) then
#    echo "Verifiying co results against the RFM file."
#    ./verification.x -rRFM_SMALLSUBSET_RESULTS/$co_rfm_reference_file -o$co_verification_results ../RESULTS/$co_output_file
#endif

#Methane
#if ( -f "RFM_SMALLSUBSET_RESULTS/${ch4_rfm_reference_file}" ) then
#    echo "Verifiying ch4 results against the RFM file."
#    ./verification.x -rRFM_SMALLSUBSET_RESULTS/$ch4_rfm_reference_file -o$ch4_verification_results ../RESULTS/$ch4_output_file
#endif

#Oxygen
#if ( -f "RFM_SMALLSUBSET_RESULTS/${o2_rfm_reference_file}" ) then
#    echo "Verifiying o2 results against the RFM file."
#    ./verification.x -rRFM_SMALLSUBSET_RESULTS/$o2_rfm_reference_file -o$o2_verification_results ../RESULTS/$o2_output_file
#endif

#5 Gases
if ( -f "RFM_SMALLSUBSET_RESULTS/${gas5_rfm_reference_file}" ) then
    echo "Verifiying 5 gas results against the RFM file."
    ./verification.x -rRFM_SMALLSUBSET_RESULTS/$gas5_rfm_reference_file -o$gas5_verification_results ../RESULTS/$gas5_output_file
endif

#Move output files to the RESULTS directory.
mv $h2o_verification_results ./RESULTS/
mv $co2_verification_results ./RESULTS/
mv $o3_verification_results ./RESULTS/
#mv $n2o_verification_results ./RESULTS/
#mv $co_verification_results ./RESULTS/
#mv $ch4_verification_results ./RESULTS/
#mv $o2_verification_results ./RESULTS/
mv $gas5_verification_results ./RESULTS/

#Move the outputted ".gnuplot" files to the plots directory.
mv "${h2o_verification_results}.gnuplot" ./plots/
mv "${co2_verification_results}.gnuplot" ./plots/
mv "${o3_verification_results}.gnuplot" ./plots/
#mv "${n2o_verification_results}.gnuplot" ./plots/
#mv "${co_verification_results}.gnuplot" ./plots/
#mv "${ch4_verification_results}.gnuplot" ./plots/
#mv "${o2_verification_results}.gnuplot" ./plots/
mv "${gas5_verification_results}.gnuplot" ./plots/

#Write out that the verifications have finished.
echo "Verifications for ${test_type} test suite complete ..."

#Create the plots and gifs.

#Change to the plots directory.
cd ./plots

#Get rid of any old verification binaries and executables.
echo "Creating plots ..."

#Run the python script to make the plots.
./create_plots.py -d -g -f

#Print all done.
echo "All done."
