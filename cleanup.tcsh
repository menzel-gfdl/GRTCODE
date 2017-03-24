#!/bin/tcsh -f

#Clean up after the test suite is run.

#Handle command line arguments.
if ($#argv == 0) then
    echo "No arguments given, using default values."
else
    echo "Error: No command line arguments supported."
    exit
endif

#Print out starting test suite message.
echo "Cleaning up ..."

#Change into the build directory.
cd build
if ($?) then
    echo "Error: build directory does not exist or path is wrong."
    exit 1
endif

#Get rid of any old GRTcode binaries and executables.
echo "Removing old GRTcode binary files and executables ..."
make clean
if ($?) then
    echo "Error: make clean for grtcode.x failed."
    exit 1
endif

#Change to the run directory.
cd ../run
if ($?) then
    echo "Error: run directory does not exist or path is wrong."
    exit 1
endif

#Remove the executable from the run directory.
rm grtcode.x
if ($?) then
    echo "Error: attempt to remove grtcode.x failed."
    exit 1
endif

#Change to the RESULTS directory.
cd RESULTS
if ($?) then
    echo "Error: RESULTS directory does not exist or path is wrong."
    exit 1
endif

#Remove results files from the RESULTS directory.
rm *.nc
if ($?) then
    echo "Error: attempt to remove resulting files failed."
    exit 1
endif

#Change to the verification directory.
cd ../verification
if ($?) then
    echo "Error: verification directory does not exist or path is wrong."
    exit 1
endif

#Change into the build directory.
cd build
if ($?) then
    echo "Error: build directory does not exist or path is wrong."
    exit 1
endif

#Get rid of any old verification binaries and executables.
echo "Removing old verification binary files and executables ..."
make clean
if ($?) then
    echo "Error: make clean failed for verification.x."
    exit 1
endif

#Change to the run directory.
cd ../run
if ($?) then
    echo "Error: run directory does not exist or path is wrong."
    exit 1
endif

#Remove the executable from the run directory.
rm verification.x
if ($?) then
    echo "Error: attempt to remove grtcode.x failed."
    exit 1
endif

#Change to the RESULTS directory.
cd RESULTS
if ($?) then
    echo "Error: RESULTS directory does not exist or path is wrong."
    exit 1
endif

#Remove results files from the RESULTS directory.
rm *.nc.*
if ($?) then
    echo "Error: attempt to remove resulting files failed."
    exit 1
endif

#Change to the plots directory.
cd ../../plots
if ($?) then
    echo "Error: plots directory does not exist or path is wrong."
    exit 1
endif

#Remove .gnuplot files.
rm *.gnuplot
if ($?) then
    echo "Error: attempt to remove .gnuplot files failed."
    exit 1
endif

#Remove *errors_above* files.
rm *errors_above*
if ($?) then
    echo "Error: attempt to remove *errors_above* files failed."
    exit 1
endif

#Remove .png files.
rm *.png
if ($?) then
    echo "Error: attempt to remove .png files failed."
    exit 1
endif

#Remove .gif files.
rm *.gif
if ($?) then
    echo "Error: attempt to remove .gif files failed."
    exit 1
endif

#Print all done.
echo "All done."
