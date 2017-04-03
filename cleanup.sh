#!/bin/bash -f

#Clean up after the test suite is run.

echo "Cleaning up ..."

#Change into the build directory.
cd build
if [ $? -ne 0 ]; then
    echo "Error: build directory does not exist or path is wrong."
    exit 1
fi

#Get rid of any old GRTcode binaries and executables.
echo "Removing old GRTcode binary files and executables ..."
make clean
if [ $? -ne 0 ]; then
    echo "Error: make clean for grtcode.x failed."
    exit 1
fi

#Change to the run directory.
cd ../run
if [ $? -ne 0 ]; then
    echo "Error: run directory does not exist or path is wrong."
    exit 1
fi

#Remove the executable from the run directory.
rm -f grtcode.x
if [ $? -ne 0 ]; then
    echo "Error: attempt to remove grtcode.x failed."
    exit 1
fi

#Change to the RESULTS directory.
cd RESULTS
if [ $? -ne 0 ]; then
    echo "Error: RESULTS directory does not exist or path is wrong."
    exit 1
fi

#Remove results files from the RESULTS directory.
rm -f *.nc
if [ $? -ne 0 ]; then
    echo "Error: attempt to remove resulting files failed."
    exit 1
fi

#Change to the verification directory.
cd ../verification
if [ $? -ne 0 ]; then
    echo "Error: verification directory does not exist or path is wrong."
    exit 1
fi

#Change into the build directory.
cd build
if [ $? -ne 0 ]; then
    echo "Error: build directory does not exist or path is wrong."
    exit 1
fi

#Get rid of any old verification binaries and executables.
echo "Removing old verification binary files and executables ..."
make clean
if [ $? -ne 0 ]; then
    echo "Error: make clean failed for verification.x."
    exit 1
fi

#Change to the run directory.
cd ../run
if [ $? -ne 0 ]; then
    echo "Error: run directory does not exist or path is wrong."
    exit 1
fi

#Remove the executable from the run directory.
rm -f verification.x
if [ $? -ne 0 ]; then
    echo "Error: attempt to remove grtcode.x failed."
    exit 1
fi

#Change to the RESULTS directory.
cd RESULTS
if [ $? -ne 0 ]; then
    echo "Error: RESULTS directory does not exist or path is wrong."
    exit 1
fi

#Remove results files from the RESULTS directory.
rm -f *.nc.*
if [ $? -ne 0 ]; then
    echo "Error: attempt to remove resulting files failed."
    exit 1
fi

#Change to the plots directory.
cd ../../plots
if [ $? -ne 0 ]; then
    echo "Error: plots directory does not exist or path is wrong."
    exit 1
fi

#Remove .gnuplot files.
rm -f *.gnuplot
if [ $? -ne 0 ]; then
    echo "Error: attempt to remove .gnuplot files failed."
    exit 1
fi

#Remove *errors_above* files.
rm -f *errors_above*
if [ $? -ne 0]; then
    echo "Error: attempt to remove *errors_above* files failed."
    exit 1
fi

#Remove .png files.
rm -f *.png
if [ $? -ne 0 ]; then
    echo "Error: attempt to remove .png files failed."
    exit 1
fi

#Remove .gif files.
rm -f *.gif
if [ $? -ne 0 ]; then
    echo "Error: attempt to remove .gif files failed."
    exit 1
fi

echo "All done."
