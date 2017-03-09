#!/usr/bin/env python
import re
import subprocess
import optparse
import os

#------------------------------------------------------------------------------
def createHeatMap(f):
    """
    Plot the contents of the inputted file as a heat map.
    """

    mols = f.split("_")
    mol = (mols[0].strip()).upper()
    outFile = mol + "_optical_depth_differences.png"
    title = mol + " (RFM-GRT)/RMF Optical Depth Percent Differences"
    gnuplotCommandString = ('gnuplot -e "' +
                            "reset;\n" + 
                            "set terminal png size 1000,800;\n" +
                            "set output '" + outFile + "';\n" +
                            "set title '" + title + "';\n" +
                            "unset key;\n" +
                            "set tic scale 0;\n" +
                            "set palette defined (0 'white'," +
                                ".00000000001 'blue', 1 'red');\n" +
                            "set xlabel 'Wavenumber (1/cm)';\n" +
                            "set xrange [1:3000];\n" +
                            "set ylabel 'Pressure level';\n" +
                            "set yrange [0:47];\n" +
                            "plot '" + f + "' using 5:4:9 with" +
                                " image;\n" +
                            '"')
    tmp = subprocess.Popen(gnuplotCommandString,
                           shell=True,
                           stdout=subprocess.PIPE)
    tmp.wait()
    if tmp.returncode != 0:
        print "Error: gnuplot failed."
        exit(1)

    return

#------------------------------------------------------------------------------
def findHighErrors(f,diffThreshold):
    """
    Find all the percent differences between the GRT and RFM optical depths
    that are greater than the inputted diffThreshold.
    """

    heightList = []
    currentHeight = -1
    fp = open(f,"r")
    for line in fp:
        cols = line.split()
        height = int(cols[3].strip())
        wavenumber = float(cols[4].strip())
        tauGRT = float(cols[5].strip())
        tauRFM = float(cols[6].strip())
        absoluteDiff = float(cols[7].strip())
        percentDiff = float(cols[8].strip())
        if height != currentHeight:
            if currentHeight >= 0:
                heightList.append(waveList)
            waveList = []
            currentHeight = height
        if (percentDiff >= float(diffThreshold)):
            data = (height,wavenumber,tauGRT,tauRFM,absoluteDiff,percentDiff)
            waveList.append(data)
    fp.close()
    basename = f.rstrip(".gnuplot") + "_errors_above_" + str(diffThreshold)
    fp = open(basename,"w")
    fp.write(basename + "(height,wavenumber,tauGRT,tauRFM," +
             "absoluteDiff,percentDiff)\n")
    for err in heightList:
        for val in err:
            for dum in val:
                fp.write(str(dum) + " ")

            fp.write("\n")
    fp.close()

    return heightList

#------------------------------------------------------------------------------
def plotSpectraForHeight(f,height,outFile):
    """
    Plot the spectra from an inputted file for a given height.
    """

    wavenumber = []
    tauGRT = []
    tauRFM = []
    fp = open(f,"r")
    for line in fp:
        cols = line.split()
        fheight = int(cols[3].strip())
        if fheight > height:
            break
        if fheight == int(height):
            wavenumber.append(float(cols[4].strip()))
            tauGRT.append(float(cols[5].strip()))
            tauRFM.append(float(cols[6].strip()))
    fp.close()
    tmpFile = "tmp.txt"
    fp = open(tmpFile,"w")
    for i in range(0,len(wavenumber)):
        fp.write(str(wavenumber[i]) + " " + str(tauGRT[i]) + " " +
                     str(tauRFM[i]) + "\n")
    fp.close()
    mols = f.split("_")
    mol = (mols[0].strip()).upper()
    title = mol + " Optical Depth at height " + str(height)
    gnuplotCommandString = ('gnuplot -e "' +
                            "reset;\n" + 
                            "set terminal png size 1000,800;\n" +
                            "set output '" + outFile + "';\n" +
                            "set title '" + title + "';\n" +
                            "set key;\n" +
                            "set tic scale 0;\n" +
                            "set xlabel 'Wavenumber (1/cm)';\n" +
                            "set ylabel 'Optical Depth';\n" +
                            "set log y;\n" +
                            "plot '" + tmpFile + "' using 1:2 lc rgb 'red'" +
                                " title 'GRT' with line, " +
                                "'" + tmpFile + "' using 1:3 lc rgb 'green'" +
                                " title 'RFM' with line;" +
                            '"')
    tmp = subprocess.Popen(gnuplotCommandString,
                           shell=True,
                           stdout=subprocess.PIPE)
    tmp.wait()
    if tmp.returncode != 0:
        print "Error: gnuplot failed."
        exit(1)
    os.remove(tmpFile)

    return

#------------------------------------------------------------------------------
def makeSpectraGif(f,startHeight,endHeight,gifName,cleanup):
    """
    Create of Gif of spectra for various heights.  Use "display gifName"
    to view the gif on a GFDL workstation.
    """

    framesList = []
    frameImages = ""
    for i in range(startHeight,endHeight):
        frameName = "frame_" + str(i) + ".png"
        framesList.append(frameName)
        plotSpectraForHeight(f,i,frameName)
        frameImages = frameName + " " + frameImages
    tmp = subprocess.Popen("convert -delay 50 -loop 0 " + frameImages + " "
                               + gifName,
                           shell=True,
                           stdout=subprocess.PIPE)
    tmp.wait()
    if tmp.returncode != 0:
        print "Error: convert failed."
        exit(1)
    if cleanup:
        for fname in framesList:
            os.remove(fname)

    return

#------------------------------------------------------------------------------
#Main part of the script.

#Parse command line arguments.
parser = optparse.OptionParser()
parser.add_option("-d",
                  "--diffs",
                  dest="calcDiffs",
                  action="store_true")
parser.add_option("-g",
                  "--gif",
                  dest="makeGifs",
                  action="store_true")
parser.add_option("-f",
                  "--finderrs",
                  dest="findErrors",
                  action="store_true")
options,args = parser.parse_args()

if (not options.calcDiffs and not options.makeGifs and not
    options.findErrors):
    options.calcDiffs = True

#Create the appropriate output for each file with a .gnuplot extentsion
#in the current directory.
files = os.listdir(os.getcwd())
patternGnuplot = re.compile(r'\.gnuplot')
for f in files:
    if patternGnuplot.search(f):
        basename = f.rstrip(".gnuplot")
        if options.calcDiffs:
            createHeatMap(f)
        if options.makeGifs:
            makeSpectraGif(f,0,48,basename+".gif",True)
        if options.findErrors:
            HighErrors = findHighErrors(f,100.0)
