#!/usr/bin/env python
from os import chdir, getcwd, linesep, listdir, remove, unlink
from re import compile, escape
from subprocess import Popen, PIPE
from time import time
import errno
from extract_single_line import create_lines_specific_hitran_file
from hitran_utils import hitranDict, hitranDictKeyString
from utils import copy_file, move_file, run_executable, run_make

lineShapeDict = {"voigt"   : "VOI",
                 "lorentz" : "LOR",
                 "doppler" : "DOP"}

lineShapeDictKeyString = ""
for key in lineShapeDict:
    lineShapeDictKeyString += "\t" + key + "\n"

def create_hitbin_from(parfiles, out="mypar.bin", header="TEMPORARY FILE: MODIFY AT YOUR OWN RISK"):
    parout = out + ".par"
    # Concatenate parfiles
    with open(parout, "w") as fout:
        for f in parfiles:
            with open(f) as fin:
                fout.write(fin.read())
    # Get rid of the old bin file. sorry old bin file that you may have needeed...
    try:
        unlink(out)
    except OSError as e:
        if e.errno != errno.ENOENT:
            raise
    # Convert combined parfile to bin
    pinput = [parout,
              "", # wavenumbers, use default for now
              out,
              header]
    Popen(['./hitbin'], stdin=PIPE).stdin.write(linesep.join(pinput) + linesep)

class Layer(object):
    """
    A container to hold layer data.
    """

    def __init__(self):
        self.layerID = 0
        self.temperature = 0
        self.pressure = 0
        self.deltaz = 0
        self.rh2o = 0
        self.ro3 = 0

    def print_conds(self):
        print("")
        print("layer id:    ",self.layerID)
        print("temperature: ",self.temperature)
        print("pressure:    ",self.pressure)
        print("deltaz:      ",self.deltaz)
        print("rh2o:        ",self.rh2o)
        print("ro3:         ",self.ro3)
        print("")

def get_layer_conditions(filename):
    """
    Read in the layer conditions from a file.
    """

    layers = []
    f = open(filename,
             "r")
    next(f)
    for line in f:
        foo = Layer()
        vals = line.split()
        foo.layerID = vals[0]
        foo.temperature = vals[1]
        foo.pressure = str(float(vals[2])/100.)
        foo.deltaz = str(float(vals[3])/1000.)
        foo.rh2o = str(float(vals[4])*(29.0/18.0)*1000000.0)
        foo.ro3 = str(float(vals[5])*(29.0/48.0)*1000000.0)
        layers.append(foo)
    f.close()

    return layers

def modify_line_after_identifiers(identifier,
                                  string,
                                  input_file):
    with open(input_file,"r") as f:
        lines = f.readlines()
    with open(input_file,"w") as f:
        for i,line in enumerate(lines):
            if identifier in lines[max(i-1,0)]:
                f.write(string + "\n")
            else:
                f.write(line)

def modify_rfm_output(fileName,
                      resDir):
    """
    Rewrite RFM output in a format that is more gnuplot friendly.
    """

    f = open(fileName,"r")
    f1 = open(resDir + "/" + fileName,"w")

    for i in range(0,4):
        line = f.readline()
        f1.write(line)

    for line in f:
        vals = line.split()
        for val in vals:
            f1.write(val.strip() + "\n")

    f.close()
    f1.close()

def run_rfm(mols,
            layerFile,
            baseDir,
            lineShape,
            minFreq,
            maxFreq,
            freqRes,
            lines=[],
            forceBuild=False):
    """
    Build (if necessary) and run the RFM model.
    """

    #Check mols input.
    if not isinstance(mols,list) and not isinstance(mols,set):
        raise TypeError("the molecules input must be a list or set.\n")
    for m in mols:
        if not isinstance(m,str):
            raise TypeError("the inputted molecule (" + repr(m) +
                                ") must be a string.\n")
        tmp = (m.strip()).lower()
        if not tmp in hitranDict:
            raise ValueError("the inputted molecule (" + tmp + ") must be" +
                                 " one of:\n" + hitranDictKeyString)

    #Check layerFile input.
    if not isinstance(layerFile,str):
        raise TypeError("the inputted layer file (" + repr(layerFile) +
                            ") must be a string.\n")

    #Check baseDir input.
    if not isinstance(baseDir,str):
        raise TypeError("the inputted grtcode base directory (" +
                            repr(baseDir) + ") must be a string.\n")

    #Check lineShape input.
    if not isinstance(lineShape,str):
        raise TypeError("the inputted line shape (" + repr(lineShape) +
                            ") must be a string.\n")
    if not lineShape in lineShapeDict:
        raise ValueError("the inputted line shape (" + lineShape +
                             ") must be one of:\n" + lineShapeDictKeyString)

    #Check minFreq input.
    if not isinstance(minFreq,int):
        raise TypeError("the inputted frequency lower bound(" +
                            repr(minFreq) + ") must be an int.\n")

    #Check maxFreq input.
    if not isinstance(maxFreq,int):
        raise TypeError("the inputted frequency upper bound(" +
                            repr(maxFreq) + ") must be an int.\n")

    #Check freqRes input.
    if not isinstance(freqRes,int) and not isinstance(freqRes,float):
        raise TypeError("the inputted frequency resolution(" +
                            repr(freqRes) + ") must be an int or a float.\n")

    #Store the current directory.
    pwd = getcwd()

    #Change into the inputted GRT base directory.  Store necessary paths.
    chdir(baseDir)
    grtHitDir = getcwd() + "/run/HITFILES"
    rfmHomeDir = getcwd() + "/run/verification/rfm_layer_test"
    rfmBuildDir = rfmHomeDir + "/build"
    rfmRunDir = rfmHomeDir + "/run"
    rfmResultsDir = rfmRunDir + "/RESULTS"
    rfmHitBinDir = getcwd() + "/tests"

    #Check that the inputted layer file exists in the correct directory.
    if layerFile not in listdir(rfmRunDir):
        raise ValueError("the inputted layer file (" + layerFile +
                             ") does not exist in the directory " +
                             rfmRunDir + ".\n")

    #Check that the rfm.drv file exists in the correct directory.
    if "rfm.drv" not in listdir(rfmRunDir):
        raise ValueError("the required rfm.drv file does not exist in the " +
                             "directory " + rfmRunDir + ".\n")

    #Check that the test.atm file exists in the correct directory.
    if "test.atm" not in listdir(rfmRunDir):
        raise ValueError("the required test.atm file does not exist in the " +
                             "directory " + rfmRunDir + ".\n")

    #Remove any duplicates from the mols list.
    molecules = set(mols)

    #Check whether the necessary hitran files exist.
    chdir(grtHitDir)
    for m in molecules:
        tmp = (m.strip()).lower()
        if not hitranDict[tmp] in listdir("./"):
            raise ValueError("the hitran file (" + hitranDict[tmp] +
                                 ") does not exist in the directory " +
                                 grtHitDir + ".\n")

    #Remove any duplicates from the lines list.
    spectralLines = set(lines)

    #If a specific set of lines will be used, then generate the necessary
    #hitran file for each molecule.
    if len(lines) > 0:
        newHitranFiles = {}
        for m in molecules:
            tmp = (m.strip()).lower()
            hitranPath = grtHitDir + "/" + hitranDict[tmp]
            newHitranFiles[tmp] = create_lines_specific_hitran_file(hitranPath,
                                                                    spectralLines)
        runHitDict = newHitranFiles
    else:
        runHitDict = {}
        for key in hitranDict:
            runHitDict[key] = grtHitDir + "/" + hitranDict[key]

    #Build hitbin from hitbin.f
    chdir(rfmHitBinDir)
    run_executable("make",
                   ["FC=gfortran", "hitbin"])

    #Run hitbin to make the RFM .bin file and move it to the run directory.
    hitbinf = "mypar.bin"
    create_hitbin_from(runHitDict.values(),
                       out=hitbinf)
    move_file(hitbinf,
              rfmRunDir)

    #Build RFM and move the executable to the run directory.
    executable = "rfm.x"
    if forceBuild or executable not in listdir(rfmRunDir):

        #Run make clean and make.
        run_make(rfmBuildDir,
                 "clean")
        run_make(rfmBuildDir)

        #Copy the executable to the run directory.
        copy_file(rfmBuildDir + "/" + executable,
                  rfmRunDir)

    #Read in the layer properties from the file layer_cond.
    chdir(rfmRunDir)
    layers = get_layer_conditions(layerFile)

    #Set names of rfm input files.
    rfmDriverFile = "rfm.drv"
    atmosFile = "test.atm"
    modify_line_after_identifiers("ATM",
                                  atmosFile,
                                  rfmDriverFile)

    #Put the inputted molecules into the RFM driver input file.
    modify_line_after_identifiers("GASs",
                                  (" ".join(map(str,molecules))).upper(),
                                  rfmDriverFile)

    #Put the inputted frequency range and frequency resolution into the
    #RFM driver input file.
    modify_line_after_identifiers("SPC",
                                  str(minFreq) + " " + str(maxFreq) + " " +
                                      str(freqRes),
                                  rfmDriverFile)

    #Put the inputted line shape into the RFM driver input file.
    modify_line_after_identifiers("SHP",
                                  lineShapeDict[lineShape] + "   *",
                                  rfmDriverFile)

    #Put the created HITRAN binary file into the RFM driver input file.
    modify_line_after_identifiers("HIT",
                                  hitbinf,
                                  rfmDriverFile)

    #Set a base name for the RFM output files.  One output file will
    #be generated per layer.
    rfmOutBaseName = (("_".join(map(str,molecules))).lower() + "." +
                          lineShape.lower() + "." + str(minFreq) + "_" +
                          str(maxFreq) + "_" + str(freqRes) +
                          ".rfm.out.layer.")

    #Loop through the layers.
    timing = 0.0
    for layer in layers:

        #Append the layer id onto the end of the RFM output file name.
        rfmOutFile = rfmOutBaseName + str(layer.layerID)

        #Put the output file name in the RFM driver input file.
        modify_line_after_identifiers("OPT",
                                      rfmOutFile,
                                      rfmDriverFile)

        #Put the correct layer thickness in the RFM driver input file.
        modify_line_after_identifiers("TAN",
                                      str(layer.deltaz),
                                      rfmDriverFile)

        #Put the correct layer pressure in the RFM atmosphere input file.
        modify_line_after_identifiers("Pre",
                                      str(layer.pressure),
                                      atmosFile)

        #Put the correct layer temperature in the RFM atmosphere input file.
        modify_line_after_identifiers("tem",
                                      str(layer.temperature),
                                      atmosFile)

        #Put the correct molecular concentrations in the RFM atmosphere input
        #file.
        m = (" ".join(map(str,molecules))).lower()
        if "h2o" in m:
            modify_line_after_identifiers("h2o",
                                          str(layer.rh2o),
                                          atmosFile)
        if "o3" in m:
            modify_line_after_identifiers("o3",
                                          str(layer.ro3),
                                          atmosFile)

        #Run the executable.  Time how long the executable takes to run.
        start = time()
        run_executable(rfmRunDir + "/" + executable)
        timing += time() - start

    #Modify the RFM output files so that each optical depth value is on its
    #own line and move the output files into the RFM results directory.
    patternString = r'' + escape(rfmOutBaseName)
    patternRFM = compile(patternString)
    contents = listdir(rfmRunDir)
    for item in contents:
        if patternRFM.search(item.strip()):
            modify_rfm_output(item.strip(),
                              rfmResultsDir)
            remove(item.strip())

    #Change back to the directory you started in.
    chdir(pwd)

    return timing
