#!/usr/bin/env python
import os
import re
import subprocess
import sys
import time
from hitran_utils import hitranDict, hitranDictKeyString

lineShapeDict = {"voigt"   : "VOI",
                 "lorentz" : "LOR",
                 "doppler" : "DOP"}

lineShapeDictKeyString = ""
for key in lineShapeDict:
    lineShapeDictKeyString += "\t" + key + "\n"

def run_rfm(mols,
            layerFile,
            baseDir,
            lineShape,
            minFreq,
            maxFreq,
            freqRes,
            lines=[],
            forceBuild=False)
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
    pwd = os.getcwd()

    #Change into the inputted GRT base directory.  Store necessary paths.
    os.chdir(baseDir)
    grtHitDir = os.getcwd() + "/run/HITFILES"
    rfmHomeDir = os.getcwd() + "/run/verification/rfm_layer_test"
    rfmBuildDir = rfmHomeDir + "/build"
    rfmRunDir = rfmHomeDir + "/run"

    #Check that the inputted layer file exists in the correct directory.
    if not layerFile in os.listdir(rfmRunDir):
        raise ValueError("the inputted layer file (" + layerFile +
                             ") does not exist in the directory " +
                             rfmRunDir + ".\n")

    #Check that the rfm.drv file exists in the correct directory.
    if not "rfm.drv" in os.listdir(rfmRunDir):
        raise ValueError("the required rfm.drv file does not exist in the " +
                             "directory " + rfmRunDir + ".\n")

    #Check that the test.atm file exists in the correct directory.
    if not "test.atm" in os.listdir(rfmRunDir):
        raise ValueError("the required test.atm file does not exist in the " +
                             "directory " + rfmRunDir + ".\n")

    #Remove any duplicates from the mols list.
    molecules = set(mols)

    #Check whether the necessary hitran files exist.
    os.chdir(grtHitDir)
    for m in molecules:
        tmp = (m.strip()).lower()
        if not hitranDict[tmp] in os.listdir("./"):
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

    #Create the necessary hitran binary file, as required by the RFM model.




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def perform_sed(sedExp,
                sfile):
    """
    Use sed to modify a file.
    """

    sedString = "sed -i " + sedExp + " " + sfile
    tmp = subprocess.Popen(sedString,
                           shell=True)
    tmp.wait()
    if tmp.returncode != 0:
        print("Error: sed failed with expression " + sedExp + " and file " +
                  sfile + ".\n")
        sys.exit(1)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def run_rfm_model():
    """
    Run the RFM model.
    """

    tmp = subprocess.Popen(["./rfm.x"])
    tmp.wait()
    if tmp.returncode != 0:
        print("Error: rfm model failed.\n")
        sys.exit(1)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def update_hitran_symlink(fileName):
    """
    Create/point a symlink called HITRAN2012.bin to the inputted file.
    """

    symString = "ln -sf " + fileName + " HITRAN2012.bin"
    tmp = subprocess.Popen(symString,
                           shell=True)
    tmp.wait()
    if tmp.returncode != 0:
        print("Error: ln failed with file " + fileName + ".\n")
        sys.exit(1)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def main(molecule,
         lineshape,
         hitranFile,
         w,
         W,
         res):
    """
    Execute this routine if the script is not imported.
    """

    #Initialize variables.
    atmosFile = "test.atm"
    rfmDriverFile = "rfm.drv"
    timings = []

    #Create a symlink to the inputted HITRAN file.
    update_hitran_symlink(hitranFile)

    #Set rfm to run only the inputted molecule.
    sedString = "'/GASs/{n;s/.*/" + molecule.upper() + "/}'"
    perform_sed(sedString,
                rfmDriverFile)

    #Set rfm to run the inputted frequency range and frequency resolution.
    freqRangeString = str(w) + " " + str(W) + " " + str(res)
    sedString = "'/SPC/{n;s/.*/" + freqRangeString + "/}'"
    perform_sed(sedString,
                rfmDriverFile)

    #Read in the layer properties from the file layer_cond.
    layers = get_layer_conditions("layer_cond")

    #For each layer, make the appropriate changes to the rfm input
    #file and run rfm.
    rfmOutBaseName = (molecule + "." + lineshape + "." + str(w) + "_" + 
                          str(W) + "_" + str(res) + ".rfm.out.layer.")
    for layer in layers:
        rfmOutFile = rfmOutBaseName + layer.layerID
        sedString = "'/Pre/{n;s/.*/" + layer.pressure + "/}'"
        perform_sed(sedString,
                    atmosFile)
        sedString = "'/tem/{n;s/.*/" + layer.temperature + "/}'"
        perform_sed(sedString,
                    atmosFile)
        if molecule == "o3":
            sedString = "'/o3/{n;s/.*/" + layer.ro3 + "/}'"
        elif molecule == "h2o":
            sedString = "'/h2o/{n;s/.*/" + layer.rh2o + "/}'"
        perform_sed(sedString,
                    atmosFile)
        sedString = "'/TAN/{n;s/.*/" + layer.deltaz + "/}'"
        perform_sed(sedString,
                    rfmDriverFile)
        sedString = "'/OPT/{n;s/.*/" + rfmOutFile + "/}'"
        perform_sed(sedString,
                    rfmDriverFile)
        if lineshape == "voigt":
            sedString = "'/SHP/{n;s/.*/VOI   \*/}'"
        elif lineshape == "lorentz":
            sedString = "'/SHP/{n;s/.*/LOR   \*/}'"
        elif lineshape == "doppler":
            sedString = "'/SHP/{n;s/.*/DOP   \*/}'"
        perform_sed(sedString,
                    rfmDriverFile)
        start = time.time()
        run_rfm_model()
        timings.append(time.time()-start)

    #Move all of the outputed files into the appropriate directory.
    patternString = r'' + re.escape(rfmOutBaseName)
    patternRFM = re.compile(patternString)
    contents = os.listdir("./")
    for item in contents:
        if patternRFM.search(item.strip()):
            modify_rfm_output(item.strip(),
                              "./RESULTS")
            os.remove(item.strip())

    #Write out the maximum, minimum, average, and total timings for the
    #rfm runs.
    maxTime = 0.0
    minTime = 1.e8
    totalTime = 0.0
    averageTime = 0.0
    for t in timings:
        if t > maxTime:
            maxTime = t
        if t < minTime:
            minTime = t
        totalTime += t
    averageTime = totalTime/len(timings)
    timingFileName = (molecule + "." + lineshape + "." + str(w) + "_" +
                          str(W) + "_" + str(res) + ".rfm.out.timing_summary")
    f = open(timingFileName,
             "w")
    f.write("Timing Summary:\n")
    f.write("\nmolecule  = " + molecule)
    f.write("\nlineshape = " + lineshape)
    f.write("\n")
    f.write("\nMax time (s):     " + str(maxTime))
    f.write("\nMin time (s):     " + str(minTime))
    f.write("\nAverage time (s): " + str(averageTime))
    f.write("\nTotal time (s):   " + str(totalTime))
    f.write("\n")
    f.close()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if __name__ == "__main__":

    #Parser command line arguments.
    parser = optparse.OptionParser()
    parser.add_option("-m",
                      "--molecule",
                      dest="molecule",
                      action="store",
                      type="string")
    parser.add_option("-l",
                      "--lineshape",
                      dest="lineshape",
                      action="store",
                      type="string")
    parser.add_option("-H",
                      "--hitran",
                      dest="hitranFile",
                      action="store",
                      type="string")
    parser.add_option("-w",
                      "--lowfreq",
                      dest="w",
                      action="store",
                      type="int")
    parser.add_option("-W",
                      "--highfreq",
                      dest="W",
                      action="store",
                      type="int")
    parser.add_option("-r",
                      "--resolution",
                      dest="res",
                      action="store",
                      type="float")
    options,args = parser.parse_args()

    #Set valid ranges for the inputs.
    validMols = []
    validMols.append("h2o")
    validMols.append("o3")
    validLineShapes = []
    validLineShapes.append("voigt")
    validLineShapes.append("lorentz")
    validLineShapes.append("doppler")
    minFreq = 1
    maxFreq = 1800
    minRes = 10
    maxRes = 0.0001

    #Set usage messages.
    molOpts = "\nMolecules:\n"
    for mol in validMols:
        molOpts += "\t" + mol + "\n"
    lineOpts = "\nLine Shapes:\n"
    for line in validLineShapes:
        lineOpts += "\t" + line + "\n"
    wOpts = ("\nFrequency range:\n\t" + str(minFreq)  + " - " + str(maxFreq) +
                 ", inclusive\n")
    resOpts = ("\nResolution range:\n\t" + str(minRes) + " - " + str(maxRes) +
                   "\n")
    usageMesg = ("\n\nUsage: ./run_rfm.py -m<molecule> -l<lineshape> " +
                     "-H<hitranfile> -w<lowestfreq> -W<highestfreq> " +
                     "-r<resolution>\n" + molOpts + lineOpts + wOpts +
                     resOpts)

    #Check inputs.
    if options.molecule == None:
        sys.stderr.write("\nError: missing molecule input." + usageMesg)
        sys.exit(1)
    elif options.molecule not in validMols:
        sys.stderr.write("\nError: invalid molecule input (" +
                             options.molecule + ")." + usageMesg)
        sys.exit(1)
    if options.lineshape == None:
        sys.stderr.write("\nError: missing line shape input." + usageMesg)
        sys.exit(1)
    elif options.lineshape not in validLineShapes:
        sys.stderr.write("\nError: invalid line shape input (" +
                             options.lineshape + ")." + usageMesg)
        sys.exit(1)
    if options.hitranFile == None:
        sys.stderr.write("\nError: missing hitran file." + usageMesg)
        sys.exit(1)
    elif options.hitranFile not in os.listdir("./"):
        sys.stderr.write("\nError: hitran file (" + options.hitranFile + 
                             ") does not exist in the "
                             "current directory." + usageMesg)
        sys.exit(1)
    if options.w == None:
        sys.stderr.write("\nError: mising lower frequency bound."
                             + usageMesg)
        sys.exit(1)
    elif options.w < minFreq or options.w > maxFreq:
        sys.stderr.write("\nError: invalid lower frequency bound (" +
                             str(options.w) + ")." + usageMesg)
        sys.exit(1)
    if options.W == None:
        sys.stderr.write("\nError: missing upper frequency bound."
                             + usageMesg)
        sys.exit(1)
    elif options.W < minFreq or options.W > maxFreq:
        sys.stderr.write("\nError: invalid upper frequency bound (" +
                             str(options.W) + ")." + usageMesg)
        sys.exit(1)
    elif options.W < options.w:
        sys.stderr.write("\nError: upper frequency bound (" + str(options.W) +
                             ") cannot be less than the lower frequency " +
                             "bound (" + str(options.w) + ")." + usageMesg)
        sys.exit(1)
    if options.res == None:
        sys.stderr.write("\nError: missing frequency resolution."
                             + usageMesg)
        sys.exit(1)
    elif options.res > minRes or options.res < maxRes:
        sys.stderr.write("\nError: invalid frequency resolution (" +
                             str(options.res) + ")." + usageMesg)
        sys.exit(1)

    main(options.molecule,
         options.lineshape,
         options.hitranFile,
         options.w,
         options.W,
         options.res)
