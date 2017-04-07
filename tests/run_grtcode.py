import os
from extract_single_line import create_lines_specific_hitran_file
from hitran_utils import hitranDict, hitranDictKeyString
from utils import copy_file, move_file, run_executable, run_make

ppmvDict = {"h2o" : "-1a",
            "co2" : "-2400",
            "o3"  : "-3a",
            "n2o" : "-40.32",
            "co"  : "-50.001",
            "ch4" : "-61.7",
            "o2"  : "-7200000"}

grtExecDict = {"voigt"     : "grtcode.x",
               "voigt_ida" : "grtcodeIdaVoigt.x",
               "lorentz"   : "grtcodeLorentz.x",
               "doppler"   : "grtcodeGauss.x"}

grtExecDictKeyString = ""
for key in grtExecDict:
    grtExecDictKeyString += "\t" + key + "\n"

def run_grtcode(mols,
                atmosFile,
                baseDir,
                lineShape,
                minFreq,
                maxFreq,
                freqRes,
                lines=[],
                forceBuild=False):
    """
    Build (if necessary) and run grtcode.
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

    #Check atmosFile input.
    if not isinstance(atmosFile,str):
        raise TypeError("the inputted atmosphere file (" + repr(atmosFile) +
                            ") must be a string.\n")

    #Check baseDir input.
    if not isinstance(baseDir,str):
        raise TypeError("the inputted grtcode base directory (" +
                            repr(baseDir) + ") must be a string.\n")

    #Check lineShape input.
    if not isinstance(lineShape,str):
        raise TypeError("the inputted line shape (" + repr(lineShape) +
                            ") must be a string.\n")
    if not lineShape in grtExecDict:
        raise ValueError("the inputted line shape (" + lineShape +
                             ") must be one of:\n" + grtExecDictKeyString)

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
    grtHomeDir = os.getcwd()
    grtBuildDir = grtHomeDir + "/build"
    grtRunDir = grtHomeDir + "/run"
    grtInputDir = grtRunDir + "/INPUT"
    grtHitDir = grtRunDir + "/HITFILES"
    grtResultsDir = grtRunDir + "/RESULTS"

    #Check that the inputted atmosphere file exists in the correct directory.
    if not atmosFile in os.listdir(grtInputDir):
        raise ValueError("the inputted atmosphere file (" + atmosFile +
                             ") does not exist in the directory " +
                             grtInputDir + ".\n")

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
    if lines:
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

    #Build the executable if necessary.
    executable = grtExecDict[lineShape]
    if forceBuild or not executable in os.listdir(grtRunDir):

        #Change to the build directory.
        os.chdir(grtBuildDir)

        #Run make clean and make all_grtcode.
        run_make(".",
                 "clean")
        run_make(".",
                 "all_grtcode")

        #Copy the executable to the run directory.
        copy_file(executable,
                  grtRunDir)

    #Run the executable.
    grtOutputFile = "foo"
    args = ["-a" + grtInputDir + "/" + atmosFile,
            "-o" + grtOutputFile,
            "-w" + str(minFreq),
            "-W" + str(maxFreq),
            "-r" + str(freqRes)]
    for m in molecules:
        tmp = (m.strip()).lower()
        args.append(ppmvDict[tmp])
        args.append(runHitDict[tmp])
    run_executable(grtRunDir + "/" + executable,
                   args)

    #Move the output into the results directory.
    move_file(grtOutputFile,
              grtResultsDir)

    #Change back to the directory you started in.
    os.chdir(pwd)

    return grtResultsDir + grtOutputFile
