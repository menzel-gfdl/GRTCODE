from os import chdir, getcwd, listdir
from time import time
from extract_single_line import create_lines_specific_hitran_file
from hitran_utils import hitranDict, hitranDictKeyString
from utils import copy_file, move_file, run_executable, run_make

#Dictionary used for running the grtcode executable.
gfdl_ppmvDict = {"h2o" : "-1a",
                 "co2" : "-2400",
                 "o3"  : "-3a",
                 "n2o" : "-40.32",
                 "co"  : "-50.001",
                 "ch4" : "-61.7",
                 "o2"  : "-7200000"}

rfmip_ppmvDict = {"h2o" : "-1a",
                  "co2" : "-2a",
                  "o3"  : "-3a",
                  "n2o" : "-4a",
                  "co"  : "-5a",
                  "ch4" : "-6a",
                  "o2"  : "-7a"}

atmosFileTypeList = ["gfdl",
                     "rfmip"]
atmosFileTypeListString = ""
for item in atmosFileTypeList:
    atmosFileTypeListString += "\t" + item + "\n"

#Dictionary used for running the correct grtcode executable.
grtExecDict = {"voigt"     : "grtcode.x",
               "voigt_ida" : "grtcodeIdaVoigt.x",
               "lorentz"   : "grtcodeLorentz.x",
               "doppler"   : "grtcodeGauss.x"}

#String containing the names of all valid line shape strings.
grtExecDictKeyString = ""
for key in grtExecDict:
    grtExecDictKeyString += "\t" + key + "\n"

archTypeList = ["gpu",
                "cpu",
                "cpu_openmp"]
archTypeListString = ""
for item in archTypeList:
    archTypeListString += "\t" + item + "\n"

def run_grtcode(architecture,
                mols,
                atmosFile,
                atmosFileType,
                baseDir,
                lineShape,
                minFreq,
                maxFreq,
                freqRes,
                lines,
                skipBuild=False):
    """
    Build and run grtcode.  Return the path of the output
    file and the time it took to run the executable.
    """

    #Check architecture input.
    if not isinstance(architecture,str):
        raise TypeError("the inputted architecture (" + repr(architecture) +
                            ") must be a string.\n")
    if architecture.lower() not in archTypeList:
        raise ValueError("the inputted architecture (" + architecture +
                             ") must be one of:\n" + archTypeListString)

    #Check mols input.
    if not isinstance(mols,list) and not isinstance(mols,set):
        raise TypeError("the molecules input must be a list or set.\n")
    for m in mols:
        if not isinstance(m,str):
            raise TypeError("the inputted molecule (" + repr(m) +
                                ") must be a string.\n")
        tmp = (m.strip()).lower()
        if tmp not in hitranDict:
            raise ValueError("the inputted molecule (" + tmp + ") must be" +
                                 " one of:\n" + hitranDictKeyString)

    #Check atmosFile input.
    if not isinstance(atmosFile,str):
        raise TypeError("the inputted atmosphere file (" + repr(atmosFile) +
                            ") must be a string.\n")

    #Check atmosFileType input.
    if not isinstance(atmosFileType,str):
        raise TypeError("the inputted atmosphere file (" + repr(atmosFileType)
                            + ") must be a string.\n")
    tmp = (atmosFileType.strip()).lower()
    if tmp not in atmosFileTypeList:
        raise ValueError("the inputted atmosphere file type (" + tmp +
                             ") must be one of:\n" +
                             atmosFileTypeListString)

    #Check baseDir input.
    if not isinstance(baseDir,str):
        raise TypeError("the inputted grtcode base directory (" +
                            repr(baseDir) + ") must be a string.\n")

    #Check lineShape input.
    if not isinstance(lineShape,str):
        raise TypeError("the inputted line shape (" + repr(lineShape) +
                            ") must be a string.\n")
    if lineShape not in grtExecDict:
        raise ValueError("the inputted line shape (" + lineShape +
                             ") must be one of:\n" + grtExecDictKeyString)

    #Check minFreq input.
    if isinstance(minFreq,str):
        try:
            minFreq = int(minFreq.strip())
        except:
            raise ValueError("the inputted frequency lower bound (" +
                                 repr(minFreq) + ") cannot be converted" +
                                 " to an int.\n")
    else:
        if not isinstance(minFreq,int):
            raise TypeError("the inputted frequency lower bound (" +
                                repr(minFreq) + ") must be an int.\n")

    #Check maxFreq input.
    if isinstance(maxFreq,str):
        try:
            maxFreq = int(maxFreq.strip())
        except:
            raise ValueError("the inputted frequency upper bound (" +
                                 repr(maxFreq) + ") cannot be converted" +
                                 " to an int.\n")
    else:
        if not isinstance(maxFreq,int):
            raise TypeError("the inputted frequency upper bound(" +
                                repr(maxFreq) + ") must be an int.\n")

    #Check freqRes input.
    if isinstance(freqRes,str):
        try:
            freqRes = float(freqRes.strip())
        except:
            raise ValueError("the inputted frequency resolution (" +
                                 repr(freqRes) + ") cannot be converted" +
                                 " to a float.")
    else:
        if not isinstance(freqRes,int) and not isinstance(freqRes,float):
            raise TypeError("the inputted frequency resolution (" +
                                repr(freqRes) + ") must be an int or a " +
                                "float.\n")

    #Check lines input.
    use_all_lines = False
    spectralLines = []
    if not isinstance(lines,list) and not isinstance(lines,set):
        raise TypeError("the inputted lines must be a list or a set.\n")
    if len(lines) == 1 and lines[0].strip() == "all":
        use_all_lines = True
    else:
        for line in lines:
            if isinstance(line,str):
                try:
                    spectralLines.append(float(line.strip()))
                except:
                    raise TypeError("the inputted line frequency (" +
                                        repr(line) + ") cannot be converted" +
                                        " to a float.")
            elif not isinstance(lines,int) and not isinstance(lines,float):
                raise TypeError("the inputted line frequency (" +
                                    repr(line) + ") must be an int or a " +
                                    "float.\n")

    #Store the current directory.
    pwd = getcwd()

    #Change into the inputted GRT base directory.  Store necessary paths.
    chdir(baseDir)
    grtHomeDir = getcwd()
    grtBuildDir = grtHomeDir + "/build"
    grtRunDir = grtHomeDir + "/run"
    grtInputDir = grtRunDir + "/INPUT"
    grtHitDir = grtRunDir + "/HITFILES"
    grtResultsDir = grtRunDir + "/RESULTS"

    #Check that the inputted atmosphere file exists in the correct directory.
    if atmosFile not in listdir(grtInputDir):
        raise ValueError("the inputted atmosphere file (" + atmosFile +
                             ") does not exist in the directory " +
                             grtInputDir + ".\n")

    #Remove any duplicates from the mols list.
    molecules = set(mols)

    #Check whether the necessary hitran files exist.
    for m in molecules:
        tmp = (m.strip()).lower()
        if hitranDict[tmp] not in listdir(grtHitDir):
            raise ValueError("the hitran file (" + hitranDict[tmp] +
                                 ") does not exist in the directory " +
                                 grtHitDir + ".\n")

    #If a specific set of lines will be used, then generate the necessary
    #hitran file for each molecule.
    if not use_all_lines:
        spectralLines = set(spectralLines)
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
    if not skipBuild:

        #Run make clean and make all_grtcode.
        if architecture.lower() == "gpu":
            makefile = "Makefile"
        else:
            makefile = "Makefile.gnu"

        if architecture.lower() == "cpu_openmp":
            opts = "OPENMP=on"
        else:
            opts = ""

        run_make(grtBuildDir,
                 makefile,
                 "clean")

        run_make(grtBuildDir,
                 makefile,
                 "all_grtcode " + opts)

        #Copy the executable to the run directory.
        copy_file(grtBuildDir + "/" + executable,
                  grtRunDir)

    #Run the executable.  Time how long the executable takes to run.
    grtOutputFile = "foo"
    tmp = (atmosFileType.strip()).lower()
    if tmp == "gfdl":
        ppmvDict = gfdl_ppmvDict
        input_file_format = "gfdl"
    else:
        ppmvDict = rfmip_ppmvDict
        input_file_format = "rfmip"
    args = ["-a" + grtInputDir + "/" + atmosFile,
            "-o" + grtOutputFile,
            "-w" + str(minFreq),
            "-W" + str(maxFreq),
            "-r" + str(freqRes),
            "-f" + str(input_file_format)]
    if architecture.lower() != "gpu":
        args.append("-h")
    for m in molecules:
        tmp = (m.strip()).lower()
        args.append(ppmvDict[tmp])
        args.append(runHitDict[tmp])
    start = time()
    run_executable(grtRunDir + "/" + executable,
                   args)
    timing = time() - start

    #Move the output into the results directory.
    move_file(grtOutputFile,
              grtResultsDir)

    #Change back to the directory you started in.
    chdir(pwd)

    return (grtResultsDir + "/" + grtOutputFile),timing
