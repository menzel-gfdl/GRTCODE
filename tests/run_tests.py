#!/usr/bin/env python

import math
import optparse
import os
import re
import sys
from run_grtcode import run_grtcode
from run_rfm import run_rfm

MIN_RES = 10.
MAX_RES = 0.001
MIN_LINE_FREQ = 1
MAX_LINE_FREQ = 3000

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class testParams(object):
    """
    Container to hold test parameters.
    """

    def __init__(self,
                 configFile,
                 testName):
        """
        Initialize a testParams object by parsing a config file.
        """

        #Check for the config file in the current directory.
        if not configFile in os.listdir("./"):
            raise ValueError("config file (" + configFile + ") does not " +
                                 "exist in the current directory.")

        #Initialize members of the class.
        self.testName = testName
        self.configFile = configFile
        self.mols = []
        self.lineShape = ""
        self.lines = []
        self.lowFreq = -1
        self.highFreq = -1
        self.res = -1.0
        self.doAllLines = False

        #Open the config file and look for the test name.
        f = open(configFile,
                 "r")
        testFound = False
        molsFound = False
        lineShapeFound = False
        linesFound = False
        lowFreqFound = False
        highFreqFound = False
        resFound = False

        testRegex = r'^\s*@test\s*=\s*' + re.escape(testName) + r'\s*$'
        patternTestName = re.compile(testRegex)
        patternTestEnd = re.compile(r'^\s*@endtest\s*$')
        patternMols = re.compile(r'^\s*mols\s*=')
        patternLineShapes = re.compile(r'^\s*lineshapes\s*=')
        patternLines = re.compile(r'^\s*lines\s*=')
        patternLowFreq = re.compile(r'^\s*lowFreq\s*=')
        patternHighFreq = re.compile(r'^\s*highFreq\s*=')
        patternRes = re.compile(r'^\s*resolution\s*=')
        patternAlpha = re.compile(r'^\s*([A-Z]|[a-z])+\s*$')
        patternAlphaNum = re.compile(r'^\s*([A-Z]|[a-z]|[0-9])+(\s|[A-Z]|[a-z]|[0-9])*$')
        patternInt = re.compile(r'^\s*([0-9])+\s*$')
        patternFloat = re.compile(r'^\s*([0-9])+\.([0-9])*\s*$')
        patternFloats = re.compile(r'^\s*([0-9])+\.([0-9])*(\s|(\s([0-9])+\.([0-9])*\s))*$')

        for line in f:
            if patternTestName.search(line):
                #Test is found in the file.
                testFound = True
                sys.stdout.write("\nTest '" + testName + "' found in file '" +
                                     configFile + "'.\n")

            #Parse out the test parameters.
            if testFound:
                if patternTestEnd.search(line):
                    #End of test found, break out of loop.
                    break

                elif patternMols.search(line):
                    #Make sure at least one valid molecule exists.
                    if not molsFound:
                        molsFound = True
                    else:
                        raise ValueError("More than one mols field exists" +
                                             " in this test.\n")
                    errMesg1 = ("Error on line:\n" + line +
                                    "\nMolecule must be one " +
                                    "of:\n\th2o\n\tco2\n\to3" +
                                    "\n\tn2o\n\tco\n\tch4\n\t" +
                                    "o2\n")
                    tmpString = patternMols.split(line)
                    if not patternAlphaNum.search(tmpString[1]):
                        raise ValueError(errMesg1)
                    tmp = tmpString[1].split()
                    for m in tmp:
                        m = m.lower()
                        if (m != "h2o" and m != "co2" and m != "o3" and
                                m != "n2o" and m != "co" and m != "ch4" and
                                m != "o2"):
                            raise ValueError(errMesg1)
                        else:
                            self.mols.append(m)
                    self.mols = set(self.mols)

                elif patternLineShapes.search(line):
                    #Make sure a valid line shape exists.
                    if not lineShapeFound:
                        lineShapeFound = True
                    else:
                        raise ValueError("More than one lineshape field" +
                                             " exists in this test.\n")
                    errMesg1 = ("Error on line:\n" + line +
                                    "\nLine shape must be one of:" +
                                    "\n\tvoigt\n\tvoigt_ida\n\t" +
                                    "lorentz\n\tdoppler\n")
                    tmpString = patternLineShapes.split(line)
                    if not patternAlpha.search(tmpString[1]):
                        raise ValueError(errMesg1)
                    tmp = tmpString[1].split()
                    tmp2 = tmp[0].lower()
                    if (tmp2 != "voigt" and tmp2 != "voigt_ida" and
                            tmp2 != "lorentz" and tmp2 != "doppler"):
                        raise ValueError(errMesg1)
                    else:
                        self.lineShape = tmp2

                elif patternLines.search(line):
                    #Make sure that valid line frequencies exist.
                    if not linesFound:
                        linesFound = True
                    else:
                        raise ValueError("More than one lines field" +
                                             " exists in this test.\n")
                    errMesg1 = ("Error on line:\n" + line +
                                    "\nLines must be float " +
                                    "values or 'all'.\n")
                    tmpString = patternLines.split(line)
                    if patternAlpha.search(tmpString[1]):
                        tmp = tmpString[1].split()
                        tmp2 = tmp[0].lower()
                        if tmp2 != "all":
                            raise ValueError(errMesg1)
                        else:
                            self.lines.append(tmp2)
                            self.doAllLines = True
                    elif patternFloats.search(tmpString[1]):
                        tmp = tmpString[1].split()
                        for m in tmp:
                            self.lines.append(float(m))
                        self.lines = set(self.lines)
                    else:
                        raise ValueError(errMesg1)

                elif patternLowFreq.search(line):
                    #Make sure that a valid lower frequency bound exists.
                    if not lowFreqFound:
                        lowFreqFound = True
                    else:
                        raise ValueError("More than one lowFreq field" +
                                             " exists in this test.\n")
                    errMesg1 = ("Error on line:\n" + line +
                                    "\nLower frequency bound must be one" +
                                    " integer.\n")
                    tmpString = patternLowFreq.split(line)
                    if not patternInt.search(tmpString[1]):
                        raise ValueError(errMesg1)
                    tmp = tmpString[1].split()
                    self.lowFreq = int(tmp[0])

                elif patternHighFreq.search(line):
                    #Make sure that a valid higher frequency bound exists.
                    if not highFreqFound:
                        highFreqFound = True
                    else:
                        raise ValueError("More than one highFreq field" +
                                             " exists in this test.\n")
                    errMesg1 = ("Error on line:\n" + line +
                                    "\nHigher frequency bound must be one" +
                                    " integer.\n")
                    tmpString = patternHighFreq.split(line)
                    if not patternInt.search(tmpString[1]):
                        raise ValueError(errMesg1)
                    tmp = tmpString[1].split()
                    self.highFreq = int(tmp[0])

                elif patternRes.search(line):
                    #Make sure that a valid frequency resolution exists.
                    if not resFound:
                        resFound = True
                    else:
                        raise ValueError("More than one res field" +
                                             " exists in this test.\n")
                    errMesg1 = ("Error on line:\n" + line +
                                    "\nFrequency resolution must be one" +
                                    " float.\n")
                    errMesg2 = ("Error on line:\n" + line +
                                    "\nFrequency resolution must be <= " +
                                    str(MIN_RES) + " and >= " +
                                    str(MAX_RES) + ".\n")
                    tmpString = patternRes.split(line)
                    if not patternFloat.search(tmpString[1]):
                        raise ValueError(errMesg1)
                    tmp = tmpString[1].split()
                    self.res = float(tmp[0])
                    if self.res > MIN_RES or self.res < MAX_RES:
                        raise ValueError(errMesg2)
        f.close()

        #Make sure that the test was found in the config file.
        if not testFound:
            raise ValueError("test (" + testName + ") does not exist in " +
                                 " the config file (" + configFile + ").\n")

        #Make sure that all the necessary fields were found.
        if not molsFound:
            raise ValueError("mols field does not exist in the test.\n")
        if not lineShapeFound:
            raise ValueError("lineshape field does not exist in the test.\n")
        if not linesFound:
            raise ValueError("lines field does not exist in the test.\n")
        if not lowFreqFound:
            raise ValueError("lowFreq field does not exist in the test.\n")
        if not highFreqFound:
            raise ValueError("highFreq field does not exist in the test.\n")
        if not resFound:
            raise ValueError("res field does not exist in the test.\n")

        #Make sure the inputted frequency ranges are valid.
        maxLine = -1
        minLine = 100000000
        if not self.doAllLines:
            for line in self.lines:
                if line < float(MIN_LINE_FREQ) or line > float(MAX_LINE_FREQ):
                    raise ValueError("inputted line (" + str(line) +
                                         ") must be >= " + str(MIN_LINE_FREQ)
                                         + " and <= " + str(MAX_LINE_FREQ) +
                                         ".\n")
                if math.ceil(line) > maxLine:
                    maxLine = math.ceil(line)
                if math.floor(line) < minLine:
                    minLine = math.floor(line)
        else:
            maxLine = MAX_LINE_FREQ
            minLine = MIN_LINE_FREQ
        if maxLine < minLine:
            raise ValueError("maxLine (" + str(maxLine) + ") cannot be" +
                                 " < minLine (" + str(minLine) + ").\n")
        if maxLine > self.highFreq:
            raise ValueError("frequency upper bound (" + str(self.highFreq) +
                                 ") must be >= the highest line (" +
                                 str(maxLine) + ").\n")
        if minLine < self.lowFreq:
            raise ValueError("frequency lower bound (" + str(self.lowFreq) +
                                 ") must be <= the lowest line (" +
                                 str(minLine) + ").\n")
        if self.highFreq < self.lowFreq:
            raise ValueError("highFreq (" + str(self.highFreq) + ") cannot"
                             " be < lowFreq (" + str(self.lowFreq) + ").\n")

    def show(self):
        """
        Print out the members of self.
        """

        print self.testName
        print self.configFile
        print self.mols
        print self.lineShape
        print self.lines
        print self.lowFreq
        print self.highFreq
        print self.res

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#def run_rfm(testObj):

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if __name__ == "__main__":

    #Parser command line arguments.
    parser = optparse.OptionParser()
    parser.add_option("-c",
                      "--configfile",
                      dest="configFile",
                      action="store",
                      type="string")
    parser.add_option("-t",
                      "--testname",
                      dest="testName",
                      action="store",
                      type="string")
    options,args = parser.parse_args()

    #Check inputs.
    if options.configFile == None:
        sys.stdout.write("\nUsing default config file (tests.config).\n")
        options.configFile = "tests.config"
    if options.testName == None:
        sys.stdout.write("\nPreforming @test = standard.\n")
        options.testName = "standard"

    #Parse the config file to get the test parameters.
    testObject = testParams(options.configFile,
                            options.testName)
    testObject.show()

    #Run the test using GRTcode.
    if testObject.doAllLines:
        tLines = []
    else:
        tLines = testObject.lines
    grt_output_file, grt_timing = run_grtcode(testObject.mols,
                                              "smallSubset_2t.nc",
                                              "../",
                                              testObject.lineShape,
                                              testObject.lowFreq,
                                              testObject.highFreq,
                                              testObject.res,
                                              lines=tLines,
                                              forceBuild=False)

    #Run the test using rfm.
    rfm_timing = run_rfm(testObject.mols,
                         "layer_cond",
                         "../",
                         testObject.lineShape,
                         testObject.lowFreq,
                         testObject.highFreq,
                         testObject.res,
                         lines=tLines,
                         forceBuild=False)

    #Write out timings to stdout.
    sys.stdout.write("\nTimings: \nGRTcode runtime (s): " + str(grt_timing) +
                     "\nRFM runtime (s):     " +  str(rfm_timing) + "\n")

