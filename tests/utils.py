import os
import subprocess

def run_make(buildDir,
             target):
    """
    Run make target in the inputted build directory.
    """

    #Check input types.
    if not isinstance(buildDir,str):
        raise TypeError("inputted build directory (" + repr(buildDir) +
                            ") must be a string.\n")
    if not isinstance(target,str):
        raise TypeError("inputted make target (" + repr(target) +
                            ") must be a string.\n")

    #Store the current directory.
    pwd = os.getcwd()

    #Change to the build directory.
    os.chdir(buildDir)

    #Make sure a Makefile exists in the build directory.
    if not "Makefile" in os.listdir("."):
        raise ValueError("no Makefile exists in the build directory " +
                             buildDir + ".\n")

    #Run make target.
    makeArgs = ["make",
                target]
    tmp = subprocess.Popen(makeArgs)
    tmp.wait()
    if tmp.returncode != 0:
        raise ValueError("make " + target + " failed and returned code " +
                             str(tmp.returncode) + ".\n")

    #Change back to the directory you started in.
    os.chdir(pwd)

def copy_file(fileName,
              destDir):
    """
    Copy the inputted file to the inputted destination directory.
    """

    #Check input types.
    if not isinstance(fileName,str):
        raise TypeError("inputted file (" + repr(fileName) +
                            ") must be a string.\n")
    if not isinstance(destDir,str):
        raise TypeError("inputted destination directory (" + repr(destDir) +
                            ") must be a string.\n")

    #Copy the file.
    cpArgs = ["cp",
              fileName,
              destDir]
    tmp = subprocess.Popen(cpArgs)
    tmp.wait()
    if tmp.returncode != 0:
        raise ValueError("cp " + fileName + " " + destDir + " failed and " +
                             "returned code " + str(tmp.returncode) + ".\n")

def run_executable(executable,
                   args):
    """
    Run the inputted executable with the inputted arguments.
    """

    #Check input types."
    if not isinstance(executable,str):
        raise TypeError("inputted executable (" + repr(executable) +
                            ") must be a string.\n")
    if not isinstance(args,list) and not isinstance(args,set):
        raise TypeError("inputted arguments must be a list or set.\n")
    for arg in args:
        if not isinstance(arg,str):
            raise TypeError("inputted argument (" + repr(args) +
                                ") must be a string.\n")

    #Remove duplicates from args.
    tmpSet = set(args)

    #Run the executable.
    execArgs = []
    execArgs.append(executable)
    for arg in tmpSet:
        execArgs.append(arg)
    tmp = subprocess.Popen(execArgs)
    tmp.wait()
    if tmp.returncode != 0:
        execString = ""
        for arg in execArgs:
            execString += arg + " "
        raise ValueError(execString + " failed and returned code " +
                             str(tmp.returncode) + ".\n")

def move_file(fileName,
              destDir):
    """
    Move the inputted file to the inputted destination directory.
    """

    #Check input types.
    if not isinstance(fileName,str):
        raise TypeError("inputted file (" + repr(fileName) +
                            ") must be a string.\n")
    if not isinstance(destDir,str):
        raise TypeError("inputted destination directory (" + repr(destDir) +
                            ") must be a string.\n")

    #Copy the file.
    mvArgs = ["mv",
              fileName,
              destDir]
    tmp = subprocess.Popen(mvArgs)
    tmp.wait()
    if tmp.returncode != 0:
        raise ValueError("mv " + fileName + " " + destDir + " failed and " +
                             "returned code " + str(tmp.returncode) + ".\n")
