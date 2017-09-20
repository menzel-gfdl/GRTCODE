import os
import subprocess

def run_make(buildDir,
             makefile="Makefile",
             target=[]):
    """
    Run make target in the inputted build directory.
    """

    #Check input types.
    if not isinstance(buildDir,str):
        raise TypeError("inputted build directory (" + repr(buildDir) +
                            ") must be a string.\n")
    if isinstance(target,str):
        target = target.split()

    #Store the current directory.
    pwd = os.getcwd()

    #Make sure a Makefile exists in the build directory.
    if not makefile in os.listdir(buildDir):
        raise ValueError(makefile + "does not exist in the build directory " +
                             buildDir + ".\n")

    #Change to the build directory.
    os.chdir(buildDir)

    #Run make target.
    makeArgs = ["make"] + ["-j12"] + ["-f"] + [makefile] + target
    tmp = subprocess.Popen(makeArgs)
    tmp.wait()
    if tmp.returncode != 0:
        raise ValueError("make -f " + makefile + " " + 
                             " ".join(map(str,target)) +
                             " failed and returned code " +
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

def run_executable(args):
    """
    Run the inputted executable with the inputted arguments.
    """

    #Check input types."
    if not isinstance(args,list):
        raise TypeError("inputted command (" + repr(args) +
                            ") must be a list.\n")

    for arg in args:
        if not isinstance(arg,str):
            raise TypeError("inputted argument (" + repr(arg) +
                            ") must be a string.\n")

    #Run the executable.
    tmp = subprocess.Popen(args)
    tmp.wait()
    if tmp.returncode != 0:
        execString = " ".join(map(str,args))
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
