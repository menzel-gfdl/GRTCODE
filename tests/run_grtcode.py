from get_test_params import testParams
from os import chdir, getcwd, path
from sys import stdout
from time import time
from utils import copy_file, move_file, run_executable, run_make

#Dictionary containing required hitran files.
hitran_files = {"h2o" : "01_hit12.par",
                "co2" : "02_hit12.par",
                "o3"  : "03_hit12.par",
                "n2o" : "04_hit08.par",
                "co"  : "05_hit12.par",
                "ch4" : "06_hit12.par",
                "o2"  : "07_hit12.par"}

#Dictionary used for running the correct grtcode executable.
exec_name = {"voigt" : "grtcode.x",
             "voigt_ida" : "grtcodeIdaVoigt.x",
             "lorentz" : "grtcodeLorentz.x",
             "doppler" : "grtcodeGauss.x"}

#Dictionary used for running make on the correct makefile.
makefile_name = {"gpu_devbox_gpu" : "Makefile",
                 "theta_cpu_openmp" : "Makefile.theta",
                 "gaea.c3_cpu_openmp" : "Makefile.gaea",
                 "gaea.c4_cpu_openmp" : "Makefile.gaea"}

#Dictionaries used for running the grtcode executable.
gfdl_ppmv_flag = {"h2o" : "-1a",
                  "co2" : "-2400",
                  "o3" : "-3a",
                  "n2o" : "-40.32",
                  "co" : "-50.001",
                  "ch4" : "-61.7",
                  "o2" : "-7200000"}

rfmip_ppmv_flag = {"h2o" : "-1a",
                   "co2" : "-2a",
                   "o3" : "-3a",
                   "n2o" : "-4a",
                   "co" : "-5a",
                   "ch4" : "-6a",
                   "o2" : "-7a"}

def run_grtcode(params,
                base_dir,
                skip_build=False):
    """
    Build and run grtcode.  Return the path of the output
    file and the time it took to run the executable.
    """

    #Store the expected paths.
    base_dir = path.abspath(base_dir)
    build_dir = base_dir + "/build"
    run_dir = base_dir + "/run"
    input_dir = run_dir + "/INPUT"
    hitran_dir = run_dir + "/HITFILES"
    results_dir = run_dir + "/RESULTS"

    #Check that the inputted atmosphere file exists in the correct directory.
    atmos_file = input_dir + "/" + params.params_dict["atmos_input_file"].strip()
    if not path.isfile(atmos_file):
        raise ValueError("the inputted atmosphere file (" + atmos_file +
                         ") does not exist.\n")

    #Check whether the necessary hitran files exist.
    for m in params.params_dict["mols"]:
        hfile = hitran_dir + "/" + hitran_files[m]
        if not path.isfile(hfile):
            raise ValueError("the hitran file (" + hfile +
                             ") does not exist.\n")

    #Build the executable.
    executable = exec_name[params.params_dict["lineshape"]]

    if not skip_build:

        #Determine which makefile to use based on the platform/architecture.
        tmp = params.params_dict["platform"] + "_" + \
              params.params_dict["architecture"]
        makefile = makefile_name[tmp]

        #Add make options.
        if params.params_dict["architecture"] == "cpu_openmp":
            opts = "OPENMP=on"
        else:
            opts = ""

        #Run make clean.
        run_make(build_dir,
                 makefile,
                 "clean")

        #Run make all_grtcode.
        run_make(build_dir,
                 makefile,
                 "all_grtcode " + opts)

        #Copy the executable to the run directory.
        copy_file(build_dir + "/" + executable,
                  run_dir)

    #Setup command line flags for the executable.
    if params.params_dict["platform"] == "gpu_devbox":
        args = []
    elif params.params_dict["platform"] == "gaea.c3":
        args = ["aprun","-n","1","-d","32","-cc","depth"]
    elif params.params_dict["platform"] == "gaea.c4":
        args = ["aprun","-n","1","-d","36","-cc","depth"]
    elif params.params_dict["platform"] == "theta":
        args = ["aprun","-n","1","-d","256","-j","4","-cc","depth"]

    #Add the executable.
    args += [run_dir + "/" + executable]

    #Add the auxiliary flags.
    output_file = params.params_dict["output_file"].strip()
    args += ["-a" + atmos_file,
             "-f" + params.params_dict["atmos_input_file_type"],
             "-o" + output_file,
             "-w" + str(params.params_dict["low_freq"]),
             "-W" + str(params.params_dict["high_freq"]),
             "-r" + str(params.params_dict["resolution"]),
             "-t" + str(params.params_dict["time_begin"]),
             "-T" + str(params.params_dict["time_end"]),
             "-y" + str(params.params_dict["lat_begin"]),
             "-Y" + str(params.params_dict["lat_end"]),
             "-x" + str(params.params_dict["lon_begin"]),
             "-X" + str(params.params_dict["lon_end"])]

    if params.params_dict["continuum"] == "y" or \
           params.params_dict["continuum"] == "yes":

        #Add the continuum command line flag.
        args.append("-C")

    if params.params_dict["architecture"] != "gpu":

        #Add the host only flag.
        args.append("-h")

    #Add the ppmv flags and hitran files.
    if params.params_dict["atmos_input_file_type"] == "rfmip":
        ppmv_flags = rfmip_ppmv_flag
    elif params.params_dict["atmos_input_file_type"] == "gfdl":
        ppmv_flags = gfdl_ppmv_flag

    for m in params.params_dict["mols"]:
        args.append(ppmv_flags[m])
        args.append(hitran_dir + "/" + hitran_files[m])

    #Run the executable.  Time how long the executable takes to run.  You
    #must be in the run directory to run the executable with the continuum
    #turned on.
    stdout.write("\nRun command:\n\n" + " ".join(map(str,args)) + "\n\n")
    pwd = getcwd()
    chdir(run_dir)
    start = time()
    run_executable(args)
    timing = time() - start

    #Move the output into the results directory.
    move_file(output_file,
              results_dir)
    chdir(pwd)

    return (results_dir + "/" + output_file),timing
