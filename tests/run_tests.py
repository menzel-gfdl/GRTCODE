#!/usr/bin/env python

from optparse import OptionParser
from sys import stdout
from get_test_params import testParams
from run_grtcode import run_grtcode
from run_rfm import run_rfm
#from verify_results import verify_results

if __name__ == "__main__":

    #Parser command line arguments.
    parser = OptionParser()
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
    parser.add_option("-s",
                      "--skipbuild",
                      dest="skipbuild",
                      action="store_true",
                      default=False)
    options,args = parser.parse_args()

    #Check inputs.
    if options.configFile == None:
        stdout.write("\nUsing default config file (tests.config).\n")
        options.configFile = "tests.config"
    if options.testName == None:
        stdout.write("\nPreforming @test = standard.\n")
        options.testName = "standard"

    #Parse the config file to get the test parameters.
    testObject = testParams(options.configFile,
                            options.testName)
    testObject.show()

    #Run the test using GRTcode.
    grt_output_file, grt_timing = run_grtcode(testObject,
                                              "..",
                                              skip_build=options.skipbuild)

    #Write out timing results.
    num_columns = (testObject.params_dict["lat_end"] - \
                  testObject.params_dict["lat_begin"] + 1)* \
                  (testObject.params_dict["lon_end"] - \
                  testObject.params_dict["lon_begin"] + 1)* \
                  (testObject.params_dict["time_end"] - \
                  testObject.params_dict["time_begin"] + 1)
    with open("grtcode.timings","a") as f:
        f.write(testObject.params_dict["platform"] + "," +
                testObject.params_dict["architecture"] + "," +
                str(num_columns) + "," +
                str(testObject.params_dict["resolution"]) + "," +
                str(grt_timing) + "\n")

    stdout.write("\nOutput file located at: " + grt_output_file + "\n")
    stdout.write("\nGRTcode runtime (s): " + str(grt_timing) + "\n")
    exit()

    #Run the test using RFM.
    rfm_timing,rfm_output_files = run_rfm(testObject.mols,
                                          "layer_cond",
                                          "../",
                                          testObject.lineShape,
                                          testObject.lowFreq,
                                          testObject.highFreq,
                                          testObject.res,
                                          lines=tLines,
                                          forceBuild=False)

    #Verify the results.
    max_abs_diff, max_rel_diff, out_files = verify_results(grt_output_file,
                                                           rfm_output_files,
                                                           True,
                                                           "optical_depths")

    #Write out differences and timings to stdout.
    stdout.write("\nMax absolute difference: " + str(max_abs_diff) +
                     "\nMax relative difference: " + str(max_rel_diff) +
                     "\n"
                     "\nTimings: \nGRTcode runtime (s): " + str(grt_timing) +
                     "\nRFM runtime (s):     " +  str(rfm_timing) + "\n")
