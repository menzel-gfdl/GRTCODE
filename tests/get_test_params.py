import errno
import re
import sys

valid_params = {"architecture" : ["gpu","cpu","cpu_openmp"],
                "atmos_input_file_type" : ["gfdl","rfmip"],
                "continuum" : ["y","n","yes","no"],
                "lineshape" : ["voigt","voigt_ida","lorentz","doppler"],
                "platform" : ["gpu_devbox","gaea.c3","gaea.c4","theta"]}

valid_mols = ["h2o","co2","o3","n2o","co","ch4","o2"]

class testParams(object):
    """
    Container to hold test parameters.
    """

    def __init__(self,
                 config_file,
                 test_name):
        """
        Initialize a testParams object by parsing a config file.
        """

        #Initialize members of the class.
        self.test_name = test_name
        self.config_file = config_file
        self.params_dict = {"architecture" : "",
                            "atmos_input_file" : "",
                            "atmos_input_file_type" : "",
                            "continuum" : "",
                            "high_freq" : "",
                            "lat_begin" : "",
                            "lat_end" : "",
                            "lines" : "",
                            "lineshape" : "",
                            "lon_begin" : "",
                            "lon_end" : "",
                            "low_freq" : "",
                            "mols" : "",
                            "output_file" : "",
                            "platform" : "",
                            "resolution" : "",
                            "time_begin" : "",
                            "time_end" : ""}

        #Open the config file.
        try:
            f = open(config_file,
                     "r")
        except IOError as err:
            sys.stderr.write("Error: " + str(err.errno) + " " + err.strerror +
                                 "\n")
            raise

        #Find the test parameters in the file.
        test_found = False
        for i,line in enumerate(f):
            if not test_found:
                test_found = re.match('^\s*@test\s*=\s*' + re.escape(test_name) + '\s*$',
                                      line)
            else:
                test_end = re.match('^\s*@endtest\s*$',
                                    line)
                if test_end:
                    break
                else:
                    key_on_line = False
                    for key in self.params_dict:
                        if not self.params_dict[key]:
                            param_match = re.match('^\s*' + re.escape(key) + '\s*=\s*(.*)\s*$',
                                                   line)
                            if param_match:
                                self.params_dict[key] = param_match.group(1)
                                key_on_line = True
                                break
                    if not key_on_line:
                        all_keys_str = ""
                        for key in self.params_dict:
                            all_keys_str += "\n\t" + str(key)
                        raise ValueError("Error: line " + str(i+1) + " in" +
                                             " file " + config_file + ":\n\n" +
                                             line + "\nInside a test block,"
                                             + " lines must have the form:\n\n"
                                             + "<key> =  <value>\n\nwhere key"
                                             + " is one of:" + all_keys_str +
                                             "\n\nEach key may only appear" +
                                             " once in a test block.\n")

        #Close the file.
        f.close()

        #Make sure that the test was found in the config file.
        if not test_found:
            raise ValueError("Test " + test_name + " not found in config" +
                             " file " + config_file + ".\n")

        #Check the test configuration to make sure it is valid.
        self.check_config()

    def check_config(self):
        """
        Make sure that a valid set of test parameters were passed in.
        """

        for key in self.params_dict:

            #Make sure that the key was given a value.
            if not self.params_dict[key]:
                raise ValueError("value for " + key + " is missing from" +
                                     " test " + self.test_name + " in file " +
                                     self.config_file + "\n")

            if key in valid_params:
                if (self.params_dict[key].strip()).lower() in valid_params[key]:
                    self.params_dict[key] = (self.params_dict[key].strip()).lower()

                else:
                    all_vals_str = ""
                    for item in valid_params[key]:
                        all_vals_str += "\t" + str(item) + "\n"

                    raise ValueError("value for " + key + " (" +
                                         self.params_dict[key] + ") must" +
                                         " be one of:\n" + all_vals_str)
            else:
                if key == "high_freq" or key == "lat_end" or \
                       key == "lat_begin" or key == "lon_end" or \
                       key == "lon_begin" or key == "low_freq" or \
                       key == "time_begin" or key == "time_end":

                    try:
                        self.params_dict[key] = int(self.params_dict[key].strip())

                    except:
                        raise ValueError("the value for " + key + "(" +
                                         repr(self.params_dict[key]) +
                                         ") cannot be converted" +
                                         " to an int.\n")

                elif key == "lines":

                    try:
                        self.params_dict[key] = float(self.params_dict[key].strip())

                    except:
                        if (self.params_dict[key].strip()).lower() == "all":
                            self.params_dict[key] = (self.params_dict[key].strip()).lower()
                        else:
                            raise ValueError("the value for " + key + "(" +
                                             repr(self.params_dict[key]) +
                                             ") must be either the string" +
                                             " 'all' or a string that can be" +
                                             " converted to a float.\n")

                elif key == "mols":

                    m = self.params_dict[key].split()
                    self.params_dict[key] = []
                    for item in m:
                        tmp = (item.strip()).lower()
                        if tmp not in valid_mols:
                            all_mols_str = ""
                            for item in valid_mols:
                                all_mols_str += "\t" + str(item) + "\n"
                            raise ValueError("value for " + key + " (" +
                                             (m.strip()).lower() + ") must" +
                                             " be one of:\n" + all_mols_str)
                        else:
                            if tmp not in self.params_dict[key]:
                                self.params_dict[key].append(tmp)

                elif key == "resolution":

                    try:
                        self.params_dict[key] = float(self.params_dict[key].strip())

                    except:
                        raise ValueError("the value for " + key + "(" +
                                         repr(self.params_dict[key]) +
                                         ") cannot be converted to a float.\n")

        #Check for incompatible platform/architecture combinations.
        if self.params_dict["architecture"] == "gpu" and \
               self.params_dict["platform"] != "gpu_devbox":
            raise ValueError("architecture = gpu can only be run on" +
                             " platform = gpu_devbox.\n")

    def show(self):
        """
        Print out the members of self.
        """

        sys.stdout.write("\nTest properties:\n")
        sys.stdout.write("Test name: " + self.test_name + "\n")
        sys.stdout.write("Config file: " + self.config_file + "\n")
        for key in self.params_dict:
            sys.stdout.write(key + ": " + str(self.params_dict[key]) + "\n")
        sys.stdout.write("\n")
