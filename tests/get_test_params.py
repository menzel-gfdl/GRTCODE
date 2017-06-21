import errno
import re
import sys

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
        self.params_dict = {"atmos_input_file"      : "",
                            "atmos_input_file_type" : "",
                            "mols"                  : "",
                            "lineshape"             : "",
                            "lines"                 : "",
                            "low_freq"              : "",
                            "high_freq"             : "",
                            "resolution"            : ""}

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

        #Make sure that no keys are missing values.
        for key in self.params_dict:
            if not self.params_dict[key]:
                raise ValueError("value for " + key + " is missing from" +
                                     " test " + test_name + " in file " +
                                     config_file + "\n")

    def show(self):
        """
        Print out the members of self.
        """

        sys.stdout.write("\nTest properties:\n")
        sys.stdout.write("Test name: " + self.test_name + "\n")
        sys.stdout.write("Config file: " + self.config_file + "\n")
        for key in self.params_dict:
            sys.stdout.write(key + ": " + self.params_dict[key] + "\n")
        sys.stdout.write("\n")
