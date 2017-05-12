import netCDF4 as nc
from numpy import float32, zeros_like

def verify_results(grtcode_file_path,
                   rfm_files,
                   print_out_file=True,
                   out_file_name=""):
    """
    Calculate the absolute and relative differences between GRTcode and
    RFM generated optical depth values.
    """

    #Open the netCDF GRTcode file.
    grtcode_vals = nc.Dataset(grtcode_file_path,
                              "r")

    #Get the GRTcode optical depth dimensions.
    for dim in grtcode_vals.dimensions.values():
        if dim.name == "pfull":
            grt_pfull_dim_len = dim.size
        elif dim.name == "wavenumber":
            grt_wavenumber_dim_len = dim.size

    #Make sure that the RFM optical depth dimensions match the GRTcode
    #dimensions.
    rfm_pfull_dim_len = len(rfm_files)
    if rfm_pfull_dim_len != grt_pfull_dim_len:
        raise TypeError("the number of GRTcode layers (" +
                            str(grt_pfull_dim_len)
                            + ") does not match the number of RFM layers (" +
                            str(rfm_pfull_dim_len) + ").")

    with open(rfm_files[0],"r") as f:
        for _ in range(3):
            f.readline()
        line = f.readline()
        rfm_wavenumber_dim_len = int((line.split())[0])
        if rfm_wavenumber_dim_len != grt_wavenumber_dim_len:
            raise TypeError("the number of GRTcode wavenumbers (" +
                                str(grt_wavenumber_dim_len)
                                + ") does not match the number of RFM" +
                                " wavenumbers (" +
                                str(rfm_wavenumber_dim_len) + ").")

    #Get the grtcode optical depths.
    tmp = grtcode_vals.variables["OpticalDepth"]
    grt_opt_depth = tmp[0][0][0][:][:]
    grtcode_vals.close()

    #Get the rfm optical depths.
    rfm_opt_depth = zeros_like(grt_opt_depth)
    for f in rfm_files:
        layer_num = int((f.split("."))[-1].strip())
        with open(f,"r") as fobj:

            #Skip the first 4 lines.
            for _ in range(4):
                next(fobj)

            for i,line in enumerate(fobj):
                rfm_opt_depth[layer_num-1][i-1] = float32(line.strip())

    #Calculate the absolute and relative differences between the
    #GRTcode and RFM optical depths.
    max_abs_diff = 0.0
    max_rel_diff = 0.0
    out_file_list = []

    if print_out_file:
        if not out_file_name:
            raise ValueError("Please supply a name for the output file.")

    for i in range(grt_pfull_dim_len):

        #Open the output file if necessary.
        if print_out_file:
            out_file_list.append(out_file_name + ".layer." + str(i+1))
            f = open(out_file_list[i],
                     "w")

        for j in range(grt_wavenumber_dim_len):

            #Calculate the absolute and relative differences.
            abs_diff = grt_opt_depth[i][j] - rfm_opt_depth[i][j]
            if rfm_opt_depth[i][j] != 0.0:
                rel_diff = 100.0*abs(abs_diff)/rfm_opt_depth[i][j]
            else:
                rel_diff = "N/A"

            #Update the maximum absolute and relative differences.
            if abs(abs_diff) > max_abs_diff:
                max_abs_diff = abs(abs_diff)
            if rel_diff != "N/A":
                if rel_diff > max_rel_diff:
                    max_rel_diff = rel_diff

            #Write the values to the output file.
            if print_out_file:
                f.write(str(j+1) + " " + str(grt_opt_depth[i][j]) + " " +
                            str(rfm_opt_depth[i][j]) + " " + str(abs_diff) +
                            " " + str(rel_diff) + "\n")

        #Close the output file.
        if print_out_file:
            f.close()

    return max_abs_diff, max_rel_diff, out_file_list
