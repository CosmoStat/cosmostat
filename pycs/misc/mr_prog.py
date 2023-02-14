#  @file mr_filter.py
#
#  WAVELET FILTERING ROUTINE
#
#  Functions for image denoising using mr_filter c++ binary
#
#  @author Samuel Farrens
#  @version 1.0
#  @date 2015
#

import numpy as np
from os import remove
from subprocess import check_call
from subprocess import call
from datetime import datetime
from astropy.io import fits
import shlex
from pycs.misc.cosmostat_init import *
from pycs.misc.cosmostat_init import writefits


##
#  Function that calls mr_filter to perform a wavelet filtering on the
#  input data.
#
#  @param[in] data: 2D Input array.
#  @param[in] opt: List of additonal mr_transform options.
#  @param[in] path: Path for output files.
#  @param[in] remove_files: Option to remove output files.
#
#  @return Results of wavelet transform (and mr file name).
#
# %load_ext autoreload
# %autoreload 2
def mr_prog(
    data,
    prog="mr_filter",
    opt=None,
    path="./",
    remove_files=True,
    verbose=False,
    FileOut=None,
):
    # Create a unique string using the current date and time.
    # print('mr_filter ', opt)
    unique_string = datetime.now().strftime("%Y.%m.%d_%H.%M.%S")
    result = 0
    # Set the ouput file names.
    file_name = path + "mr_temp_" + unique_string
    file_fits = file_name + ".fits"
    if FileOut is not None:
        file_out = FileOut
    else:
        file_out = file_name + "_out.fits"

    # Write the input data to a fits file.
    writefits(file_fits, data)

    # print("PROG: ", prog)
    cmd = prog

    if isinstance(opt, type(None)):
        optF = " "
    else:
        optF = opt
    if verbose:
        optF = optF + " -v "

    cmd = cmd + " " + optF + " " + file_fits + " " + file_out
    if verbose:
        print("CMD = ", cmd)

    args = shlex.split(cmd)
    # print('args ', args)
    call(args)

    # Retrieve wavelet filtered data.
    result = readfits(file_out)

    # Return the mr_transform results (and the output file names).
    if remove_files:
        remove(file_fits)
        remove(file_out)
        return result
    else:
        return result
