#  @file mr_gmca.py
#
#  GMCA WRAPPER ROUTINE
#
#  Functions for blind source separation using mr_gmca c++ binary
#
#  @author Jean-Luc Starck
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
#  Function that calls mr_gmca to perform blind source separation on the
#  input data.
#
#  @param[in] data: 2D Input array.
#  @param[in] opt: List of additonal mr_gmca options.
#  @param[in] path: Path for output files.
#  @param[in] remove_files: Option to remove output files.
#
#  @return Results of wavelet transform (and mr file name).
#
# %load_ext autoreload
# %autoreload 2
def mr_gmca(data, opt=None, path="./", remove_files=True, verbose=False, FileOut=None):
    # Create a unique string using the current date and time.
    # print('mr_filter ', opt)
    prog = "mr_gmca"
    unique_string = datetime.now().strftime("%Y.%m.%d_%H.%M.%S")
    result = 0
    # Set the ouput file names.
    file_name = path + "mr_temp_" + unique_string
    file_fits = file_name + ".fits"
    if FileOut is not None:
        file_out = FileOut
    else:
        file_out = file_name + "_out"

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
    file_out_source = file_out + ".fits"
    file_out_mat = path + "xx_EstMixmat.fits"
    file_out_invmat = path + "xx_InvMixingMat.fits"

    result = readfits(file_out_source)
    est_mixmat = readfits(file_out_mat)
    est_mixmat = est_mixmat.T
    est_invmixmat = readfits(file_out_invmat)
    est_invmixmat = est_invmixmat.T

    # Return the mr_transform results (and the output file names).
    if remove_files:
        remove(file_fits)
        remove(file_out_source)
        remove(file_out_mat)
        remove(file_out_invmat)

        return result, est_mixmat, est_invmixmat
    else:
        return result, est_mixmat, est_invmixmat


# Main
if __name__ == "__main__":
    # to run the test, need to install scikit-image and import the two following pakages
    from skimage import data, color
    from skimage.transform import resize
    import matplotlib.pyplot as plt
    from pycs.sparsity.sparse2d.bss_eval import *

    sources = load_source_images()
    # mixed_images, A = mix_sources_images(sources)
    #  np.random.seed(0)
    mixed_images, A = mix_sources_images_noise(sources, noise_level=0.1)
    info(sources)

    print("Mixing Matrix (A):\n", A)
    print("Sources shape:", sources.shape)  # (2, H, W)
    print("Mixed images shape:", mixed_images.shape)  # (3, H, W)

    optF = "-S2 -K3 -t14 -n5"  # with bi-orthogonal WT abd final denoising at 3sigma
    optStarlet = "-S2 -K3 -t2 -n5"  # with starlet  and final denoising at 3sigma
    optCurvelet = (
        "-E1 -S2 -K3  -t28 -n5 "  # with curvelet abd final denoising at 3sigma
    )

    corrPerm = False
    verbose = False
    SRec, Emat, Eimat = mr_gmca(
        mixed_images, opt=optF, remove_files=False, verbose=verbose
    )
    SRec = reorder_and_fix_sign(sources, SRec)
    # print(" ==> Bi-Orth Wavelet Source err = ", compute_sdr(sources, SRec))
    # print("A_true shape:", A.shape)
    # print("A_est  shape:", Emat.shape)
    # error = amari_error(A, Emat)
    # print(" ==> Mixing matrix err = ", amari_error(A, Emat))
    # print("Mixing Matrix (A):\n", Emat)
    CA, NMSE = evaluate(A, sources, Emat, SRec, corrPerm=corrPerm)
    print(" ==> Bi-Orth Wavelet Source err: CA = %.4f | NMSE= %.4f" % (CA, NMSE))

    # Visualization
    fig, axs = plt.subplots(1, 5, figsize=(15, 5))
    axs[0].imshow(sources[0], cmap="gray")
    axs[0].set_title("Source Image 1")
    axs[1].imshow(sources[1], cmap="gray")
    axs[1].set_title("Source Image 2")
    for i in range(3):
        axs[i + 2].imshow(mixed_images[i], cmap="gray")
        axs[i + 2].set_title(f"Mixed Image {i+1}")
    fig, axs = plt.subplots(1, 2, figsize=(15, 5))
    axs[0].imshow(SRec[0], cmap="gray")
    axs[0].set_title("7/9 WT GMCA Image 1")
    axs[1].imshow(SRec[1], cmap="gray")
    axs[1].set_title("7/9 WT GMCA Image 2")
    for ax in axs:
        ax.axis("off")
    plt.tight_layout()
    plt.show()

    # ----- STARLET ------
    SRec, Emat, Eimat = mr_gmca(
        mixed_images, opt=optStarlet, remove_files=False, verbose=verbose
    )
    SRec = reorder_and_fix_sign(sources, SRec)
    CA, NMSE = evaluate(A, sources, Emat, SRec, corrPerm=corrPerm)
    print("==> Starlet Source err: CA = %.4f | NMSE =  %.4f" % (CA, NMSE))

    # Visualization
    fig, axs = plt.subplots(1, 2, figsize=(15, 5))
    axs[0].imshow(SRec[0], cmap="gray")
    axs[0].set_title("Starlet GMCA Image 1")
    axs[1].imshow(SRec[1], cmap="gray")
    axs[1].set_title("Starlet GMCA Image 2")
    for ax in axs:
        ax.axis("off")
    plt.tight_layout()
    plt.show()

    # ----- CURVELET ------
    SRec, Emat, Eimat = mr_gmca(
        mixed_images, opt=optCurvelet, remove_files=False, verbose=verbose
    )
    SRec = reorder_and_fix_sign(sources, SRec)
    CA, NMSE = evaluate(A, sources, Emat, SRec, corrPerm=corrPerm)
    print("==> Curvelet Source err: CA = %.4f | NMSE =  %.4f" % (CA, NMSE))

    # Visualization
    fig, axs = plt.subplots(1, 2, figsize=(15, 5))
    axs[0].imshow(SRec[0], cmap="gray")
    axs[0].set_title("Curvelet GMCA Image 1")
    axs[1].imshow(SRec[1], cmap="gray")
    axs[1].set_title("Curvelet GMCA Image 2")
    for ax in axs:
        ax.axis("off")
    plt.tight_layout()
    plt.show()
