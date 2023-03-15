#! /usr/bin/env Python
"""
Created on Sept 11 2020

@authors:  Jean-Luc Starck  
"""


import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys

sys.path.insert(0, "/Users/starck/Main/python/data/wlens/COSMOS/Data")

# Set data path
import pycs
from pycs.misc.cosmostat_init import *
from pycs.misc.mr_prog import *
from pycs.misc.im1d_tend import *
from pycs.misc.utilHSS import *
from pycs.misc.stats import *
from pycs.sparsity.mrs.mrs_tools import *


from sys import getsizeof

from pycs.sparsity.mrs.mrs_tools import *

# cmd="/Users/starck/git/cosmostat/cosmostat/src/cxx/build/wls_mcalens -F2 -i100 -N covmat.fits -v -s5 -n5  -m4 -S mica_kappa_map_cl.fits mice_g1_noisy.fits mice_g2_noisy.fits res_mice_100"


def make_healpix_map(ra, dec, weights, nside):
    pixels = hp.ang2pix(nside, theta=0.5 * np.pi - np.deg2rad(dec), phi=np.deg2rad(ra))
    bincount = np.bincount(pixels, minlength=hp.nside2npix(nside))
    bincount_weighted = np.bincount(
        pixels, minlength=hp.nside2npix(nside), weights=weights
    )
    return np.where(bincount > 0.5, bincount_weighted / bincount, hp.UNSEEN)


def get_bincount(ra, dec, nside):
    pixels = hp.ang2pix(nside, theta=0.5 * np.pi - np.deg2rad(dec), phi=np.deg2rad(ra))
    bincount = np.bincount(pixels, minlength=hp.nside2npix(nside))
    return bincount


def shape_noise_realisation(ra, dec, e1_orig, e2_orig, nside):
    gamma1_shuffle = np.copy(e1_orig)
    gamma2_shuffle = np.copy(e2_orig)
    random.shuffle(gamma1_shuffle)
    random.shuffle(gamma2_shuffle)
    e1_noise = make_healpix_map(ra, dec, gamma1_shuffle, nside)
    e2_noise = make_healpix_map(ra, dec, gamma2_shuffle, nside)
    return e1_noise, e2_noise


def shear2healpix(g1, g2, ra, dec, nside, Average):
    pixels = hp.ang2pix(nside, theta=0.5 * np.pi - np.deg2rad(dec), phi=np.deg2rad(ra))
    bincount = np.bincount(pixels, minlength=hp.nside2npix(nside))
    bincount_weighted_g1 = np.bincount(
        pixels, minlength=hp.nside2npix(nside), weights=g1
    )
    bincount_weighted_g2 = np.bincount(
        pixels, minlength=hp.nside2npix(nside), weights=g2
    )
    if Average:
        gamma1_map = np.where(bincount > 0.5, bincount_weighted_g1 / bincount, 0)
        gamma2_map = np.where(bincount > 0.5, bincount_weighted_g2 / bincount, 0)
    else:
        gamma1_map = bincount_weighted_g1
        gamma2_map = bincount_weighted_g2
    return bincount, gamma1_map, gamma2_map


# count, m_g1, m_g2 = shear2healpix(g1, g2, ra, dec, nside, Average=False)


#     a = np.random.normal(loc=0.0, scale=10.0, size=[200])


def wls_massmapping(
    mg1,
    mg2,
    covmat=None,
    SigmaGaussianNoise=None,
    NSigma=None,
    NScale=None,
    Niter=None,
    FirstDetectionScale=None,
    Signal_PowSpec=None,
    Mask=None,
    SigmaArcmin=None,
    remove_files=True,
    verbose=True,
    Prefix="xx_ks",
    opt=None,
):
    """
    Parameters
    ----------
    datag1 : TYPE
        DESCRIPTION.
    datag2 : TYPE
        DESCRIPTION.
    opt : TYPE, optional
        DESCRIPTION. The default is None.
    path : TYPE, optional
        DESCRIPTION. The default is './'.
    remove_files : TYPE, optionalAAAQ2
        DESCRIPTION. The default is True.
    verbose : TYPE, optional
        DESCRIPTION. The default is False.
    FileOut : TYPE, optional
        DESCRIPTION. The default is None.

    Returns
    -------
    None.
    """

    prog = "wls_mcalens"
    # prog_path = "/Users/starck/git/cosmostat/cosmostat/src/cxx/build/"

    # Create a unique string using the current date and time.
    # print('mr_filter ', opt)
    unique_string = "_" + datetime.now().strftime("%Y.%m.%d_%H.%M.%S")
    result = 0
    # Set the ouput file names.
    # file_name = prog_path + 'wls_temp_' + unique_string
    file_name_g1 = "wls_temp_g1" + unique_string
    file_name_g2 = "wls_temp_g2" + unique_string
    file_fits_g1 = file_name_g1 + ".fits"
    file_fits_g2 = file_name_g2 + ".fits"

    if Prefix is not None:
        file_out = Prefix
    else:
        file_out = "xx_out_ks"

    # Write the input data to a fits file.
    mrs_write(file_fits_g1, mg1)
    mrs_write(file_fits_g2, mg2)

    # print("PROG: ", prog)
    cmd = prog

    if isinstance(opt, type(None)):
        optF = " "
    else:
        optF = opt
    if verbose:
        optF = optF + " -v "
    if not isinstance(NSigma, type(None)):
        optF = optF + " -s " + str(NSigma) + " "
    if not isinstance(NScale, type(None)):
        optF = optF + " -n " + str(NScale) + " "
    if not isinstance(Niter, type(None)):
        optF = optF + " -i " + str(Niter) + " "
    if not isinstance(FirstDetectionScale, type(None)):
        optF = optF + " -F " + str(FirstDetectionScale) + " "
    if not isinstance(SigmaGaussianNoise, type(None)):
        optF = optF + " -g " + str(SigmaGaussianNoise) + " "

    if not isinstance(covmat, type(None)):
        file_name_cov = "wls_temp_cov" + unique_string
        file_fits_cov = file_name_cov + ".fits"
        mrs_write(file_fits_cov, covmat)
        optF = optF + " -N " + file_fits_cov + " "

    if not isinstance(Mask, type(None)):
        file_name_mask = "wls_temp_mask" + unique_string
        file_fits_mask = file_name_mask + ".fits"
        mrs_write(file_fits_mask, Mask)
        optF = optF + " -M " + file_fits_mask + " "

    if not isinstance(Signal_PowSpec, type(None)):
        file_name_ps = "wls_temp_sigps" + unique_string
        file_fits_ps = file_name_ps + ".fits"
        writefits(file_fits_ps, Signal_PowSpec)
        optF = optF + " -S " + file_fits_ps + " "

    cmd = cmd + " " + optF + " " + file_fits_g1 + " " + file_fits_g2 + " " + file_out
    if verbose:
        print("CMD = ", cmd)

    args = shlex.split(cmd)
    # print('args ', args)
    call(args)

    # Retrieve wavelet filtered data.
    file_out_e = file_out + "_e.fits"
    file_out_b = file_out + "_b.fits"
    result_e = mrs_read(file_out_e)
    result_b = mrs_read(file_out_b)

    if SigmaArcmin is not None:
        SigmaRad = SigmaArcmin / (60.0 * 360.0) * (np.pi * 2)
        result_e = hp.smoothing(result_e, sigma=SigmaRad, pol=False)
        result_b = hp.smoothing(result_b, sigma=SigmaRad, pol=False)

    # Remove all created intermediate files
    if remove_files:
        remove(file_fits_g1)
        remove(file_fits_g2)
        remove(file_out_e)
        remove(file_out_b)
        if not isinstance(covmat, type(None)):
            remove(file_fits_cov)
        if not isinstance(Signal_PowSpec, type(None)):
            remove(file_fits_ps)
        if not isinstance(Mask, type(None)):
            remove(file_fits_mask)

    # Return the mr_transform results (and the output file names).
    return result_e, result_b


def wls_kaiser_squires(
    mg1,
    mg2,
    SigmaArcmin=None,
    remove_files=True,
    verbose=True,
    Prefix="xx_ks",
    opt=None,
):
    optF = " -m1 "
    if not isinstance(opt, type(None)):
        optF = optF + opt
    result_e, result_b = wls_massmapping(
        mg1,
        mg2,
        SigmaArcmin=SigmaArcmin,
        remove_files=remove_files,
        verbose=verbose,
        Prefix=Prefix,
        opt=optF,
    )
    return result_e, result_b


def wls_wiener(
    mg1,
    mg2,
    covmat=None,
    SigmaGaussianNoise=None,
    Signal_PowSpec=None,
    Niter=None,
    Mask=None,
    SigmaArcmin=None,
    remove_files=True,
    verbose=True,
    Prefix="xx_wiener",
    opt=None,
):
    optF = " -m2 "
    if not isinstance(opt, type(None)):
        optF = optF + opt
    result_e, result_b = wls_massmapping(
        mg1,
        mg2,
        covmat=covmat,
        SigmaGaussianNoise=SigmaGaussianNoise,
        Signal_PowSpec=Signal_PowSpec,
        Niter=Niter,
        Mask=Mask,
        SigmaArcmin=SigmaArcmin,
        remove_files=remove_files,
        verbose=verbose,
        Prefix=Prefix,
        opt=optF,
    )
    return result_e, result_b


def wls_sparse(
    mg1,
    mg2,
    covmat=None,
    SigmaGaussianNoise=None,
    NSigma=3,
    NScale=None,
    Niter=None,
    FirstDetectionScale=None,
    Mask=None,
    SigmaArcmin=None,
    remove_files=True,
    verbose=True,
    Prefix="xx_sparse_",
    opt=None,
):
    optF = " -m3 "
    if not isinstance(opt, type(None)):
        optF = optF + opt
    result_e, result_b = wls_massmapping(
        mg1,
        mg2,
        covmat=covmat,
        SigmaGaussianNoise=SigmaGaussianNoise,
        NSigma=NSigma,
        NScale=NScale,
        Niter=Niter,
        Mask=Mask,
        FirstDetectionScale=FirstDetectionScale,
        SigmaArcmin=SigmaArcmin,
        remove_files=remove_files,
        verbose=verbose,
        Prefix=Prefix,
        opt=optF,
    )
    return result_e, result_b


class MCAlens:
    map_e = 0
    map_b = 0
    map_sparse_e = 0
    map_sparse_b = 0
    map_wiener_e = 0
    map_wiener_b = 0
    filter_wiener_e = 0
    filter_wiener_b = 0
    map_active_coef_e = 0
    map_active_coef_b = 0
    active_coef_e = 0
    active_coef_b = 0

    def __init__(self):
        self.map_e = 0


def wls_mcalens(
    mg1,
    mg2,
    covmat=None,
    SigmaGaussianNoise=None,
    Signal_PowSpec=None,
    NSigma=3,
    NScale=None,
    Niter=None,
    FirstDetectionScale=None,
    Mask=None,
    SigmaArcmin=None,
    remove_files=True,
    verbose=True,
    Prefix="xx_mcalens_",
    opt=None,
):
    optF = " -m4 "
    if not isinstance(opt, type(None)):
        optF = optF + opt

    result_e, result_b = wls_massmapping(
        mg1,
        mg2,
        covmat=covmat,
        SigmaGaussianNoise=SigmaGaussianNoise,
        Signal_PowSpec=Signal_PowSpec,
        NSigma=NSigma,
        NScale=NScale,
        Niter=Niter,
        Mask=Mask,
        FirstDetectionScale=FirstDetectionScale,
        SigmaArcmin=SigmaArcmin,
        remove_files=remove_files,
        verbose=verbose,
        Prefix=Prefix,
        opt=optF,
    )
    C = MCAlens()
    C.map_e = result_e
    C.map_b = result_b

    file_out_e = Prefix + "_sparse_e.fits"
    file_out_b = Prefix + "_sparse_b.fits"
    C.map_sparse_e = mrs_read(file_out_e)
    C.map_sparse_b = mrs_read(file_out_b)
    C.map_wiener_e = C.map_e - C.map_sparse_e
    C.map_wiener_b = C.map_b - C.map_sparse_b
    if remove_files:
        remove(file_out_e)
        remove(file_out_b)

    file_out_e = Prefix + "_filter_wiener_e.fits"
    file_out_b = Prefix + "_filter_wiener_b.fits"
    print("map_wiener_e = ", file_out_e)
    C.filter_wiener_e = readfits(file_out_e)
    print("map_wiener_b = ", file_out_b)
    C.filter_wiener_b = readfits(file_out_b)
    if remove_files:
        remove(file_out_e)
        remove(file_out_b)

    file_out_e = Prefix + "_active_coef_e.fits"
    file_out_b = Prefix + "_active_coef_b.fits"
    C.active_coef_e = readfits(file_out_e)
    C.active_coef_b = readfits(file_out_b)

    C.map_active_coef_e = np.sum(C.active_coef_e, axis=0)
    C.map_active_coef_b = np.sum(C.active_coef_b, axis=0)

    if remove_files:
        remove(file_out_e)
        remove(file_out_b)

    return C


# se,sb = wls_sparse(mg1, mg2, covmat, Nsigma=3, opt=None, remove_files=RemFiles,verbose=True, Prefix="xx_sparse_3sig")
#    ks_e = wls_massmapping(mg1,mg2, opt='-m1 -v', path='./', remove_files=RemFiles, verbose=True, FileOut=None):
#    ks_e,ks_b =  wls_kaiser_squires(mg1, mg2, Sigma=None,remove_files=RemFiles,verbose=True, Prefix="xx_ks_")
