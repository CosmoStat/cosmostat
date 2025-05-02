import numpy as np
import random

import os, sys
from scipy import ndimage
import healpy as hp
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from astropy.io import fits
from importlib import reload
from pycs.misc.cosmostat_init import *
from pycs.misc.mr_prog import *


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


def mrs_read(FN):
    return hp.read_map(FN)


def mrs_write(FN, mapin):
    hp.write_map(FN, mapin, overwrite=True)


def rims(FN):
    return hp.read_map(FN)


def mrs_resize(mapin, nsideout):
    k = hp.ud_grade(mapin, nsideout)
    return k


# smoothing with sigma in arcmin
def smooth(map, sigma):
    s = hp.smoothing(mapin, sigma=sigma / (360.0 * 60.0) * (np.pi * 2), pol=False)


#        lut='rainbow'  #  'inferno'   'gist_stern'
def tvs(
    mapin, min=None, max=None, title=None, sigma=None, lut=None, filename=None, dpi=100
):
    if sigma is None:
        hp.mollview(mapin, max=max, min=min, title=title, cmap=lut)
    else:
        s = hp.smoothing(mapin, sigma=sigma / (360.0 * 60.0) * (np.pi * 2), pol=False)
        hp.mollview(s, max=max, min=min, title=title, cmap=lut)
    if not isinstance(filename, type(None)):
        print("Filename = ", filename)
        plt.savefig(filename, dpi=dpi)


def get_nside(Npix):
    return hp.npix2nside(Npix)


def gnside(data):
    npix = data.shape[0]
    nside = hp.npix2nside(npix)
    return nside


def pixel_size(nside):
    # Return the pixel size of a healpix map in arc minutes
    # SKI_SURFACE IN SQUARE DEGREES =  4. * !PI * (360. / (2*!PI))^2 = 41253
    psize = 41253.0 / (float(nside) ** 2.0 * 12.0) * 60.0**2.0
    return np.sqrt(psize)


def l2amin(l):
    a = 1.0 / l
    a = a * 180.0 * 60.0 / np.pi
    return a


def amin2l(a):
    ar = a / (180.0 * 60.0) * np.pi
    l = 1.0 / ar
    return l


def g2eb(g1, g2):
    nside = gnside(g1)
    (ae, ab) = hp.map2alm_spin((g1, g2), 2)
    ke = hp.alm2map(ae, nside, pol=False)
    kb = hp.alm2map(ab, nside, pol=False)
    return ke, kb


def g2k(g1, g2):
    nside = gnside(g1)
    (ae, ab) = hp.map2alm_spin((g1, g2), 2)
    ke = hp.alm2map(ae, nside, pol=False)
    return ke

def k2g(ke):
    nside = gnside(ke)
    ae = hp.map2alm(ke, 1, pol=False)
    ab = np.copy(ae) * 0.0
    (g1, g2) = hp.alm2map_spin((ae, ab), 2, lmax=lmax)
    return g1, g2


# it seems that hp.alm2map_spin crashes.
def eb2g(ke, kb):
    nside = gnside(ke)
    lmax = nside * 3 - 1
    ae = hp.map2alm(ke, 1, pol=False)
    ab = hp.map2alm(kb, 1, pol=False)
    (g1, g2) = hp.alm2map_spin((ae, ab), nside, 2, lmax)
    return g1, g2


def mrs_prog(
    data,
    prog="mrs_powspec",
    opt=None,
    path="./",
    remove_files=True,
    verbose=False,
    FileOut=None,
    InputFormatisHealpix=True,
    OutputFormatisHealpix=True,
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
    if InputFormatisHealpix:
        mrs_write(file_fits, data)
    else:
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
    if OutputFormatisHealpix:
        result = mrs_read(file_out)
    else:
        result = readfits(file_out)

    # Return the mr_transform results (and the output file names).
    if remove_files:
        remove(file_fits)
        remove(file_out)
        return result
    else:
        return result


def mrs_powspec(map, verbose=False):
    p = mrs_prog(map, prog="mrs_powspec", verbose=verbose, OutputFormatisHealpix=False)
    return p


def mrs_smooth(map, opt=None, verbose=False):
    p = mrs_prog(
        map, prog="mrs_smooth", verbose=verbose, opt=opt, OutputFormatisHealpix=True
    )
    return p


def mrs_almtrans(map, lmax=None, opt=None, verbose=False):
    optParam = " -T "
    if opt is not None:
        optParam = " -T " + opt
    if lmax is not None:
        optParam = " -l " + str(lmax) + optParam
    p = mrs_prog(
        map,
        prog="mrs_almtrans",
        verbose=verbose,
        opt=optParam,
        OutputFormatisHealpix=False,
    )
    return p


def mrs_almrec(map, opt=None, verbose=False, nside=None):
    optParam = " -T "
    if opt is not None:
        optParam = " -T " + opt
    if nside is not None:
        optParam = " -n " + str(nside) + optParam
    p = mrs_prog(
        map,
        prog="mrs_almrec",
        verbose=verbose,
        opt=optParam,
        InputFormatisHealpix=False,
        OutputFormatisHealpix=True,
    )
    return p


def tol(map, lmax_amin, amin=False):
    ns = gnside(map)
    lmax = lmax_amin
    if amin is True:
        lmax = amin2l(lmax_amin)
    a = mrs_almtrans(map, lmax=lmax)
    b = mrs_almrec(a, nside=ns)
    return b


# Spherical harmonic transform : Code from Remi Carloni Gertosio

def map2alm(maps, lmax=None, iter=3):
    """Computes the alm of a Healpix map.

    Parameters
    ----------
    maps: np.ndarray
        (p,) or (n,p) float array, map or stack of n maps in Healpix representation
    lmax: int
        maximum l of the alm (default: 3*nside)
    iter: 3
        number of iterations

    Returns
    -------
    np.ndarray
        (t,) or (n,t) complex array, alm or stack of n alms
    """

    if len(np.shape(maps)) == 1:
        if lmax is None:
            lmax = 3*hp.get_nside(maps)
        return hp.sphtfunc.map2alm(maps, lmax=lmax, iter=iter)

    n = np.shape(maps)[0]
    if lmax is None:
        lmax = 3*hp.get_nside(maps[0, :])
    return np.array([hp.sphtfunc.map2alm(maps[i, :], lmax=lmax, iter=iter) for i in range(n)])


def alm2map(alms, nside):
    """Computes a Healpix map given the alm.

    Parameters
    ----------
    alms: np.ndarray
        (t,) or (n,t) complex array, alm or stack of n alms
    nside: int
        nside of the output Healpix maps

    Returns
    -------
    np.ndarray
        (p,) or (n,p) float array, map or stack of n maps in Healpix representation
    """

    if len(np.shape(alms)) == 1:
        return hp.alm2map(alms, nside)

    n = np.shape(alms)[0]
    return np.array([hp.sphtfunc.alm2map(alms[i, :], nside) for i in range(n)])


def alm_product(alms, filters):
    """Apply an isotropic filter on an alm.

    Parameters
    ----------
    alms: np.ndarray
        (t,) or (n,t) complex array, alm or stack of n alms
    filters: np.ndarray
        (lmax+1,) or (n,lmax+1) float array, isotropic filter or stack of n isotropic filters (one filter per source) in
         spherical harmonic domain

    Returns
    -------
    np.ndarray
        (t,) or (n,t) complex array, filtered alm or stack of n filtered alms
    """

    dim_filters = len(np.shape(filters))
    dim_alms = len(np.shape(alms))

    if dim_filters == 1 and dim_alms == 1:
        return hp.sphtfunc.smoothalm(alms, beam_window=filters, inplace=False)

    n = np.shape(alms)[0]

    if dim_filters == 1:
        return np.array([hp.sphtfunc.smoothalm(alms[i, :], beam_window=filters, verbose=False, inplace=False)
                         for i in range(n)])

    return np.array([hp.sphtfunc.smoothalm(alms[i, :], beam_window=filters[i, :], verbose=False, inplace=False)
                     for i in range(n)])


def convolve(maps, filters, lmax=None, nside=None):
    """Convolve maps with filters.

    Parameters
    ----------
    maps: np.ndarray
        (p,) or (n,p) float array, map or stack of n maps in Healpix representation
    filters: np.ndarray
        (lmax+1,) or (n,lmax+1) float array, isotropic filter or stack of n isotropic filters (one filter per source)
    lmax: int
        maximum l of the filtering (default: deduced from filters)
    nside: int
        nside of the output Healpix maps (default: deduced from maps)

    Returns
    -------
    maps: np.ndarray
        (p,) or (n,p) float array, convolved map or stack of n convolved maps in Healpix representation
    """

    if lmax is not None:
        if len(np.shape(filters)) == 1:
            lmax = len(filters) - 1
        else:
            lmax = np.shape(filters)[1] - 1

    alms = map2alm(maps, lmax=lmax)

    alms = alm_product(alms, filters)

    if nside is None:
        nside = hp.get_nside(maps)

    return alm2map(alms, nside=nside)


def anafast(maps, lmax=None, iter=3):
    """Computes the angular power spectrum of a Healpix map.

    Parameters
    ----------
    maps: np.ndarray
        (p,) or (n,p) float array, map or stack of n maps in Healpix representation
    lmax: int
        maximum l of the angular power spectrum (default: 3*nside of maps)
    iter: 3
        number of iterations

    Returns
    -------
    np.ndarray
        (lmax+1,) or (n,lmax+1) float array, angular power spectrum or stack of n angular power spectra
    """

    if len(np.shape(maps)) == 1:
        if lmax is None:
            lmax = 3*hp.get_nside(maps)
        return hp.sphtfunc.anafast(maps, lmax=lmax, iter=iter)

    n = np.shape(maps)[0]
    if lmax is None:
        lmax = 3 * hp.get_nside(maps[0, :])
    return np.array([hp.sphtfunc.anafast(maps[i, :], lmax=lmax) for i in range(n)])


def alm2cl(alms):
    """Computes the angular power spectrum from an alm.

    Parameters
    ----------
    alms: np.ndarray
        (t,) or (n,t) complex array, alm or stack of n alms

    Returns
    -------
    np.ndarray
        (lmax+1,) or (n,lmax+1) float array, angular power spectrum or stack of n angular power spectra
    """

    if len(np.shape(alms)) == 1:
        return hp.sphtfunc.alm2cl(alms)

    n = np.shape(alms)[0]
    return np.array([hp.sphtfunc.alm2cl(alms[i, :]) for i in range(n)])


# Alm index computation

def getsize(lmax):
    """Returns the size of the array needed to store alm up to lmax.

    Parameters
    ----------
    lmax: int
        maximum l of the alm

    Returns
    -------
    int
        size of the array needed to store alm up to lmax

    """

    return hp.Alm.getsize(lmax)


def getlm(lmax):
    """Get the mapping of an alm.

    Parameters
    ----------
    lmax: int
        maximum l of the alm

    Returns
    -------
    (np.ndarray,np.ndarray)
        l to index map,
        m to index map
    """

    return hp.Alm.getlm(lmax)


def npix2nside(npix):
    """
    Give the nside parameter for the given number of pixels.

    Parameters
    ----------
    npix: int
        number of pixels

    Returns
    -------
    nside: int
        nside
    """

    return hp.npix2nside(npix)

 
# Plots

def mollview(maps, maps2=None, log=False, unit='', title='', minimum=None, maximum=None, cbar=True):
    """Plot one or more Healpix maps in Mollweide projection.

    Parameters
    ----------
    maps: np.ndarray
        (p,) or (n,p) float array, map or stack of n maps
    maps2: np.ndarray
        (p,) or (n,p) float array, second map or stack of n maps, optional
    log: bool
        logarithmic scale
    unit: str
        unit of the data
    title: str
        title of the plots
    minimum: float
        minimum range value (default: min(maps, maps2))
    maximum: float
        maximum range value (default: max(maps, maps2))
    cbar: bool
        show color bar

    Returns
    -------
    None
    """

    if len(np.shape(maps)) == 1:
        maps = np.expand_dims(maps, axis=0)
        if maps2 is not None:
            maps2 = np.expand_dims(maps2, axis=0)
    if minimum is None:
        minimum = np.min(maps)
        if maps2 is not None:
            minimum = np.min([minimum, np.min(maps2)])
    if maximum is None:
        maximum = np.max(maps)
        if maps2 is not None:
            maximum = np.max([maximum, np.max(maps2)])
    if not log:
        def f(x): return x
    else:
        def f(x): return np.log10(x - minimum + 1)
    for i in range(np.shape(maps)[0]):
        hp.mollview(f(maps[i, :]), fig=None, unit=unit, title=title, min=f(minimum), max=f(maximum), cbar=cbar)
        if maps2 is not None:
            hp.mollview(f(maps2[i, :]), fig=None, unit=unit, title=title, min=f(minimum), max=f(maximum), cbar=cbar)

def view_spec(inputs, lmax=None, alm_in=False):
    """Plot the angular power spectrum of one or several maps.

    Parameters
    ----------
    inputs: np.ndarray
        (p,) or (n,p) float array, map or stack of n maps / if alm_in, (t,) or (n,t) complex array, alm or stack
        of n alms
    lmax: int
        maximum l (default: 3*nside / if alm_in, deduced from inputs)
    alm_in: bool
        inputs is alm

    Returns
    -------
    None
    """

    if len(np.shape(inputs)) == 1:
        inputs = np.expand_dims(inputs, axis=0)

    if not alm_in:
        cls = anafast(inputs, lmax=lmax)
    else:
        cls = alm2cl(inputs)

    plt.figure()
    for i in range(np.shape(inputs)[0]):
        plt.semilogy(cls[i, :], label='Source '+str(i+1))
    plt.xlabel('$l$')
    plt.ylabel('$c_l$')
    if np.shape(inputs)[0] != 1:
        plt.legend()

# Miscellaneous

def getidealbeam(lmax, cutmin=None, cutmax=None):
    """Compute a beam, with value 1 until a first cutoff frequency and 0 after a second cutoff frequency. The transition
    is computed with a spline.

    Parameters
    ----------
    lmax: int
        maximum l
    cutmin: int
        frequency below which filter is 1 (default: int((lmax+1)/4))
    cutmax: int
        frequency above which filter is 0 (default: int((lmax+1)/2))

    Returns
    -------
    np.ndarray
        (lmax+1,) float array, filter
    """

    if cutmin is None:
        cutmin = np.int64((lmax+1)/4)
    if cutmax is None:
        cutmax = np.int64((lmax+1)/2)
    bl = np.zeros(lmax+1)
    bl[0:cutmin] = 1
    bl[cutmin:cutmax] = spline2(cutmax-cutmin-1, 1, 1)
    return bl


def getbeam(fwhm=100, lmax=512):
    """Get a spherical Gaussian-shaped beam.

    Parameters
    ----------
    fwhm: float
        full width at half maximum in the harmonic space (in terms of l)
    lmax: int
        maximum l

    Returns
    -------
    np.ndarray
        (lmax+1,) float array, Gaussian-shaped beam
    """
    
    tor = 0.0174533
    if len(np.shape(fwhm)) == 1:
        fwhm = np.expand_dims(fwhm, axis=1)
    F = fwhm / 60 * tor
    l = np.arange(0, lmax+1)
    ell = l * (l + 1)
    bl = np.exp(-ell * F * F / 16 / np.log(2))
    return bl
