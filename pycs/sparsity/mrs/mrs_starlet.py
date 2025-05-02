import numpy as np
import random

import os, sys
from scipy import ndimage
import healpy as hp
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.io import fits
from importlib import reload
from pycs.misc.cosmostat_init import *
from pycs.misc.mr_prog import *
from pycs.sparsity.mrs.mrs_tools import *

def mrs_starlet(map, nscale=None, lmax=None):
    nside = gnside(map)
    if nscale is None:
        Ns = np.int64(np.log2(nside) - 2)
    else:
        Ns = nscale

    npix = map.shape[0]
    w = wt_trans(map, lmax=lmax,nscales=Ns-1)
    trans = w.T
    return trans

def mrs_istarlet(trans):
    r = np.sum(trans, axis=0)
    return r


def mrs_uwttrans(map, nscale=None, lmax=None, opt=None, verbose=False, path="./", progpath=None, cxx=False):
    nside = gnside(map)
    if nscale is None:
        Ns = np.log2(nside) - 2
    else:
        Ns = nscale

    if cxx:
        optParam = " "
        if opt is not None:
            optParam = " " + opt
        if lmax is not None:
            optParam = " -l " + str(lmax) + optParam
        if nscale is not None:
            optParam = " -n " + str(nscale) + optParam
        if progpath is None:
            prog = "mrs_uwttrans"
        else:
            prog = progpath + "mrs_uwttrans"
        p = mrs_prog(
            map,
            prog=prog,
            verbose=verbose,
            opt=optParam,
            OutputFormatisHealpix=False,
            path=path,
        )
    else:
        npix = map.shape[0]
        w = wt_trans(map, lmax=lmax,nscales=Ns-1)
        p = np.zeros(Ns, npix)
        for j in range(Ns):
            print(j+1)
            p[j,:] = w[:,j]

    return p


def mrs_uwtrecons(Tmap, lmax=None, opt=None, verbose=False, path="./", progpath=None):
    optParam = " "
    if opt is not None:
        optParam = " " + opt
    if lmax is not None:
        optParam = " -l " + str(lmax) + optParam
    if progpath is None:
        prog = "mrs_uwttrans"
    else:
        prog = progpath + "mrs_uwttrans -r "
    p = mrs_prog(
        Tmap,
        prog=prog,
        verbose=verbose,
        opt=optParam,
        InputFormatisHealpix=False,
        OutputFormatisHealpix=True,
        path=path,
    )
    return p



# Wavelet filtering

def spline2(size, l, lc):
    """
    Compute a non-negative decreasing spline, with value 1 at index 0.

    Parameters
    ----------
    size: int
        size of the spline
    l: float
        spline parameter
    lc: float
        spline parameter

    Returns
    -------
    np.ndarray
        (size,) float array, spline
    """

    res = np.arange(0, size+1)
    res = 2*l*res/(lc*size)
    res = (3/2) * 1/12 * (abs(res-2)**3 - 4*abs(res-1)**3 + 6*abs(res)**3 - 4*abs(res+1)**3 + abs(res+2)**3)
    return res


def compute_h(size, lc):
    """
    Compute a low-pass filter.

    Parameters
    ----------
    size: int
        size of the filter
    lc: float
        cutoff parameter

    Returns
    -------
    np.ndarray
        (size,) float array, filter
    """

    tab1 = spline2(size, 2*lc, 1)
    tab2 = spline2(size, lc, 1)
    h = tab1/(tab2+1e-6)
    h[np.int64(size/(2*lc)):size] = 0
    return h


def compute_g(size, lc):
    """
    Compute a high-pass filter.

    Parameters
    ----------
    size: int
        size of the filter
    lc: float
        cutoff parameter

    Returns
    -------
    np.ndarray
        (size,) float array, filter
    """

    tab1 = spline2(size, 2*lc, 1)
    tab2 = spline2(size, lc, 1)
    g = (tab2-tab1)/(tab2+1e-6)
    g[np.int64(size/(2*lc)):size] = 1
    return g


def get_wt_filters(lmax, nscales):
    """Compute wavelet filters.

    Parameters
    ----------
    lmax: int
        maximum l
    nscales: int
        number of wavelet detail scales

    Returns
    -------
    np.ndarray
        (lmax+1,nscales+1) float array, filters
    """

    wt_filters = np.ones((lmax+1, nscales+1))
    wt_filters[:, 1:] = np.array([compute_h(lmax, 2**scale) for scale in range(nscales)]).T
    wt_filters[:, :nscales] -= wt_filters[:, 1:(nscales+1)]
    return wt_filters


def wt_trans(inputs, nscales=3, lmax=None, alm_in=False, nside=None, alm_out=False):
    """Wavelet transform an array.

    Parameters
    ----------
    inputs: np.ndarray
        (p,) or (n,p) float array, map or stack of n maps / if alm_in, (t,) or (n,t) complex array, alm or stack
        of n alms
    nscales: int
        number of wavelet detail scales
    lmax: int
        maximum l (default: 3*nside / if alm_in, deduced from inputs)
    alm_in: bool
        inputs is alm
    nside: int
        nside of the output Healpix maps (default: deduced from maps)
    alm_out: bool
        output is alm

    Returns
    -------
    np.ndarray
        (p,nscales+1) or (n,p,scales+1) float array, wavelet transform of the input array or stack of the wavelet
        transforms of the n input arrays / if alm_out, (t,nscales+1) or (n,t,scales+1) complex array, alm of the
        wavelet transform of the input array or stack of the alms of the wavelet transforms of the n input arrays
    """
    dim_inputs = len(np.shape(inputs))
    maps = None  # to remove warnings

    if alm_in:
        alms = inputs
        if nside is None and not alm_out:
            raise ValueError("nside is missing")
        if not alm_out:
            maps = alm2map(alms, nside)
        if lmax is None:
            lmax = hp.Alm.getlmax(np.shape(alms)[-1])

    else:
        maps = inputs
        if dim_inputs == 1:
            nside = hp.get_nside(maps)
        else:
            nside = hp.get_nside(maps[0, :])
        if lmax is None:
            lmax = 3 * nside
        alms = map2alm(maps, lmax=lmax)

    if not alm_out:
        l_scale = maps.copy()
        if dim_inputs == 1:
            npix = len(maps)
            wts = np.zeros((npix, nscales + 1))
        else:
            npix = np.shape(maps)[1]
            wts = np.zeros((np.shape(maps)[0], npix, nscales + 1))
    else:
        l_scale = alms.copy()
        if dim_inputs == 1:
            npix = np.size(alms)
            wts = np.zeros((npix, nscales + 1), dtype='complex')
        else:
            npix = np.shape(alms)[1]
            wts = np.zeros((np.shape(maps)[0], npix, nscales + 1), dtype='complex')

    scale = 1
    for j in range(nscales):
        h = compute_h(lmax, scale)
        if not alm_out:
            m = alm2map(alm_product(alms, h), nside)
        else:
            m = alm_product(alms, h)
        h_scale = l_scale - m
        l_scale = m
        if dim_inputs == 1:
            wts[:, j] = h_scale
        else:
            wts[:, :, j] = h_scale
        scale *= 2

    if dim_inputs == 1:
        wts[:, nscales] = l_scale
    else:
        wts[:, :, nscales] = l_scale

    return wts


def wt_rec(wts):
    """Reconstruct a wavelet decomposition.

    Parameters
    ----------
    wts: np.ndarray
        (p,nscales+1) or (n,p,scales+1) float array, wavelet transform of a map or stack of the wavelet transforms of n
        maps

    Returns
    -------
    np.ndarray
        (p,) or (n,p,) float array, reconstructed map or stack of n reconstructed maps
    """

    return np.sum(wts, axis=-1)


# Plots

def mrs_tv(maps, log=False, unit='', title='', minimum=None, maximum=None, cbar=True):
    """Plot one or more Healpix maps in Mollweide projection.

    Parameters
    ----------
    maps: np.ndarray
        (p,) or (n,p) float array, map or stack of n maps
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

    if minimum is None:
        minimum = np.min(maps)
        
    if maximum is None:
        maximum = np.max(maps)
        
    if not log:
        def f(x): return x
    else:
        def f(x): return np.log10(x - minimum + 1)
    for i in range(np.shape(maps)[0]):
        if title:
            tit = title + ": Scale " + str(i+1)
        else:
            tit = "Scale " + str(i+1)
        hp.mollview(f(maps[i, :]), fig=None, unit=unit, title=tit, min=f(minimum), max=f(maximum), cbar=cbar)


