
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
from pycs.tools.cosmostat_init import *

def make_healpix_map(ra, dec, weights, nside):
    pixels= hp.ang2pix(nside,theta = 0.5*np.pi - np.deg2rad(dec), phi = np.deg2rad(ra))
    bincount = np.bincount(pixels, minlength = hp.nside2npix(nside))
    bincount_weighted = np.bincount(pixels, minlength = hp.nside2npix(nside), weights=weights)
    return np.where(bincount>0.5, bincount_weighted/bincount, hp.UNSEEN)

def get_bincount(ra, dec, nside):
    pixels= hp.ang2pix(nside,theta = 0.5*np.pi - np.deg2rad(dec), phi = np.deg2rad(ra))
    bincount = np.bincount(pixels, minlength = hp.nside2npix(nside))
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
    s= hp.smoothing(mapin, sigma=sigma/(360.*60.) * (np.pi*2),pol=False)

def tvs(mapin,min=None,max=None,title=None,sigma=None):
    if sigma is None:
        hp.mollview(mapin,max=max,min=min, title=title)
    else:
       s= hp.smoothing(mapin, sigma=sigma/(360.*60.) * (np.pi*2),pol=False)
       hp.mollview(s,max=max,min=min, title=title)

def get_nside(Npix):
    return hp.npix2nside(Npix)

def gnside(data):
    npix = data.shape[0]
    nside = hp.npix2nside(npix)
    return nside

def pixel_size(nside):
    # Return the pixel size of a healpix map in arc minutes
    # SKI_SURFACE IN SQUARE DEGREES =  4. * !PI * (360. / (2*!PI))^2 = 41253
    psize = 41253. / (float(nside)**2.*12.) * 60.**2.
    return np.sqrt(psize)

def l2amin(l):
    a = 1. / l
    a  =  a * 180.* 60. / np.pi
    return a

def amin2l(a):
    ar =  a / (180.* 60.) * np.pi
    l = 1. / ar
    return l

def g2eb(g1,g2):
    nside = gnside(g1)
    (ae,ab) = hp.map2alm_spin((g1,g2), 2)
    ke= hp.alm2map(ae, nside, pol=False)
    kb= hp.alm2map(ab, nside, pol=False)
    return ke,kb

def g2k(g1,g2):
    nside = gnside(g1)
    (ae,ab) = hp.map2alm_spin((g1,g2), 2)
    ke= hp.alm2map(ae, nside, pol=False)
    return ke

def k2g(ke):
    nside = gnside(ke)
    ae = hp.map2alm(ke, 1,pol=False)
    ab = np.copy(ae) * 0.
    (g1,g2) = hp.alm2map_spin((ae,ab), 2, lmax=lmax)
    return g1,g2

# it seems that hp.alm2map_spin crashes.
def eb2g(ke,kb):
    nside = gnside(ke)
    lmax=nside*3 - 1
    ae = hp.map2alm(ke, 1, pol=False)
    ab = hp.map2alm(kb, 1, pol=False)
    (g1,g2) = hp.alm2map_spin( (ae,ab), nside, 2, lmax)
    return g1,g2


