
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

#        lut='rainbow'  #  'inferno'   'gist_stern'
def tvs(mapin,min=None,max=None,title=None,sigma=None,lut=None,filename=None,dpi=100):
    if sigma is None:
        hp.mollview(mapin,max=max,min=min, title=title,cmap=lut)
    else:
       s= hp.smoothing(mapin, sigma=sigma/(360.*60.) * (np.pi*2),pol=False)
       hp.mollview(s,max=max,min=min, title=title,cmap=lut)
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


def mrs_prog(data, prog="mrs_powspec", opt=None, path='./', remove_files=True, verbose=False, FileOut=None, InputFormatisHealpix=True, OutputFormatisHealpix=True):

    # Create a unique string using the current date and time.
    # print('mr_filter ', opt)
    unique_string = datetime.now().strftime('%Y.%m.%d_%H.%M.%S')
    result=0
    # Set the ouput file names.
    file_name = path + 'mr_temp_' + unique_string
    file_fits = file_name + '.fits'
    if FileOut is not None:
        file_out = FileOut
    else:
        file_out = file_name + '_out.fits'

    # Write the input data to a fits file.
    if InputFormatisHealpix:
        mrs_write(file_fits, data)
    else:
        writefits(file_fits, data)


    # print("PROG: ", prog)
    cmd = prog

    if isinstance(opt, type(None)):
        optF=' '
    else:
        optF= opt
    if verbose:
        optF = optF + " -v "

    cmd = cmd + " " + optF + " "  + file_fits + " "   + file_out
    if verbose:
        print ('CMD = ', cmd)

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
    p = mrs_prog(map, prog="mrs_smooth", verbose=verbose, opt=opt, OutputFormatisHealpix=True)
    return p

def mrs_almtrans(map, lmax=None, opt=None, verbose=False):
    optParam = ' -T '
    if opt is not None:
        optParam = ' -T ' + opt
    if lmax is not None:
        optParam = ' -l ' + str(lmax) + optParam
    p = mrs_prog(map, prog="mrs_almtrans", verbose=verbose, opt=optParam, OutputFormatisHealpix=False)
    return p

def mrs_almrec(map, opt=None, verbose=False,nside=None):
    optParam = ' -T '
    if opt is not None:
        optParam = ' -T ' + opt
    if nside is not None:
        optParam = ' -n ' + str(nside) + optParam
    p = mrs_prog(map, prog="mrs_almrec", verbose=verbose, opt=optParam, InputFormatisHealpix=False, OutputFormatisHealpix=True)
    return p

def tol(map,lmax_amin,amin=False):
    ns= gnside(map)
    lmax=lmax_amin
    if amin is True:
        lmax=amin2l(lmax_amin)
    a = mrs_almtrans(map, lmax=lmax)
    b = mrs_almrec(a, nside=ns)
    return  b

def mrs_uwttrans(map, lmax=None, opt=None, verbose=False, path='./',progpath=None):
    optParam = ' '
    if opt is not None:
        optParam = ' ' + opt
    if lmax is not None:
        optParam = ' -l ' + str(lmax) + optParam
    if progpath is None:
        prog="mrs_uwttrans"
    else:
        prog=progpath+"mrs_uwttrans"
    p = mrs_prog(map, prog=prog, verbose=verbose, opt=optParam, OutputFormatisHealpix=False,path=path)
    return p

def mrs_uwtrecons(Tmap, lmax=None, opt=None, verbose=False, path='./',progpath=None):
    optParam = ' '
    if opt is not None:
        optParam = ' ' + opt
    if lmax is not None:
        optParam = ' -l ' + str(lmax) + optParam
    if progpath is None:
        prog="mrs_uwttrans"
    else:
        prog=progpath+"mrs_uwttrans -r "
    p = mrs_prog(Tmap, prog=prog, verbose=verbose, opt=optParam, InputFormatisHealpix=False, OutputFormatisHealpix=True,path=path)
    return p
