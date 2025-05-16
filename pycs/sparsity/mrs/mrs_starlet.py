"""
Created on May 13, 2025

@author: ean-Luc Starck

                    CLASS FUNCTION  class MRS_starlet()

    Allow to perform a starlet transform on the sphere, manipulate it (visualisation, thresholding, statistics, etc),
    and to reconstruct.
    If pysap is installed, then the pymrs module should be available and
    the code will used C++ binding for fast calculation.
    Otherwise full python code is used.

    Details of the starlet transform on the sphere can be found in
    J.L. Starck, F. Murtagh, and J. Fadili,
    Sparse Image and Signal Processing: Wavelets and
    Related Geometric Multiscale Analysis,
    Cambridge University Press, Cambridge (GB),  2016.

    Example how to use the Class:
        CW = MRS_starlet()    # Create the class
        CW.transform(Image)   # Starlet transform of a 2D np array
        CW.stat()             # print statistics of all scales
        r = CW.recons()       #  reconstruct an image from its coefficients
    more examples are given at the end of this file.
"""

import importlib.util

spec = importlib.util.find_spec("pymrs")
if spec is None:
    # print("pymrs is available at:", spec.origin)
    MRS_CXX = False
else:
    pymrs = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pymrs)
    MRS_CXX = True


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
from pycs.misc.stats import *
from pycs.misc.mr_prog import *
from pycs.sparsity.mrs.mrs_tools import *
import getpass
import time


def test_mrs_class(Init=None):
    if Init is None:
        if getpass.getuser() == "starck":
            os.chdir("/Users/starck/Main/python/data/wlens/DES_UNION")
            print("Current directory:", os.getcwd())
            d = mrs_read("mars_topo_mola_hpx_1024.fits")
            Nside = gnside(d)
            print("Nside = ", Nside)
        else:
            Nside = 1024
            d = np.random.normal(size=(Nside**2 * 12))
        Ns = 9
        ALM_iter = 0
        C = CMRStarlet()
        C.init_starlet(Nside, nscale=Ns, ALM_iter=ALM_iter)
        print("PYTHON Code computation time:")
        start = time.time()
        C.transform(d)
        end = time.time()
        print(f"Execution time python: {end - start:.4f} seconds")
    else:
        C = Init

    # Visualization of scale 3
    j = 3
    C.tvs(j, title="Scale 3 of Mars")

    # Visualization of all scales
    C.tv(title="Mars")

    # Stat of all scales
    C.stat()

    r = C.recons()
    if Init is None:
        info(d - r, "Reconstruction error: ")

    C.threshold(SigmaNoise=20)
    r = C.recons()
    if Init is None:
        info(d - r, "Reconstruction error: ")
    tvs(r, title="Denoised map")
    if Init is None:
        tvs(d - r, title="Estimated noise map")

    # Denoising routine
    # def denoising(self, Image, SigmaNoise=0, Nsigma=3, ThresCoarse=False, hard=True):
    r = C.denoising(d, Nsigma=[4, 3, 2], ThresCoarse=True)
    tvs(r, title="Denoised map")
    if Init is None:
        tvs(d - r, title="Estimated noise map")

    C.eval_computation_time()

    return C


def pycs_env_status():
    import os
    import subprocess

    # subprocess.run(['tcsh', '-c', 'source ~/.tcshrc; echo $PATH'], check=True)
    print(os.environ.get("DYLD_LIBRARY_PATH"))
    print(os.environ.get("DYLD_FALLBACK_LIBRARY_PATH"))
    print(os.environ.get("LD_LIBRARY_PATH"))
    print(os.environ.get("PYTHONPATH"))
    env = os.environ.copy()
    env["DYLD_LIBRARY_PATH"] = "/usr/lib"
    result = subprocess.run(["mr_filter"], capture_output=True, text=True, env=env)
    print("STDOUT:", result.stdout)
    print("STDERR:", result.stderr)
    print("Return code:", result.returncode)


MRS_StarletTabNorm = np.array(
    [0.969856, 0.103676, 0.051809, 0.025798, 0.012852, 0.006446, 0.003230, 0.001725]
)


class CMRStarlet:
    """
    Class for the starlet decomposition and reconstruction
    """

    name = "wt"  # name of the class
    nx = 0  # number of pixel of the healpix map
    ns = 0  # number of scales
    coef = 0.0  # Starlet coefficients
    TabNorm = 0.0  # Coefficient normalixation table
    SigmaNoise = 1.0  # noise standard deviation
    TabNsigma = 0  # detection level per scale
    verb = False
    nside = 0
    lmax = 0
    ALM_iter = 0
    TabNameCode = ["Full python", "c++ Binding", "c++ binary"]
    TypeCode = 0  # 0 for 'Full python', '1' for 'c++ Binding' and 2 for'c++ binary'

    # __init__ is the constructor
    def __init__(self, name="wt", verb=False):
        """
        Constructor

        Parameters
        ----------
        name : string, optional
            name of transform. Used when information is printed. The default is 'wt'.
        verb : bool, optional
        Returns
        -------
        None.

        """
        self.name = name  # self.name is an object variable
        self.verb = verb

    def init_starlet(self, nside, nscale=0, lmax=0, ALM_iter=0):
        """
        Initialize the scale for a given image size and a number of scales.
        Parameters
        ----------
        nx, ny : int
            image size.
        nscale : int, optional
            Number of wavelet scales. The default is 0.
            If it is 0, the numnber of scales is fixed to
                  log( MIN([nx,ny]))
        Returns
        -------
        None.
        """
        self.nside = np.int64(nside)
        self.nx = 12 * self.nside * self.nside
        if nscale == 0:
            nscale = np.int64(np.log(nx) // 1)
        self.ns = np.int64(nscale)
        if lmax != 0:
            self.lmax = np.int64(lmax)
        else:
            self.lmax = 3 * self.nside
        if ALM_iter != 0:
            self.ALM_iter = np.int64(ALM_iter)

        if MRS_CXX and self.TypeCode == 1:
            CMRS = pymrs.MRS()
            CMRS.alloc(nside, self.ns, self.lmax, self.ALM_iter, self.verb)

        Ne = len(MRS_StarletTabNorm)
        self.TabNorm = np.zeros(self.ns, dtype=np.float64)
        if Ne >= self.ns:
            self.TabNorm = MRS_StarletTabNorm[: self.ns]
        else:
            self.TabNorm[:Ne] = MRS_StarletTabNorm
            val = MRS_StarletTabNorm[-1]
            for i in range(Ne, self.ns):
                val /= 2
                self.TabNorm[i] = val
        # print("TabNorm = ", self.TabNorm)

    def info(self):  # sound is a method (a method is a function of an object)
        """
        Print information relative to the intialisation.
        """
        print(self.name, ": Npix = ", self.nx, ", Ns = ", self.ns)
        # print("Coef TabSize = ", np.shape(self.coef))

    def stat(self):
        """
        Print Min, Max, Mean and standard deviation of all scales.
        Parameters
        ----------
        None.

        Returns
        -------
        None.

        """
        print(self.name, ": Np = ", self.nx, ", Ns = ", self.ns)
        for j in range(self.ns):
            s = (self.coef)[j]
            print(
                "%s Scale %2d: Min = %f, Max = %f, Mean = %f, std = %f"
                % (self.name, j + 1, s.min(), s.max(), s.mean(), s.std())
            )

    def transform(self, data, WTname=None, opt=None):
        """
        Apply the starlet transform to image. Coeffients are stored in
        self.coef[:,:].  self.coef[j,:] is the wavelet scale at scale j.
        See class routine put_scale to manipulate
        the coefficients.
        Parameters
        ----------
        im : 1D np.ndarray
             input image..
        opt: string
            command line options for the binary c++ code mrs_uwttrans
        WTname : string, optional
            Name given to the decomposition. The default is None.
        Returns
        -------
        None.
        """
        im = data.astype(float)
        Nx = (im.shape[-2:])[0]
        if Nx != self.nx:
            raise ValueError(f"Error: expected {self.nx}, but got {Nx}")

        if self.ns <= 1 or self.nx != Nx:
            self.init_starlet(Nx, nscale=0)
        if WTname is not None:
            self.name = WTname

        if MRS_CXX and self.TypeCode == 1:
            CMRS = pymrs.MRS()
            CMRS.alloc(self.nside, self.ns, self.nside * 3, self.ALM_iter, True)
            self.coef = CMRS.uwt(im, self.ns)
        elif self.TypeCode == 2:
            self.coef = mrs_uwttrans(
                im, self.ns, self.lmax, opt=opt, verbose=self.verb, path="./", cxx=True
            )
        else:
            self.coef = mrs_uwttrans(
                im,
                self.ns,
                self.lmax,
                opt=None,
                verbose=self.verb,
                path="./",
                cxx=False,
            )

    def recons(self):
        """
        Reconstruct an image from its calculated starlet coefficients.
        Parameters
        ----------

        Returns
        -------
        rec : 1D np.ndarray
            Reconstructed image.
        """
        return np.sum(self.coef, axis=0)

    def denoising(self, data, SigmaNoise=0, Nsigma=3, ThresCoarse=False, hard=True):
        """
        Do a denoising of the input image, by taking the wavelet decomposition,
        thresholding it, and reconstructing the denoised image.
        Parameters
        ----------
        Image : 1D np.ndarray
            DESCRIPTION.
        SigmaNoise : float, optional
            Standard deviation of the noise. Default is 0.
        Nsigma: float, optional
            Detection level. Defautl is 3  (.e. 3 SigmaNoise).
        ThresCoarse : bool, optional
            IF true the coarsest scale is removed. The default is False.
        hard : bool, optional
            Type of threshold, true for hard thresholding and false
            for soft thresholding. The default is True.
        Returns
        -------
        1D np.ndarray
            Denoised image.
        """
        Image = data.astype(float)

        if SigmaNoise == 0:
            SigmaNoise = get_noise(Image)
        self.SigmaNoise = SigmaNoise
        self.transform(Image)
        self.threshold(
            SigmaNoise=SigmaNoise, Nsigma=Nsigma, ThresCoarse=ThresCoarse, hard=hard
        )
        return self.recons()

    def put_scale(self, ScaleCoef, j):
        """
        Replace the scale j in self.coef by the 2D array ScaleCoef.
        Parameters
        ----------
        ScaleCoef : 2D np.ndarray
            New coefficients at scale j to be inserted in the class.
        j : int
             Scale number. It must be in [0:self.ns].
        Returns
        -------
        None.

        """
        self.coef[j, :] = ScaleCoef

    def tvs(
        self,
        j,
        min=None,
        max=None,
        title=None,
        sigma=None,
        lut=None,
        filename=None,
        dpi=100,
    ):
        """
        Display the scale j
        Parameters
        ----------
        j : int
             Scale number. It must be in [0:self.ns].
        Returns
        -------
        Window appearing showing scale j.
        """
        s = self.coef[j]
        tvs(
            s,
            min=min,
            max=max,
            title=title,
            sigma=sigma,
            lut=lut,
            filename=filename,
            dpi=dpi,
        )

    def tv(self, log=False, unit="", title="", minimum=None, maximum=None, cbar=True):
        """
        Display the scale j
        Parameters
        ----------
        j : int
             Scale number. It must be in [0:self.ns].
        Returns
        -------
        Window appearing showing scale j.
        """
        mrs_tv(
            self.coef,
            log=log,
            unit=unit,
            title=title,
            minimum=minimum,
            maximum=maximum,
            cbar=cbar,
        )

    def dump(self):
        """
        Print all variable and function names of the class
        Returns
        -------
        None.

        """
        print(self.__dict__)

    def get_noise(self):
        """
        Estimate the noise in the data from the first wavelet scale
        Returns
        -------
        SigmaNoise : float
            estimated noise standard deviation.
        """
        s = self.coef[0, :]
        SigmaNoise = mad(s)
        return SigmaNoise

    def get_tabsigma(self, Nsigma=3):
        """
        Create the detection table TabNsigma[0:nsale-1], for diffent type
        of calling.
        By default, it is 4 at the finest scale at 3 at the others.
        If Nsigma is an array small than the number of scales, the last value
        of Nsigma is repeated.
        exemple of call:
            print(CLASS.get_tabsigma(4))              => array([4., 3., 3., 3.])
            print(CLASS.get_tabsigma(4, Nsigma=[3,4]) => array([3, 4., 4., 4.])

        Parameters
        ----------
        Nsigma : int or 1D np.ndarray, optional
            Detect level [per scale]. The default is [4,3,..,3]

        Returns
        -------
        TabNsigma : 1D np.ndarray
            Detection level per scale.
        """
        vssig = vsize(Nsigma)
        nscale = self.ns
        if vssig[0] == 0:
            TabNsigma = np.full(nscale, Nsigma)
            TabNsigma[0] += 1
        else:
            if vssig[1] < nscale:
                extension = np.full(nscale - vssig[1], Nsigma[-1])
                TabNsigma = np.concatenate([Nsigma, extension])
            else:
                TabNsigma = Nsigma[:nscale]
        return TabNsigma

    def threshold(
        self,
        SigmaNoise=0,
        Nsigma=3,
        ThresCoarse=False,
        hard=True,
        FirstDetectScale=0,
        KillCoarse=False,
        Verbose=False,
    ):
        """
        Apply a hard or a soft thresholding on the coefficients self.coef
        Parameters
        ----------
        SigmaNoise : float, optional
            Noise standard deviation. The default is 0.
            If it is 0, it will be automatically estimated from the first scale.
        Nsigma : 1D np.ndarray, optional
            Detect level [per scale]. The default is [4,3,..,3]
        ThresCoarse : bool, optional
            If true the coarsest scale is also thresholded. The default is False.
        hard : bool, optional
            IF true, apply hard thresholding, and soft-thresholding otherwise.
            The default is True.
        FirstDetectScale : int, optional
            Remove the first FirstDetectScale scales. The default is 0.
        KillCoarse :  bool, optional
            IF true the coarsest scale is removed. The default is False.
        Verbose : TYPE, optional
            DESCRIPTION. The default is False.
        Returns
        -------
        None.

        """

        if SigmaNoise == 0:
            SigmaNoise = self.get_noise()  # scalar
        self.SigmaNoise = SigmaNoise
        if Verbose:
            print("SigmaNoise = ", SigmaNoise, ", vsize(SigmaNoise) = ", vs)

        self.TabNsigma = self.get_tabsigma(Nsigma=Nsigma)
        if Verbose:
            print("TabNsigma = ", self.TabNsigma)

        # The noise level is obtained at each scale by multiplying by self.TabNorm
        Thres = SigmaNoise * self.TabNsigma * self.TabNorm

        Last = self.ns - 1
        if ThresCoarse:
            Last = self.ns

        for j in range(Last):
            s = self.coef[j, :]
            if hard:
                mrs_hard_thresholding(s, Thres[j])
            else:
                mrs_soft_thresholding(s, Thres[j])
            if Verbose:
                print(
                    "     scale ",
                    j + 1,
                    ", % of non zeros = ",
                    np.count_nonzero(s) * 100.0 / float(self.nx),
                )

            if Verbose:
                print(
                    "     scale ",
                    j + 1,
                    ", % of non zeros = ",
                    np.count_nonzero(s) * 100.0 / float(self.nx),
                )
            self.coef[j, :] = s

        if FirstDetectScale > 0:
            self.coef[:FirstDetectScale, :] = 0.0
        if KillCoarse:
            self.coef[-1, :] = 0.0

    def copy(self, name="wt"):
        """
        Duplicate the class, making copy of the coefficients.
        Parameters
        ----------
        name : TYPE, optional
            DESCRIPTION. The default is "wt".

        Returns
        -------
        NewClass : starlet2d
            Copy of the class.
        """
        x = self
        x.name = name
        x.coef = np.zeros((x.ns, x.nx))
        x.TabNorm = np.copy(self.TabNorm)
        return x

    def eval_computation_time(self):
        import time

        # TabNameCode = ['Full python', 'c++ Binding', 'c++ binary']
        Nside = 1024
        Ns = 9
        ALMiter = 0
        self.init_starlet(Nside, nscale=Ns, ALM_iter=ALMiter)

        # d = np.random.normal(scale=100., size=self.nx)
        os.chdir("/Users/starck/Main/python/data/wlens/DES_UNION")
        print("Current directory:", os.getcwd())
        di = mrs_read("mars_topo_mola_hpx_1024.fits")
        d = di.astype(float)
        print(d.dtype)
        self.TypeCode = (
            0  # 0 for 'Full python', '1' for 'c++ Binding' and 2 for'c++ binary'
        )
        print("Use ", self.TabNameCode[self.TypeCode], " code:")

        start = time.time()
        w = self.transform(d)
        end = time.time()
        print(f"==> Execution time : {end - start:.4f} seconds")

        self.TypeCode = 1
        print("Use ", self.TabNameCode[self.TypeCode], " code:")

        start = time.time()
        w = self.transform(d)
        end = time.time()
        print(f"==> Execution time : {end - start:.4f} seconds")

        self.TypeCode = 2
        print("Use ", self.TabNameCode[self.TypeCode], " code:")
        start = time.time()
        w = self.transform(d)
        end = time.time()
        print(f"==> Execution time : {end - start:.4f} seconds")


################################  END CLASS ######################


def mrs_starlet(map, nscale=None, lmax=None):
    nside = gnside(map)
    if nscale is None:
        Ns = np.int64(np.log2(nside) - 2)
    else:
        Ns = np.int64(nscale)

    npix = map.shape[0]
    w = wt_trans(map, lmax=lmax, nscales=Ns - 1)
    trans = w.T
    return trans


def mrs_istarlet(trans):
    r = np.sum(trans, axis=0)
    return r


def mrs_uwttrans(
    map,
    nscale=None,
    lmax=None,
    opt=None,
    verbose=False,
    path="./",
    progpath=None,
    cxx=False,
):
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
        w = wt_trans(map, lmax=lmax, nscales=Ns - 1)
        p = np.zeros((Ns, npix))
        for j in range(Ns):
            p[j, :] = w[:, j]

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

    res = np.arange(0, size + 1)
    res = 2 * l * res / (lc * size)
    res = (
        (3 / 2)
        * 1
        / 12
        * (
            abs(res - 2) ** 3
            - 4 * abs(res - 1) ** 3
            + 6 * abs(res) ** 3
            - 4 * abs(res + 1) ** 3
            + abs(res + 2) ** 3
        )
    )
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

    tab1 = spline2(size, 2 * lc, 1)
    tab2 = spline2(size, lc, 1)
    h = tab1 / (tab2 + 1e-6)
    h[np.int64(size / (2 * lc)) : size] = 0
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

    tab1 = spline2(size, 2 * lc, 1)
    tab2 = spline2(size, lc, 1)
    g = (tab2 - tab1) / (tab2 + 1e-6)
    g[np.int64(size / (2 * lc)) : size] = 1
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

    wt_filters = np.ones((lmax + 1, nscales + 1))
    wt_filters[:, 1:] = np.array(
        [compute_h(lmax, 2**scale) for scale in range(nscales)]
    ).T
    wt_filters[:, :nscales] -= wt_filters[:, 1 : (nscales + 1)]
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
            wts = np.zeros((npix, nscales + 1), dtype="complex")
        else:
            npix = np.shape(alms)[1]
            wts = np.zeros((np.shape(maps)[0], npix, nscales + 1), dtype="complex")

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


def mrs_tv(maps, log=False, unit="", title="", minimum=None, maximum=None, cbar=True):
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

        def f(x):
            return x

    else:

        def f(x):
            return np.log10(x - minimum + 1)

    for i in range(np.shape(maps)[0]):
        if title:
            tit = title + ": Scale " + str(i + 1)
        else:
            tit = "Scale " + str(i + 1)
        hp.mollview(
            f(maps[i, :]),
            fig=None,
            unit=unit,
            title=tit,
            min=f(minimum),
            max=f(maximum),
            cbar=cbar,
        )


if __name__ == "__main__":
    print("Main :)")

    testclass = 1
    if testclass:
        C = test_mrs_class()
