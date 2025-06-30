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

    Example how to use the Class with 5  scales i.e. 4 wavelet scales + coarse resolution):
        CW = MRS_starlet()    # Create the class
        CW.init_starlet(Nside, nscale=5)
        CW.transform(Image)   # Starlet transform of a 2D np array
        CW.stat()             # print statistics of all scales
        r = CW.recons()       #  reconstruct an image from its coefficients
        CW.plot_filter()      # plot the filters in harmonic space  which are used in the wavelet decomposition
    more examples are given at the end of this file.
    
    Class variables are:
    nx = 0  # number of pixel of the healpix map
    ns = 0  # number of scales
    coef = 0.0  # Starlet coefficients
    TabNorm = 0.0  # Coefficient normalixation table
    SigmaNoise = 1.0  # noise standard deviation
    TabNsigma = 0  # detection level per scale
    verb = False
    nside = 0    # nside of the input map
    lmax = 0     # lmax used in spherical harmonic decomposition
    ALM_iter = 0 # numnber of iteration for the inverse spherical harmonic decomposition
    TabNameCode = ["Full python", "c++ Binding", "c++ binary"]
    TypeCode = 0  # 0 for 'Full python', '1' for 'c++ Binding' and 2 for'c++ binary'
    Tablmax =0  # lmax for each scale
    TabSigma =0 # Standard deviation of the Gaussian which fits the scaling function at every scale
    TabPhi = 0  # Scaling function for each scale
    TabPsi = 0 # Wavelet function for each scale
    Tabh = 0   # h filter  for each scale
    Tabg = 0   # g filter  for each scale
    TabResol = 0 # Resolution of each wavelet scale in arc minute
    PixelResol = 0 # pixel sizr in arc minute
    l2norm = False # if True, normlaize the coefficients (l2 normalization) such that the noise standart deviation remains constant through the scales. 
    NeedletFilter = False # If True, use needlet filters instead of spline filters
    
    Class functions are:
    def init_starlet(self, nside, nscale=0, lmax=0, ALM_iter=0, TabResolSigma=None, Needlet=None):
    def info(self):  # Print information relative to the intialisation.
    def stat(self):  # Print Min, Max, Mean and standard deviation of all scales.
    def plot_filter(self, wavelet=True, scaling=False, hfilter=False, gfilter=False):  # plot the filters which are used.
    def transform(self, data, WTname=None, opt=None): # Compute the wavelet transform
    WaveletScale = def get_scale(self, j): # returns the wavelet coefficients at a given scale : return  self.coef[j,:]
    def put_scale(self, ScaleCoef, j): # insert a scale in the wavelet transform:  self.coef[j,:] = ScaleCoef
    Rec =  recons(self): reconstruct an image from its wavelet coefficients
    DenoiseMap =  denoising(self, data, SigmaNoise=0, Nsigma=3, ThresCoarse=False, hard=True): # perform a denoising in wavelet space
    
    def tvs(self,j,min=None,max=None,title=None,sigma=None,lut=None,filename=None, dpi=100) : plot the scale j
    def tv(self, log=False, unit="", title="", minimum=None, maximum=None, cbar=True): plot all wavelet scales
    def dump(self) : print all variable values of the class
    SigmaNoise =   get_noise(self)  : estimate the noise standard deviation from the first wavelet scale
    TabNsigma = get_tabsigma(self, Nsigma=3) # Cretate the table for the detection level. Nsigma can be either a number of an array.
    def threshold(self,SigmaNoise=0,Nsigma=3,ThresCoarse=False,hard=True,FirstDetectScale=0,KillCoarse=False,Verbose=False): # Threshold the wavelet coefficient
    CopiedClass =  copy(self, name="wt"): # Return a copy of the class 
    def eval_computation_time(self): # Compare the computation time of different implementations of the wavelet transform
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
        Ns = 7
        ALM_iter = 0
        Needlet = False
        lmax = 2048
        C = CMRStarlet()
        C.init_starlet(Nside, nscale=Ns, ALM_iter=ALM_iter, lmax=lmax, Needlet=Needlet)
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
        info(d - r, "Denoising: Residual Standard Deviation: ")
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

    # Wavelet transform with non-dyadic scales
    C = CMRStarlet()
    TabResolSigma = np.array([4, 5, 7, 10, 12, 15, 20, 30, 60])
    # TabResolSigma =  np.array([7.45302012, 14.90604023,  29.81208047,  59.62416094 , 119.24832187,  238.49664375])
    C.init_starlet(Nside, TabResolSigma=TabResolSigma)

    n = np.random.normal(size=(Nside**2 * 12))
    C.transform(n)
    print("Check the theoretical norm estimation of Gaussian noise with std=1:")
    for j in range(C.ns):
        s = C.coef[j, :]
        print(
            "  Standard deviation of coefficient at scale ",
            j + 1,
            ": Sigma = ",
            s.std(),
        )
        print("                    Theoretical estimation  = ", C.TabNorm[j])

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


def sigma_filter(F_ell, nside, lmax=None, PixelWindow=False):
    """
    Calcule l'écart type théorique d'une carte HEALPix de bruit blanc gaussien
    convoluée par un filtre F_ell donné dans l'espace harmonique sphérique.

    Paramètres :
    - F_ell : array-like, filtre en harmonique sphérique (F_ell pour chaque ell)
    - sigma2_pixel : float, variance du bruit blanc par pixel (non convolué)
    - nside : int, résolution HEALPix de la carte

    Retour :
    - float, écart type théorique de la carte filtrée
    """
    sigma2_pixel = 1.0
    npix = 12 * nside**2
    if lmax is None:
        lmax = 3 * nside
    pixwin = hp.pixwin(nside, lmax=lmax)  # pixel window function
    # F_ell_eff = F_ell * pixwin

    ell = np.arange(len(F_ell))
    if PixelWindow:
        ell = ell * hp.pixwin(nside, lmax=len(F_ell) - 1)
    facteur = (2 * ell + 1) * np.abs(F_ell) ** 2
    variance_theorique = sigma2_pixel / npix * np.sum(facteur)
    return np.sqrt(variance_theorique)


MRS_StarletTabNorm = np.array(
    [0.969856, 0.103676, 0.051809, 0.025798, 0.012852, 0.006446, 0.003230, 0.001725]
)


def assert_strictly_increasing(table):
    for i in range(len(table) - 1):
        if table[i] >= table[i + 1]:
            raise ValueError(
                f"Table is not strictly increasing at index {i}: {table[i]} >= {table[i+1]}"
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
    Tablmax = 0  # lmax for each scale
    TabSigma = 0  # Standard deviation of the Gaussian which fits the scaling function at every scale
    TabPhi = 0  # Scaling function for each scale
    TabPsi = 0  # Wavelet function for each scale
    Tabh = 0  # h filter  for each scale
    Tabg = 0  # g filter  for each scale
    TabResol = 0  # Resolution of each wavelet scale in arc minute
    PixelResol = 0  # pixel sizr in arc minute
    l2norm = False
    NeedletFilter = False  # If True, use need filters instead of spline filters

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

    def init_starlet(
        self, nside, nscale=0, lmax=0, ALM_iter=0, TabResolSigma=None, Needlet=None
    ):
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
        if Needlet:
            self.NeedletFilter = True
            self.TypeCode = 0

        self.nside = np.int64(nside)
        self.nx = 12 * self.nside * self.nside
        if lmax != 0:
            self.lmax = np.int64(lmax)
        else:
            self.lmax = 3 * self.nside

        if TabResolSigma is not None:
            nscale = TabResolSigma.shape[0] + 1
            self.TabSigma = np.zeros(nscale)
            self.PixelResol = amin_pixel_size(self.nside)
            self.TabSigma[0] = splinelmax2sigma(self.lmax)
            if TabResolSigma.min() < self.TabSigma[0]:
                raise ValueError(
                    f"Smoothing Gaussian standard deviation must be larger than {self.TabSigma [0]:.2f}: {TabResolSigma.min():.2f} is not valid."
                )
            self.TabSigma[1:] = TabResolSigma
            assert_strictly_increasing(self.TabSigma)
            self.Tablmax = sigma2splinelmax(self.TabSigma, False)
        else:
            if nscale == 0:
                nscale = np.int64(np.log(self.nside) // 1) + 1
            self.TabSigma = np.zeros(nscale)

        self.ns = np.int64(nscale)

        if ALM_iter != 0:
            self.ALM_iter = np.int64(ALM_iter)

        if MRS_CXX and self.TypeCode == 1:
            CMRS = pymrs.MRS()
            CMRS.alloc(nside, self.ns, self.lmax, self.ALM_iter, self.verb)

        if Needlet is None:
            if TabResolSigma is None:
                (
                    self.Tablmax,
                    self.TabSigma,
                    self.TabPhi,
                    self.TabPsi,
                    self.Tabh,
                    self.Tabg,
                ) = get_default_filters(self.nside, self.ns)
                # print("Default TabSigma = ", self.TabSigma)
            else:
                self.TabPhi, self.TabPsi, self.Tabh, self.Tabg = get_sigmafilters(
                    self.TabSigma, self.lmax, Phi0Spline=False
                )
        else:
            filters = mrs_needlet_filters(self.lmax, NbrScale=self.ns)
            self.Tabh = filters["TabFilterH"]
            self.Tabg = filters["TabFilterG"]
            self.TabPhi = filters["TabPhi"]
            self.TabPsi = filters["TabPsi"]
            self.Tablmax = np.zeros((self.ns))
            lm = self.lmax
            for j in range(0, self.ns):
                self.Tablmax[j] = lm
                lm = lm / 2
            self.TabSigma = splinelmax2sigma(self.Tablmax)

        # print("ns = ", self.ns, self.TabSigma.shape )

        self.TabResol = np.zeros(self.ns)
        self.TabResol[0] = self.TabSigma[0]
        for j in range(1, self.ns - 1):
            self.TabResol[j] = (
                self.TabSigma[j - 1] + (self.TabSigma[j] - self.TabSigma[j - 1]) / 2.0
            )
        self.TabResol[self.ns - 1] = self.TabSigma[self.ns - 1]

        # print("TB ", self.TabResol)
        # Estimate the normalization factor for the coefficients.
        self.TabNorm = np.zeros(self.ns, dtype=np.float64)
        self.TabNorm[0] = 1.0
        PixelWindow = False
        Phi0Lmax = sigma_filter(
            self.TabPhi[:, 0], self.nside, lmax=self.lmax, PixelWindow=PixelWindow
        )
        DeltaPhi0 = np.sqrt(1.0 - Phi0Lmax**2)
        # the alm m around lmax are often not accurate, and 1 if  a conservative value, the
        # the correct value is most likely between 0.9 and 1.
        for j in range(0, self.ns):
            self.TabNorm[j] = sigma_filter(
                self.TabPsi[:, j], self.nside, lmax=self.lmax, PixelWindow=PixelWindow
            )
        if Needlet is None:
            self.TabNorm[0] = np.sqrt(
                DeltaPhi0**2
                + (
                    sigma_filter(
                        self.TabPsi[:, 0],
                        self.nside,
                        lmax=self.lmax,
                        PixelWindow=PixelWindow,
                    )
                )
                ** 2
            )

    def info(self):
        """
        Print information relative to the intialisation.
        """
        print(self.name, ": Npix = ", self.nx, ", Ns = ", self.ns)
        for j in range(self.ns - 1):
            print(
                "%s Band %2d, Scale %5.2f amin " % (self.name, j + 1, self.TabResol[j])
            )
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
                "%s Band %2d, Scale %5.2f amin: Min = %f, Max = %f, Mean = %f, std = %f"
                % (
                    self.name,
                    j + 1,
                    self.TabResol[j],
                    s.min(),
                    s.max(),
                    s.mean(),
                    s.std(),
                )
            )

    def plot_filter(self, wavelet=True, scaling=False, hfilter=False, gfilter=False):
        if wavelet:
            plot_filter(
                self.TabPsi,
                title="Spherical wavelet Filters",
                xlabel="Multipole l",
                ylabel="Amplitude",
                legend_prefix="Scale",
            )
        if scaling:
            plot_filter(
                self.TabPhi,
                title="Spherical Scaling Filters",
                xlabel="Multipole l",
                ylabel="Amplitude",
                legend_prefix="Scale",
            )
        if hfilter:
            plot_filter(
                self.Tabh,
                title="Spherical H-Filters",
                xlabel="Multipole l",
                ylabel="Amplitude",
                legend_prefix="Scale",
            )
        if gfilter:
            plot_filter(
                self.Tabg,
                title="Spherical G-Filters",
                xlabel="Multipole l",
                ylabel="Amplitude",
                legend_prefix="Scale",
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
            # print(self.ns, self.TabPhi.shape)
            if self.NeedletFilter is False:
                self.coef = wt_phi_filter_trans(im, self.TabPhi)
            else:
                self.coef = mrs_needlet_transform(im, self.TabPsi)
            # self.coef = mrs_uwttrans(im,self.ns,self.lmax,opt=None,verbose=self.verb,path="./",cxx=False)

        if self.l2norm is True:
            for j in range(self.ns):
                self.coef[j, :] = self.coef[j, :] / self.TabNorm[j]

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

        if self.l2norm is True:
            for j in range(self.ns):
                self.coef[j, :] = self.coef[j, :] * self.TabNorm[j]

        if self.NeedletFilter is False:
            rec = np.sum(self.coef, axis=0)
        else:
            rec = mrs_needlet_recons(self.coef, self.TabPsi)
        return rec

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

    def get_scale(self, j):
        """
        Return the scale j in self.coef
        Parameters
        ----------
        j : int
             Scale number. It must be in [0:self.ns].
        Returns
        -------
        None.

        """
        return self.coef[j, :]

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
        x.coef[:, :] = self.coef[:, :]
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
        print(f"==> Execution time full python code : {end - start:.4f} seconds")

        self.TypeCode = 1
        print("Use ", self.TabNameCode[self.TypeCode], " code:")

        start = time.time()
        w = self.transform(d)
        end = time.time()
        print(f"==> Execution time C++ binding : {end - start:.4f} seconds")

        self.TypeCode = 2
        print("Use ", self.TabNameCode[self.TypeCode], " code:")
        start = time.time()
        w = self.transform(d)
        end = time.time()
        print(f"==> Execution time C++ binary: {end - start:.4f} seconds")


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


# Define the spline function in l-space
def spline(l):
    return (1 / 12) * (
        np.abs(l - 2) ** 3
        - 4 * np.abs(l - 1) ** 3
        + 6 * np.abs(l) ** 3
        - 4 * np.abs(l + 1) ** 3
        + np.abs(l + 2) ** 3
    )


def spline2(size, l=1, lc=1):
    """
    Compute a non-negative decreasing spline, with value 1 at index 0.

    Parameters
    ----------
    size: int
        size of the spline
    l: float
        spline parameter. If l=1, the function goes to zero at l = size
                          If l=2, the function goes to zero at l = size / 2
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


def plot_filter(
    T,
    x=None,
    title="Spherical wavelet Filters",
    xlabel="l",
    ylabel="Amplitude",
    legend_prefix="Scale",
):
    """
    Plot a series of curves from a 2D array T.

    Parameters:
    - T: 2D numpy array of shape (n_points, n_curves)
    - x: Optional 1D array of length n_points for the x-axis values
    - title: Title of the plot
    - xlabel: Label for the x-axis
    - ylabel: Label for the y-axis
    - legend_prefix: Prefix for the legend entries
    """
    T = np.asarray(T)

    if x is None:
        x = np.arange(T.shape[0])

    if len(x) != T.shape[0]:
        raise ValueError("Length of x must match number of rows in T.")

    plt.figure(figsize=(10, 6))

    for i in range(T.shape[1]):
        z = T[:, i]
        plt.plot(x, z / z.max(), label=f"{legend_prefix} {i+1}")

    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()


def plotsig(
    T,
    x=None,
    title="Spherical wavelet Filters",
    xlabel="X",
    ylabel="T[:, i]",
    legend_prefix="Scale",
):
    """
    Plot a series of curves from a 2D array T.

    Parameters:
    - T: 2D numpy array of shape (n_points, n_curves)
    - x: Optional 1D array of length n_points for the x-axis values
    - title: Title of the plot
    - xlabel: Label for the x-axis
    - ylabel: Label for the y-axis
    - legend_prefix: Prefix for the legend entries
    """
    T = np.asarray(T)

    if x is None:
        x = np.arange(T.shape[0])

    if len(x) != T.shape[0]:
        raise ValueError("Length of x must match number of rows in T.")

    plt.figure(figsize=(10, 6))

    plt.plot(x, T, label=f"{legend_prefix}")

    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()


def ploth(
    lmax,
    x=None,
    title="h filter",
    xlabel="l",
    ylabel="T[:, i]",
    legend_prefix="H filter",
):
    """
    Plot a series of curves from a 2D array T.

    Parameters:
    - T: 2D numpy array of shape (n_points, n_curves)
    - x: Optional 1D array of length n_points for the x-axis values
    - title: Title of the plot
    - xlabel: Label for the x-axis
    - ylabel: Label for the y-axis
    - legend_prefix: Prefix for the legend entries
    """

    x = np.arange(0, lmax + 1)
    size = 1

    # tab2 = spline2(size, lc, 1)
    h = compute_h(lmax, 1)
    a = splinelmax2sigma(lmax, True) * 60.0 / 180.0 * np.pi
    g = gaussian_l(x, a / 60.0 / 180.0 * np.pi, 1)

    plt.figure(figsize=(10, 6))
    plt.plot(x, h, label="{legend_prefix}")
    plt.plot(x, g, "--", label=f"Est Gaussian fit\nσ °", lw=2)

    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()


def amin_pixel_size(nside):
    res_rad = hp.nside2resol(nside)  # en radians
    res_deg = np.degrees(res_rad)  # en degrés
    res_arcmin = res_deg * 60  # en arcminutes
    return res_arcmin


# Gaussian beam in l-space
def gaussian_l(l, sigma, A):
    return A * np.exp(-0.5 * l * (l + 1) * sigma**2)


def splinelmax2sigma(lmax, hfilter=False):
    sig3000 = 3.81594630  # (amin)) : gaussian fit done at lmax = 3000
    sig = 3000.0 / lmax * sig3000
    if hfilter:
        sig = sig * 1.771
    return sig


def sigma2splinelmax(sigma, hfilter):
    sig3000 = 3.81594630  # (amin)) : gaussian fit done at lmax = 3000
    lmax = 3000.0 / sigma * sig3000
    if hfilter:
        lmax = lmax * 1.771
    return lmax


# TabSigma=np.array([1.,2,3,4,8,12,20,30])
def safe_divide(a, b, eps=1e-8):
    b_safe = np.where(
        np.abs(b) > eps, b, np.inf
    )  # Replace near-zero with ∞ to avoid divide-by-zero
    result = a / b_safe
    result[np.abs(b) <= eps] = 0  # Explicitly set output to 0 where b was near-zero
    return result


def get_sigmafilters(TabSigma, lmax, Phi0Spline=False):
    nscales = TabSigma.shape[0]
    TabPhi = np.zeros((lmax + 1, nscales))
    TabPsi = np.zeros((lmax + 1, nscales))
    Tabh = np.zeros((lmax + 1, nscales))
    Tabg = np.zeros((lmax + 1, nscales)) + 1.0
    if Phi0Spline:
        TabPhi[:, 0] = spline2(lmax, 1)
    else:
        TabPhi[:, 0] = 1
    for j in range(1, nscales):
        lm = sigma2splinelmax(TabSigma[j], False)
        TabPhi[:, j] = spline2(lmax, lmax / lm)
        TabPsi[:, j - 1] = TabPhi[:, j - 1] - TabPhi[:, j]
        # print("Scale ", j+1, " Sigma = ",  TabSigma[j], ", lm = ", lm)
        Tabh[:, j - 1] = safe_divide(TabPhi[:, j], TabPhi[:, j - 1])
        Tabg[:, j - 1] = 1.0 - Tabh[:, j - 1]
    TabPsi[:, nscales - 1] = TabPhi[:, nscales - 1]
    Tabh[:, nscales - 1] = Tabh[:, nscales - 2]
    Tabg[:, nscales - 1] = Tabh[:, nscales - 2]
    return TabPhi, TabPsi, Tabh, Tabg


def get_lmaxfilters(Tablmax):
    nscales = Tablmax.shape[0]
    lmax = Tablmax[0]
    nl = Tablmax[0] + 1
    TabPhi = np.zeros((nl, nscales))
    TabPsi = np.zeros((nl, nscales))
    Tabh = np.zeros((nl, nscales))
    Tabg = np.zeros((nl, nscales)) + 1.0
    TabPhi[:, 0] = sigma2splinelmax(TabSigma[0], False)
    for j in range(1, nscales):
        lm = sigma2splinelmax(TabSigma[j - 1], False)
        # print("Scale ", j+1, " Sigma = ",  TabSigma[j-1], ", lm = ", lm)
        TabPhi[:, j] = spline2(lmax, lmax / lm)
        TabPsi[:, j - 1] = TabPhi[:, j - 1] - TabPhi[:, j]
        Tabh[:, j - 1] = TabPhi[:, j] / (TabPhi[:, j - 1] + 1e-6)
        Tabg[:, j - 1] = 1.0 - Tabh[:, j - 1]
    TabPsi[:, nscales - 1] = TabPhi[:, nscales - 1]
    return TabPhi, TabPsi, Tabh, Tabg


def get_default_filters(nside, nscales, lmax=None):
    if lmax is None:
        lmax = 3 * nside
    Tablmax = np.zeros((nscales))
    lm = lmax
    for j in range(0, nscales):
        Tablmax[j] = lm
        lm = lm / 2
    TabSigma = splinelmax2sigma(Tablmax)
    tl = sigma2splinelmax(TabSigma[0], True)
    TabPhi, TabPsi, Tabh, Tabg = get_sigmafilters(TabSigma, lmax, True)
    return Tablmax, TabSigma, TabPhi, TabPsi, Tabh, Tabg


def wt_hfilter_trans(Map, Tabh):
    nside = hp.get_nside(Map)
    npix = Map.shape[0]
    lmax = Tabh.shape[0] - 1
    nscales = Tabh.shape[1]
    wts = np.zeros((nscales, npix))
    h_scale = np.copy(Map)
    alms = map2alm(h_scale, lmax=lmax)
    AlmsHF = np.copy(alms)
    for j in range(nscales - 1):
        h = Tabh[:, j]
        AlmLF = alm_product(AlmsHF, h)
        # AlmLF = alm_product(alms, h)  # there a bug in the c++ code. Uncommenting this line reproduces the bug.
        l_scale = alm2map(AlmLF, nside)
        coef = h_scale - l_scale
        h_scale[:] = l_scale[:]
        AlmsHF[:] = AlmLF[:]
        wts[j, :] = coef
    wts[nscales - 1, :] = l_scale
    return wts


def wt_phi_filter_trans(Map, TabFilter):
    nside = hp.get_nside(Map)
    npix = Map.shape[0]
    lmax = TabFilter.shape[0] - 1
    nscales = TabFilter.shape[1]
    wts = np.zeros((nscales, npix))
    # print("WTS: ", nscales)
    h_scale = np.copy(Map)
    alms = map2alm(h_scale, lmax=lmax)
    for j in range(nscales - 1):
        h = TabFilter[:, j + 1]
        AlmLF = alm_product(alms, h)
        l_scale = alm2map(AlmLF, nside)
        wts[j, :] = h_scale - l_scale
        h_scale[:] = l_scale[:]
    wts[nscales - 1, :] = l_scale
    return wts


def convol_filter(Map, Filter):
    nside = hp.get_nside(Map)
    npix = Map.shape[0]
    lmax = Filter.shape[0] - 1
    alms = map2alm(Map, lmax=lmax)
    AlmCoef = alm_product(alms, Filter)
    return alm2map(AlmCoef, nside)


def test_wt_hfilter_trans():
    nside = 1024
    nscales = 5
    d = np.random.normal(size=(nside**2 * 12))
    Tablmax, TabSigma, TabPhi, TabPsi, Tabh, Tabg = get_default_filters(nside, nscales)

    w = wt_hfilter_trans(d, Tabh)
    w1 = wt_trans(d, nscales=4)
    return wts


def get_sigma_from_spline(lmax, hfilter=False):
    # Range of l values (spherical harmonics degrees)
    l_vals = np.arange(0, lmax + 1)
    # res = 2 * l * l_vals / lmax
    # s_vals = spline(res)
    if hfilter is False:
        s_vals = spline2(lmax, 1, 1)
    else:
        s_vals = compute_h(lmax, 1)

    # Fit Gaussian to spline
    popt, _ = curve_fit(gaussian_l, l_vals, s_vals, p0=(0.01, 1.0))
    sigma_est, A_est = popt

    # Print result
    print(
        f"Estimated angular σ ≈ {np.degrees(sigma_est):.8f} deg (≈ {sigma_est:.5f} rad)"
    )
    print(f"Estimated amplitude A ≈ {A_est:.4f}")
    print(
        f"Estimated angular σ ≈ {np.degrees(sigma_est)*60.:.8f} amin (≈ {sigma_est:.5f} rad)"
    )

    a = splinelmax2sigma(lmax, hfilter)
    print(f"Linear Estimated angular σ ≈ {a:.4f} amin")

    sigma_rad = sigma_est
    sigma_amin = np.degrees(sigma_rad) * 60.0

    # Compute FWHM
    fwhm_rad = np.sqrt(8 * np.log(2)) * sigma_rad
    fwhm_amin = np.degrees(fwhm_rad) * 60.0

    print(f"Sigma (rad): {sigma_rad:.6f}")
    print(f"FWHM (rad): {fwhm_rad:.6f}")
    print(f"FWHM (amin): {fwhm_amin:.3f}")

    # Plot
    plt.plot(l_vals, s_vals, label="Spline($\ell$)", lw=2)
    plt.plot(
        l_vals,
        gaussian_l(l_vals, *popt),
        "--",
        label=f"Gaussian fit\nσ ≈ {np.degrees(sigma_est):.2f}°",
        lw=2,
    )
    plt.plot(
        l_vals,
        gaussian_l(l_vals, a / 60.0 / 180.0 * np.pi, 1),
        "--",
        label=f"Est Gaussian fit\nσ ≈ {np.degrees(sigma_est):.2f}°",
        lw=2,
    )
    plt.xlabel("$\ell$")
    plt.ylabel("Amplitude")
    plt.legend()
    plt.grid(True)
    plt.title("Spline($\ell$) vs Gaussian Beam Fit")
    plt.show()


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

    t = get_wt_filters(lmax, nscales)
    t1 = np.copy(t)
    scale = 1
    alm_high = np.copy(alms)
    for j in range(nscales):
        # print("scale ", j+1)
        h = compute_h(lmax, scale)
        t1[:, j] = h
        d = h - t[:, j + 1]
        info(d)
        if not alm_out:
            alm_low = alm_product(alm_high, h)
            m = alm2map(alm_low, nside)
        else:
            m = alm_product(alm_high, h)
        h_scale = l_scale - m
        l_scale = m
        alm_high[:] = alm_low[:]
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

    testclass = 0
    if testclass:
        C = test_mrs_class()
