#! /usr/bin/env Python
"""
Created on April 9 2020

@authors: Kostas Themelis & Jean-Luc Starck & Austin Peel
"""

import numpy as np
from scipy import ndimage
from numpy import linalg as LA
from scipy.special import erf

from pycs.sparsity.sparse2d.starlet import *
from pycs.misc.cosmostat_init import *
from pycs.misc.mr_prog import *
from pycs.misc.utilHSS import *
from pycs.misc.im1d_tend import *
from pycs.misc.stats import *
from pycs.sparsity.sparse2d.dct import dct2d, idct2d
from pycs.sparsity.sparse2d.dct_inpainting import dct_inpainting
from pycs.misc.im_isospec import *


def get_ima_spectrum_map(Px, nx, ny):
    """
    Create an isotropic image from a power spectrum
        Ima[i+nx/2, j+ny/2] = Px[ sqrt(i^2 + j^j) ]

    Parameters
    ----------
    Px : : np.ndarray
        1D powspec.
    nx,ny : int
        image size to be created.

    Returns
    -------
    power_map : np.ndarray
        2D image.
    """
    Np = Px.shape[0]
    Px[0] = 0.0  # set the zero frequency to zero
    Px = np.append(Px, 0)  # set to 0 all frequencies above Np
    k1, k2 = np.meshgrid(
        np.arange(nx) - nx / 2.0, np.arange(ny) - ny / 2.0, indexing="ij"
    )
    ip = np.sqrt(k1**2 + k2**2).astype(int)  # map of frequency norms
    ip[ip > Np] = Np
    power_map = Px[ip]  # 2D power spectrum (isotropic)

    return power_map


class shear_data:
    """
    Class for input data, containing the shear components g1,g2, the covariance matrix,
    the theoretical convergence power spectrum.
    """

    g1 = 0  # shear 1st component

    def __init__(self):  # __init__ is the constructor
        self.g1 = 0

    g2 = 0  # shear 2nd component
    Ncov = 0  # diagonal noise cov mat of g = g1 + 1j g2,  of same size as g1 and g2
    # the noise cov mat relative to g1 alone is Ncov /2.   (same for g2)
    mask = 0  # mask
    ktr = 0  # true kappa (...for simulations)
    g1t = 0  # true g1
    g2t = 0  # true g2
    ps1d = 0  # theoretical convergence power spectrum
    nx = 0
    ny = 0

    # file names
    DIR_Input = 0  # dir input data
    g1_fn = 0  # g1 file name
    g2_fn = 0  # g2 file name
    ktr_fn = 0  # for sumulation only, true convergence map
    ps1d_fn = (
        0  # Convergence map 1D theoretical power spectrum used for Wiener filtering
    )
    ncov_fn = 0  # covariance filename

    def get_shear_noise(self, FillMask=False, Nrea=0, inpshape=None):
        """
        Return a noise realisation using the covariance matrix.
        If FillMask is True, the non observed area where the covariance is infinitate,
        will be filled with randon values with a variance value equals to the maximum
        variance in the observed area, i.e. where the mask is 1.

        Parameters
        ----------
        FillMask : TYPE, optional
            DESCRIPTION. The default is False.
        Nrea : int, optional


        Returns
        -------
        n1 : np.ndarray
            noise realisation for g1.
        n2 : np.ndarray
            noise realisation for g2.

        """
        Mat = np.sqrt(self.Ncov / 2.0)  # shape = (nx, ny)
        if FillMask == True:
            ind = np.where(self.mask == 1)
            MaxCov = np.max(Mat[ind])
            ind = np.where(self.mask == 0)
            Mat[ind] = MaxCov
        size = Mat.shape
        if Nrea > 0:
            size = (Nrea, *size)
        if inpshape is not None:
            size = (*inpshape, *size)
        n1 = np.random.normal(loc=0.0, scale=Mat, size=size)
        n2 = np.random.normal(loc=0.0, scale=Mat, size=size)
        return n1, n2


def massmap_get_rms_error(Res, TrueSol, Mask, sigma=0):
    if sigma > 0:
        Resi = smooth2d(TrueSol - Res, sigma) * Mask
        TS = smooth2d(TrueSol, sigma) * Mask
    else:
        Resi = (TrueSol - Res) * Mask
        TS = TrueSol * Mask
    ind = np.where(Mask != 0)
    Resi[ind] = Resi[ind] - np.mean(Resi[ind])
    TS[ind] = TS[ind] - np.mean(TS[ind])
    Err = LA.norm(Resi) / LA.norm(TrueSol * Mask) * 100.0
    return Err


class massmap2d:
    """Mass Mapping class
    This class contains the tools to reconstruct mass maps from shear measurements
    """

    kernel1 = 0  # internal variable for wiener filtering
    kernel2 = 0  # internal variable for wiener filtering
    nx = 0  # image size  (number of lines)
    ny = 0  # image size  (number of column)
    # WT=0          # Starlet wavelet Class defined in starlet.py
    Verbose = False  # Parameter to switch on/off print
    DEF_niter = 12  # Default number if iterations in the iterative methods.
    DEF_Nrea = 10  # Default number of realizations.
    DEF_Nsigma = 3.0  # Default detection level in wavelet space
    niter_debias = 0  # For space recovery using soft thresholding, a final
    # debiasing step could be useful. By default we don't any
    # debiasing.
    DEF_FirstDetectScale = 1  # default first detection scale in wavelet space.
    # very often, the noise is highly dominating the
    # the signal, and the first or the first scales
    # can be removed.
    # DEF_FirstDetectScale=2 => the  two finest scales
    # are removed
    WT_Sigma = 0  # Noise standard deviation in the wavelet space
    WT_ActiveCoef = 0  # Active wavelet coefficients
    SigmaNoise = 0  # Noise standard deviation in case of Gaussian white noise

    def __init__(self, name="mass"):  # __init__ is the constructor
        self.WT = starlet2d()  # Starlet wavelet Class defined in starlet.py

    def init_massmap(self, nx, ny, ns=0):
        """
        Initialize the class for a given image size and a number of scales ns
        to be used in the wavelet decomposition.
        If ns ==0, the number of scales is automatically calculated in the
        starlet initialization (see init_starlet, field self.WT.ns).

        Parameters
        ----------
        nx, ny : int
            Image size
        ns : int, optional
            Number of scales. The default is 0.

        Returns
        -------
        None.

        """
        self.nx = nx
        self.ny = ny
        self.WT = starlet2d(gen2=True, l2norm=True, bord=1, verb=False)
        self.WT.init_starlet(nx, ny, nscale=ns)
        self.WT.name = "WT-MassMap"
        k1, k2 = np.meshgrid(np.fft.fftfreq(ny), np.fft.fftfreq(nx))
        denom = k1 * k1 + k2 * k2
        denom[0, 0] = 1  # avoid division by 0
        self.kernel1 = (k1**2 - k2**2) / denom
        self.kernel2 = (2 * k1 * k2) / denom
        self.ker_kappa2gamma = (k1 + 1j * k2) ** 2 / denom
        self.ker_gamma2kappa = (k1 - 1j * k2) ** 2 / denom
        if self.Verbose:
            print(
                "Init Mass Mapping: Nx = ",
                nx,
                ", Ny = ",
                ny,
                ", Nscales = ",
                self.WT.ns,
            )

    def inpaint(self, kappa, mask, niter=DEF_niter):
        """
        Apply the sparse inpainting recovery technique to an image using the
        Discrete Cosine Transform.

        Parameters
        ----------
        kappa : np.ndarray
                Input data array
        mask : TYPE
            DESCRIPTION.
        niter : TYPE, optional
            DESCRIPTION. The default is self.DEF_niter.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        return dct_inpainting(kappa, mask, niter=niter, thresholding="soft")

    def get_theo_kappa_power_spectum(
        self, d, niter=None, PowSpecNoise=None, FirstFreqNoNoise=1
    ):
        """
        Estimate the theoretical convergence power spectrum from the data themselfve.
        Two methods are available:
            Method 1: Estimate inpainted ke and kb using iterative Kaiser-Squire method.
                    if PowSpecNoise==0, assume that there is no B-mode,
                    and the B-mode is used a noise power spectrum estimation.
                    powspec_Theo_E = powspec(ke) - powspec(kb)
                    powspec_Theo_B = 0
                    powspec_Theo_Noise = powspec(kb)
            Method 2: Use the input noise power spectrum
                    Then:
                        powspec_Theo_E = powspec(ke) - PowSpecNoise
                        powspec_Theo_B = powspec(kb) - PowSpecNoise

        Parameters
        ----------
        d :  Class  shear_data
            Input Class describing the obervations.
        niter : int, optional
            Number of iterations in the iKS method. The default is None.
        PowSpecNoise : np.ndarray, optional
            Noise power spectrum. The default is 0.
        FirstFreqNoNoise : int, optional
            At very low frequencies, signal is dominating and we generally prefer to not
            apply any denoising correction. So we will have :
                    powspec_Theo_Noise[0:FirstFreqNoNoise] = 0
            The default is 1.

        Returns
        -------
        pke : TYPE
            DESCRIPTION.
        pkb : TYPE
            DESCRIPTION.
        pn : TYPE
            DESCRIPTION.

        """
        k = self.iks(d.g1, d.g2, d.mask, niter=niter)
        ke = k.real
        kb = k.imag
        pe = im_isospec(ke)
        pb = im_isospec(kb)
        # fsky = mask.sum()/mask.size
        if PowSpecNoise is None:
            pn = pb
        else:
            pn = PowSpecNoise
        pn[0:FirstFreqNoNoise] = 0.0
        pke = pe - pn
        pkb = pb - pn

        # Denoise the estimated powsepc
        UseTendancyFiltering = False
        if UseTendancyFiltering is True:
            e1 = reverse(pke)
            fe1 = im1d_tend(e1)
            pke = reverse(fe1)
            b1 = reverse(pkb)
            fb1 = im1d_tend(b1)  # , opt='-T50'))
            pkb = reverse(fb1)
            pke[pke < 0] = 0
            # the smoothing does not work very well above nx/2,
            # because of the increase of the variance (frequencies at the corner)
            # we apodize the corner
            npix = pke.shape
            npix = npix[0]
            fp = int(npix / np.sqrt(2.0))
            min_end = pke[fp]
            pke[fp::] = pke[fp]
            pke[fp::] = min_end

        pkb[pkb < 0] = 0
        tv = 0
        if tv:
            plot(pke)
            plot(pkb)
            plot(pn)
            plot(d.ps1d)
        return pke, pkb, pn

    def get_tps(self, d, niter=None, Nrea=None):
        return self.get_theo_kappa_power_spectum(d, niter=niter)

    def kappa_to_gamma(self, kappa):
        """
        This routine performs computes the shear field from the  convergence
        map (no B-mode).

        Parameters
        ----------
        kappa: np.ndarray
                Input convergence data array

        Returns
        -------
        g1,g2: np.ndarray
                complext output shear field
        Notes
        -----
        """
        (Nx, Ny) = np.shape(kappa)[-2:]
        if self.nx != Nx or self.ny != Ny:
            self.init_massmap(Nx, Ny)
        out = np.fft.fft2(kappa)
        out = np.fft.ifft2(self.ker_kappa2gamma * out)
        return out.real, out.imag

    # Fast call
    def k2g(self, kappa):
        return self.kappa_to_gamma(kappa)

    def gamma_to_cf_kappa(self, g1, g2):
        """
        This routine performs a direct inversion from shear to convergence,
        it return a comlex field, with the real part being the convergence (E mode),
        the imaginary part being the B mode.

        Parameters
        ----------
        g1, g2: np.ndarray
                Input shear field

        Returns
        -------
        kappa: np.ndarray
                output complex convergence field
        Notes
        -----
        """
        if self.WT.nx == 0 or self.WT.ny == 0:
            raise NotImplementedError  # wrong arguments in self.WT.init_starlet
            nx, ny = np.shape(g1)[-2:]
            self.WT.init_starlet(nx, ny, gen2=1, l2norm=1, name="WT-MassMap")
        out = np.fft.fft2(g1 + 1j * g2)
        return np.fft.ifft2(self.ker_gamma2kappa * out)

    def gamma_to_kappa(self, g1, g2):
        """
        Same as gamma_to_cf_kappa, but returns only the E mode (convergence)

        Parameters
        ----------
        g1, g2: np.ndarray
                Input shear field

        Returns
        -------
        kappa: np.ndarray
                output convergence field
        Notes
        -----
        """
        k = self.gamma_to_cf_kappa(g1, g2)
        return k.real

    # Fast  interactive call to gamma_to_kappa
    def g2k(self, gam1, gam2):
        return self.gamma_to_kappa(gam1, gam2)

    def smooth(self, map, sigma=2.0):
        """
        Gaussian smoothing of an image.

        Parameters
        ----------
        map : 2D np.ndarray
        input image.
        sigma : float, optional
            Standard deviation of the used Gaussian kernel. The default is 2..

        Returns
        -------
         np.ndarray
            Smoother array.

        """
        return ndimage.gaussian_filter(map, sigma=sigma, axes=(-2, -1))

    def kaiser_squires(self, gam1, gam2, sigma=2.0):
        """
        This routine performs a direct inversion from shear to convergence,
        followed by a Gaussian filtering.
        This is the standard Kaiser-Squires method.

        Parameters
        ----------
        gam1, gam2: np.ndarray
                Input shear field
        sigma: float, optional
                Default is 2.
        Returns
        -------
        kappa: np.ndarray
                output convergence field
        Notes
        -----
        """
        ks = self.gamma_to_cf_kappa(gam1, gam2)
        ksg = ndimage.gaussian_filter(ks.real, sigma=sigma, axes=(-2, -1))
        return ksg

    # Fast interactive call to kaiser_squires
    def ks(self, gam1, gam2, sigma=2.0):
        return self.kaiser_squires(gam1, gam2, sigma=sigma)

    def eb_kaiser_squires(self, gam1, gam2, sigma=2.0):
        """
        Same as kaiser_squires, but return also the B-mnode.

        Parameters
        ----------
        gam1, gam2: np.ndarray
                Input shear field

        Returns
        -------
        E_kappa: np.ndarray
                output convergence field (E mode)
        B_kappa: np.ndarray
                output convergence field (B mode)
        Notes
        -----
        """
        ks = self.gamma_to_cf_kappa(gam1, gam2)
        ksg = ndimage.gaussian_filter(ks.real, sigma=sigma, axes=(-2, -1))
        ksbg = ndimage.gaussian_filter(ks.imag, sigma=sigma, axes=(-2, -1))
        return ksg, ksbg

    def H_operator_eb2g(self, ka_map, kb_map):
        """
        This routine converts (E,B) modes to shear

        Parameters
        ----------
        ka_map, kb_map : np.ndarray
            (E,B) mode

        Returns
        -------
        (g1,g2): np.ndarray
        output shear field
        None.
        """
        out = np.fft.fft2(ka_map + 1j * kb_map)
        out = np.fft.ifft2(self.ker_kappa2gamma * out)
        return out.real, out.imag

    # Fast interactice call to H_operator_eb2g
    def eb2g(self, ka_map, kb_map):
        return self.H_operator_eb2g(ka_map, kb_map)

    def H_adjoint_g2eb(self, g1_map, g2_map):
        """
        This routine reconstruct the (E,B) modes from the shear field

        Parameters
        ----------
        g1_map, g2_map : 2D np.ndarray or stack of 2D np.array
            shear field(s).

        Returns
        -------
        (E,B) modes : np.ndarray
            output convergence field
        None.
        """
        out = np.fft.fft2(g1_map + 1j * g2_map)
        out = np.fft.ifft2(self.ker_gamma2kappa * out)
        return out.real, out.imag

    # Fast interactice call to H_adjoint_g2eb
    def g2eb(self, g1_map, g2_map):
        return self.H_adjoint_g2eb(g1_map, g2_map)

    def get_wt_noise_level(self, InshearData, Nrea=DEF_Nrea):
        """
        Computes the noise standard deviation for each wavelet coefficients of
        the convergence map, using Nrea noise realisations of the shear field

        Parameters
        ----------
        InshearData : Class  shear_data
            Input Class describing the obervations.
        Nrea : int, optional
            Number of noise realisations. The default is 20.

        Returns
        -------
        WT_Sigma : 3D np.ndarray
            WT_Sigma[s,i,j] is the noise standard deviation at scale s and position
            (i,j) of the convergence.

        """
        n1, n2 = InshearData.get_shear_noise(FillMask=True, Nrea=Nrea)
        ke, _ = self.g2eb(n1, n2)
        self.WT.transform(ke)
        # Sum over noise realizations
        WT_Sigma = np.sum(
            self.WT.coef**2, axis=0
        )  # by definition the mean of wt is zero
        WT_Sigma = np.sqrt(WT_Sigma / Nrea)

        return WT_Sigma

    def get_active_wt_coef(
        self,
        InshearData,
        UseRea=False,
        SigmaNoise=1.0,
        Nsigma=None,
        Nrea=None,
        WT_Sigma=None,
        FirstDetectScale=None,
        OnlyPos=False,
        ComputeWTCoef=True,
    ):
        """
        Estimate the active set of coefficents, i.e. the coefficients of the
        convergence map with an absolute value large than Nsigma * NoiseStandardDeviation.
        It returns a cube  A[s,i,j] containing 0 or 1.
        If A[s,i,j] == 1 then we consider we have a detection at scale s and position (i,j).

        Parameters
        ----------
        InshearData : Class  shear_data
            Input Class describing the obervations.
        UseRea : bool, optional
            If true, make noise realisation to estimate the detection level in
            wavelet space. The default is False.
        Nrea : int, optional
            Number of noise realisations. The default is None.
        WT_Sigma : 3D np.ndarray, optional
            WT_Sigma[s,i,j] is the noise standard deviation at scale s and position
            (i,j) of the convergence. If it not given, the function get_wt_noise_level
            is used to calculate it.
        SigmaNoise: int, optional
            When UseRea==False, assume Gaussian nosie with standard deviation equals to SigmaNoise.
            Default is 1
        Nsigma : int, optional
            level of detection (Nsigma * noise_std). The default is None.
        FirstDetectScale: int, optional
            detect coefficients at scale < FirstDetectScale
        OnlyPos: Bool, optional
            Detect only positive wavelet coefficients. Default is no.
        ComputeWTCoef: bool, optional
            if true, recompute the wavelet coefficient  from the shear data.
            Default is true.
        Returns
        -------
        WT_Active : 3D np.ndarray
            WT_Active[s,i,j] = 1 if an coeff of the convergence map is larger
            than Nsigma * Noise std
        """

        if ComputeWTCoef:
            e, _ = self.g2eb(InshearData.g1, InshearData.g2)
            self.WT.transform(e)

        WT_Support = np.zeros_like(self.WT.coef)
        if UseRea and WT_Sigma is None:
            WT_Sigma = self.get_wt_noise_level(
                InshearData, Nrea=Nrea
            )  # shape = (ns, nx, ny)

        if Nsigma is None:
            Nsigma = self.DEF_Nsigma
        if Nrea is None:
            Nrea = self.DEF_Nrea

        wt = self.WT.coef  # shape = (ns, nx, ny) or (nimgs, ns, nx, ny)
        Nsig = Nsigma * np.ones(self.WT.ns)  # shape = (ns,)
        Nsig[0] += 1  # TODO: why?

        inp = np.abs(wt) if not OnlyPos else wt
        wsigma = WT_Sigma if UseRea else SigmaNoise
        WT_Support = (
            inp
            > wsigma
            * Nsig[:, np.newaxis, np.newaxis]
            * self.WT.TabNorm[:, np.newaxis, np.newaxis]
        ).astype(
            int
        )  # shape = ([nimgs], ns, nx, ny)

        if FirstDetectScale is not None:
            WT_Support[..., :FirstDetectScale, :, :] = 0
        WT_Support[..., -1, :, :] = 1  # TODO: why?

        return WT_Support

    def get_noise_powspec(self, CovMat, mask=None, nsimu=100, inpaint=False):
        """
        Build the noise powerspectum from the covariance map of the gamma field.
        Parameters
        ----------
        CovMat : : 2D np.ndarray
            covariance matrix of the shier field.
        mask : 2D np.ndarray, optional
            Apply a mask to the simulated noise realisation. The default is None.
        nsimu : int, optional
            Number of realisation to estimate the noise power spectrum. The default is 100.
        inpaint: Bool, optional
            Compute the power spectrum on inpainted Kaiser-Squires maps rather than on the masked
            maps. If inpaint==False, the estimated noise power spectrum is biased and
            should be corrected from .the fraction of sky not observed (i.e. fsky).
            Default is No
        Returns
        -------
        px : 1D np.ndarray
            Estimated Power spectrum from noise realizations.

        """
        if mask is None:
            m = 1.0
        else:
            m = mask
        for i in np.arange(nsimu):
            n1 = np.random.normal(loc=0.0, scale=np.sqrt(CovMat / 2.0)) * m
            n2 = np.random.normal(loc=0.0, scale=np.sqrt(CovMat / 2.0)) * m
            if mask is not None and inpaint is True:
                k = self.iks(n1, n2, mask)
            else:
                k = self.gamma_to_cf_kappa(n1, n2)
            p = im_isospec(k.real)

            if i == 0:
                Np = p.shape[0]
                TabP = np.zeros([nsimu, Np], dtype=float)
            TabP[i, :] = p
        px = np.mean(TabP, axis=0)
        return px

    def mult_wiener(self, map, WienerFilterMap):
        """ " apply one wiener step in the iterative wiener filtering"""
        # TODO: Modify `get_ima_spectrum_map` in order to avoid `np.fft.fftshift`
        # (should fix the bug with odd-sized arrays)
        return np.fft.ifft2(
            np.fft.fftshift(
                WienerFilterMap * np.fft.fftshift(np.fft.fft2(map), axes=(-2, -1)),
                axes=(-2, -1),
            )
        )

    def wiener(self, gamma1, gamma2, PowSpecSignal, PowSpecNoise):
        """
        Compute the standard wiener mass map.
        Parameters
        ----------
        gamma1,  gamma2: 2D np.ndarray
            shear fied.
        PowSpecSignal : 1D np.ndarray
            Signal theorical power spectrum.
        PowSpecNoise: 1D np.ndarray, optional
            noise theorical power spectrum.
        Returns
        -------
        TYPE  2D np.ndarray
              (E,B) reconstructed modes. Convergence = E
        """
        (nx, ny) = gamma1.shape

        if self.Verbose:
            print("Wiener filtering: ", nx, ny)

        Ps_map = get_ima_spectrum_map(PowSpecSignal, nx, ny)
        Pn_map = get_ima_spectrum_map(PowSpecNoise, nx, ny)
        Den = Ps_map + Pn_map
        ind = np.where(Den != 0)
        Wfc = np.zeros((nx, ny))
        Wfc[ind] = Ps_map[ind] / Den[ind]
        t = self.gamma_to_cf_kappa(gamma1, gamma2)  # xg + H^T(eta / Sn * (y- H * xg))
        kw = self.mult_wiener(t, Wfc)
        retr = np.zeros((nx, ny))
        reti = np.zeros((nx, ny))
        retr[:, :] = kw.real
        reti[:, :] = kw.imag

        return retr, reti

    def get_lmax_dct_inpaint(self, gamma1, gamma2):
        """return the maximum of the DCT absolute value of the convergence map"""
        eb = self.gamma_to_cf_kappa(gamma1, gamma2)
        lmax = np.max(np.abs(dct2d(eb.real, norm="ortho")), axis=(-2, -1))
        return lmax

    def step_dct_inpaint(
        self, xg, xn, mask, n, niter, lmin, lmax, InpaintAlsoImag=True
    ):
        """
        Apply one step of a iterative DCT inpainting method.
        First, we replace in xg value in mask[] == 0, obtained from the previous
        inpainting at iter n-1, and a hard thresholding in the DCT domain
        is applied to xg, on both real and imaginary parts.

        Parameters
        ----------
        xg :  2D np.cfarray
            convergence field.
        xn : 2D np.cfarray
            convergence field at the previous iteration.
        mask : 2D np.ndarray
            mask related to the observed data.
            mask[i,j] = 1 if the observed shear field has information
            at pixel (i,j)
        n : int
            iteration number in the interative method. n must be in [0,niter-1]
        niter : int
            number of iteration in the iterative method.
        lmin,max : int
            minimum and maximum absolute values of  the DCT transform
            of the shear E-mode map
        InpaintAlsoImag : bool, optional
            If true, both real and imaginary part are inpainted
            Otherwise on the real part is inpainted. The default is True.

        Returns
        -------
        xg : 2D np.cfarray
            inpainted convergence field.
        """
        lval = lmin + (lmax - lmin) * (1 - erf(2.8 * n / niter))  # exp decay
        if isinstance(lval, np.ndarray):
            lval = lval[..., np.newaxis, np.newaxis]  # broadcasting

        def _normalize(inp):
            # Enforce std. dev. constraint inside the mask
            inp_out = inp[
                ..., mask.astype(bool)
            ]  # shape = (p,) or (nimgs, p) or (nimgs, ns, p)
            inp_in = inp[
                ..., ~mask.astype(bool)
            ]  # shape = (p,) or (nimgs, p) or (nimgs, ns, p)
            std_out = inp_out.std(
                axis=-1
            )  # float or array of shape (nimgs,) or (nimgs, ns)
            std_in = inp_in.std(
                axis=-1
            )  # float or array of shape (nimgs,) or (nimgs, ns)
            if isinstance(std_in, np.ndarray):
                selectidx = (
                    std_in != 0
                )  # boolean array of shape (nimgs,) or (nimgs, ns)
                inp_in[selectidx, :] *= (std_out[selectidx] / std_in[selectidx])[
                    ..., np.newaxis
                ]  # shape = (q, p)
            else:
                if std_in != 0:
                    inp_in *= std_out / std_in

        def _step_dct_inpaint(ima, MultiScaleConstraint=False):

            alpha = dct2d(ima, norm="ortho")
            new_alpha = np.copy(alpha)  # Can we do this without copying ?
            new_alpha[np.abs(new_alpha) <= lval] = 0
            rec = idct2d(new_alpha, norm="ortho")  # shape = (nx, ny) or (nimgs, nx, ny)

            if not MultiScaleConstraint:
                _normalize(rec)
            else:  # it seems not improving the result
                self.WT.transform(rec)
                _normalize(self.WT.coef)
                rec = self.WT.recons()

            return rec

        # real part
        ima = mask * xg.real + (1 - mask) * xn.real
        rec = _step_dct_inpaint(ima, MultiScaleConstraint=False)

        # imaginary part
        if InpaintAlsoImag:
            ima = mask * xg.imag + (1 - mask) * xn.imag
            reci = _step_dct_inpaint(ima)
        else:
            reci = np.zeros_like(rec)

        return rec + 1j * reci

    def iks(self, g1, g2, mask, niter=None, dctmax=None):
        """
        Iterative Kaiser-Squires with DCT inpainting.

        Parameters
        ----------
        g1,g2 : np.ndarray
                Input convergence data array
         mask : np.ndarray
            mask  of missing data.
        niter : int, optional
            Number of iteration. The default is None.

        Returns
        -------
            2D complex np.ndarray
               output complex EB field
        """
        if niter is None:
            niter = self.DEF_niter
        lmin = 0
        assert g1.shape == g2.shape

        xg = np.zeros_like(g1, dtype=complex)  # TODO: complex or complex128?
        if dctmax is None:
            ks = self.gamma_to_cf_kappa(g1 * mask, g2 * mask)
            lmax = self.get_lmax_dct_inpaint(ks.real, ks.imag)
        else:
            lmax = dctmax

        for n in range(niter):
            t1, t2 = self.H_operator_eb2g(xg.real, xg.imag)
            r1 = mask * (g1 - t1)
            r2 = mask * (g2 - t2)
            t1, t2 = self.H_adjoint_g2eb(r1, r2)
            xg += t1 + 1j * t2  # xg + H^T(eta / Sn * (y- H * xg))
            xg = self.step_dct_inpaint(xg, xg, mask, n, niter, lmin, lmax)
        return xg

    def get_resi(
        self, xg, gamma1, gamma2, ResiWeight, mask=None, niter=None, dctmax=None
    ):
        """
        Compute the residual from an estimation of the convergence map.
        The return residual is on KappaE et KappaB.
        If a mask is given, then an inpainting iterative Kaiser-Squire
        method is used to backproject the residual from shear space to
        kappa space.

        Parameters
        ----------
        xg :  2D np.cfarray
            Estimated convergence field
        gamma1, gamma2 : 2D np.ndarray
            data: shear measurements.
        ResiWeight : 2D np.ndarray
            Weights to apply on shear residual.
        Mask:  2D np.ndarray, optional
            Weight to apply to the shear compoenents
            Default is none.
        niter: int, optional
            Number of iteration in the inpainting. Default is self.NITER
        dctmax: float, optional
            first threshold used in the inpainting. Default is automotically
            estimated.
        Returns
        -------
        resi_kappa_e : 2D np.ndarray
            residual E mode.
        resi_kappa_b : 2D np.ndarray
            residual B mode.

        """
        t1, t2 = self.H_operator_eb2g(xg.real, xg.imag)
        r1 = ResiWeight * (gamma1 - t1)  # shape = ([nimgs], [Nrea], nx, ny)
        r2 = ResiWeight * (gamma2 - t2)  # shape = ([nimgs], [Nrea], nx, ny)

        if mask is None:
            # H * xg
            r1, r2 = self.H_adjoint_g2eb(r1, r2)
        else:
            if niter is None:
                niter = self.DEF_niter
            # iterative Kaiser Squires with inpainting
            xi = self.iks(r1, r2, mask, niter=niter, dctmax=dctmax)
            r1 = xi.real
            r2 = xi.imag
        return r1, r2

    def _prepare_data(self, InshearData, msg=None, niter=None, Nsigma=None):
        gamma1 = InshearData.g1  # shape = ([nimgs], nx, ny)
        gamma2 = InshearData.g2  # shape = ([nimgs], nx, ny)

        if niter is None:
            niter = self.DEF_niter
        if Nsigma is None:
            Nsigma = self.DEF_Nsigma
        nx, ny = gamma1.shape[-2:]
        if self.Verbose:
            print(f"{msg}: ", nx, ny, ", Niter = ", niter)

        if not isinstance(InshearData.mask, np.ndarray):
            mask = (InshearData.Ncov != 0).astype(int)  # shape = (nx, ny)
        else:
            mask = InshearData.mask
        InshearData.Ncov[
            InshearData.Ncov == 0
        ] = 1e9  # infinite value for no measurement
        Ncv = InshearData.Ncov / 2.0  # shape = (nx, ny)

        # find the minimum noise variance
        ind = np.where(
            Ncv != 0
        )  # TODO: useless if we have set Ncv[mask == 0] = 1e9 before
        tau = np.min(Ncv[ind])

        # set the step size
        # eta = 1.83 * tau
        eta = tau
        # compute signal coefficient
        Esn = eta / Ncv  # shape = (nx, ny)
        Esn[
            Esn == np.inf
        ] = 0  # TODO: useless if we have set Ncv[mask == 0] = 1e9 before

        return gamma1, gamma2, nx, ny, eta, Esn, mask, ind, tau, niter, Nsigma

    def _get_Wfc(self, PowSpecSignal, nx, ny, Pn, eta):

        # calculate the wiener filter coefficients
        Px_map = get_ima_spectrum_map(PowSpecSignal, nx, ny)
        # info((Px_map + eta))
        if Pn is not None:
            Pn_map = get_ima_spectrum_map(Pn, nx, ny)
            Den = Px_map + Pn_map
            ind = np.where(Den == 0)
            Den[ind] = eta
            Wfc = Px_map / Den
        else:
            Wfc = Px_map / (Px_map + eta)
        Wfc[Wfc == np.inf] = 0

        return Wfc

    def _noise_realizations(self, InshearData, mask, **kwargs):
        n1, n2 = InshearData.get_shear_noise(**kwargs)
        gamma1 = n1 * mask  # shape = ([nimgs], [Nrea], nx, ny)
        gamma2 = n2 * mask  # shape = ([nimgs], [Nrea], nx, ny)
        return gamma1, gamma2

    def prox_wiener_filtering(
        self,
        InshearData,
        PowSpecSignal,
        Pn=None,
        niter=None,
        Inpaint=False,
        ktr=None,
        PropagateNoise=False,
        Nrea=None,
    ):
        """
        Compute the wiener mass map considering not stationary noise
        Proximal wiener method published in:
            J. Bobin, J.-L. Starck, F. Sureau, and J. Fadili,
            "CMB map restoration", Advances in Astronomy , 2012, Id703217, 2012.
        Parameters
        ----------
        InshearData : Shear Class
            Class contains the shear information.
        PowSpecSignal : 1D np.ndarray
            Signal theorical power spectrum.
        Pn: 1D np.ndarray, optional
            noise theorical power spectrum.
        niter: int
            number of iterations. Default is DEF_niter
        Inpaint: bool, optional
            if true, inpainting the missing data. Default is false.
        ktr: 2D np.ndarray, optional
            true convergence map, known in case simulated data are used.
            if given, errors are calculated at each iteration.
        PropagateNoise: Bool, optional
            if True, run the routine on a noise realization instead of the input shear field.
        Nrea: int, optional
            number of noise realizations
        Returns
        -------
        TYPE  2D np.ndarray
              (E,B) reconstructed modes. Convergence = E
        """
        gamma1, gamma2, nx, ny, eta, Esn, mask, ind, tau, niter, _ = self._prepare_data(
            InshearData, msg="Iterative Wiener filtering", niter=niter
        )

        # calculate the wiener filter coefficients
        Wfc = self._get_Wfc(PowSpecSignal, nx, ny, Pn, eta)

        if Inpaint:
            # TODO: to be placed before or after "if PropagateNoise"? Inconsistent between methods.
            lmin = 0
            lmax = self.get_lmax_dct_inpaint(gamma1, gamma2)

        if PropagateNoise:
            # Linear operator: uncertainty intervals do not depend on the input images
            gamma1, gamma2 = self._noise_realizations(
                InshearData, mask, Nrea=Nrea
            )  # shape = (Nrea, nx, ny)

        xg = np.zeros_like(gamma1)

        for n in range(niter):
            xn = np.copy(xg)
            t1, t2 = self.get_resi(xg, gamma1, gamma2, Esn)

            t = xg + (t1 + 1j * t2)  # xg + H^T(eta / Sn * (y- H * xg))
            xg = self.mult_wiener(t, Wfc)  # wiener filtering in fourier space
            if Inpaint:
                xg = self.step_dct_inpaint(xg, xn, mask, n, niter, lmin, lmax)

            if self.Verbose:
                if ktr is not None:
                    Err = massmap_get_rms_error(xg.real, ktr, mask, sigma=0)
                    print("   it. Wiener Iter ", n + 1, ", Err = ", Err)
                else:
                    print(
                        "   Wiener rec Iter: ",
                        n + 1,
                        ", std ke =  %5.4f" % (np.std(xg[..., ind[0], ind[1]] / tau)),
                    )

        return xg.real, xg.imag

    def test(self):
        print("hello")

    def prox_mse(
        self,
        InshearData,
        niter=None,
        Inpaint=True,
        sigma=None,
        ktr=None,
        PropagateNoise=False,
        Nrea=None,
    ):
        """
        Compute the Mean Square Error Estimator of the mass map
        considering not stationary noise, with an inpainting of the missing data.

        Parameters
        ----------
        InshearData : Shear Class
            Class contains the shear information.
        niter: int
            number of iterations. Default is DEF_niter
        Inpaint: bool, optional
            if true, inpainting the missing data. Default is True.
        sigma: float
            if set, a regularizing smoothing is applied to the solution.
            Default is no.
        ktr: 2D np.ndarray, optional
            true convergence map, known in case simulated data are used.
            if given, errors are calculated at each iteration.
        PropagateNoise: Bool, optional
            if True, run the routine on a noise realization instead of the input shear field.
        Nrea: int, optional
            number of noise realizations
        Returns
        -------
        TYPE  2D np.ndarray
              (E,B) reconstructed modes. Convergence = E
        """
        gamma1, gamma2, _, _, _, Esn, mask, ind, tau, niter, _ = self._prepare_data(
            InshearData, msg="Proxinal MSE estimator", niter=niter
        )

        if PropagateNoise:
            # Linear operator: uncertainty intervals do not depend on the input images
            gamma1, gamma2 = self._noise_realizations(
                InshearData, mask, Nrea=Nrea
            )  # shape = ([Nrea], nx, ny)

        xg = self.gamma_to_cf_kappa(gamma1, gamma2)  # shape = ([nimgs], nx, ny)
        if Inpaint:
            # TODO: check code: self.get_resi should be computed on shear maps, not convergence maps.
            # TODO: to be placed before or after "if PropagateNoise"? Inconsistent between methods.
            lmin = 0
            lmax = self.get_lmax_dct_inpaint(xg.real, xg.imag)

        for n in range(niter):
            t1, t2 = self.get_resi(xg, gamma1, gamma2, Esn)  # shape = ([nimgs], nx, ny)

            xg = xg + (t1 + 1j * t2)  # xg + H^T(eta / Sn * (y- H * xg))

            if sigma is not None:
                ksg = ndimage.gaussian_filter(xg.real, sigma=sigma, axes=(-2, -1))
                ksbg = ndimage.gaussian_filter(xg.imag, sigma=sigma, axes=(-2, -1))
                xg = ksg + 1j * ksbg

            if Inpaint:
                xg = self.step_dct_inpaint(xg, xg, mask, n, niter, lmin, lmax)

            if self.Verbose:
                if ktr is not None:
                    Err = massmap_get_rms_error(xg.real, ktr, mask, sigma=0)
                    print("   Prox MSE Iter ", n + 1, ", Err = ", Err)
                else:
                    print(
                        "   Prox MSE rec Iter: ",
                        n + 1,
                        ", std ke =  %5.4f" % (np.std(xg[..., ind[0], ind[1]] / tau)),
                    )

        return xg.real, xg.imag

    def sparse_wiener_filtering(
        self,
        InshearData,
        PowSpecSignal,
        niter=None,
        Nsigma=None,
        Inpaint=False,
        InpNiter=20,
        OnlyPos=True,
        FirstDetectScale=DEF_FirstDetectScale,
        Bmode=True,
        ktr=None,
        PropagateNoise=False,
        Nrea=None,
    ):
        """
        MCAlens algorithm; Estimate the complex EB mode. The solution is assumed to have
        a Gaussian component and a sparse Component.
        The algorithm estimate these two parts, for both the E and B mode.
        It returns 4 maps, estimated E and B mode and sparse E and B mode.
        The Gaussian E and B modes are obtained taking the difference between
        the estimated mode and sparse mode.
        Parameters
        ----------
        InshearData : Shear Class
            Class contains the shear information.
        PowSpecSignal : 1D np.ndarray
             Theorical Signal power spectrum.
        niter : int, optional
            number of iterations. Default is self.DEF_niter
        Nsigma : float, optional
            Detection level on wavelet coefficients. The default is self.DEF_Nsigma.
        Inpaint : Bool, optional
            if true, inpainting the missing data. Default is false.
        InpNiter : int, optional
            Number of iterationd in inpainting algorithm. The default is 20.
        OnlyPos : Bool, optional
            Only positive wavelet coefficients are detected. The default is True.
        FirstDetectScale : TYPE, optional
             No wavelet coefficient are detected in the finest wavelet scales.
        Bmode : Bool, optional
            Calculate also the B-mode. The default is True.
        ktr : 2D np.ndarray, optional
            true convergence map, known in case simulated data are used.
            if given, errors are calculated at each iteration.
        PropagateNoise : Bool, optional
            if True, run the routine on a noise realization instead of the input shear field.
        Nrea: int, optional
            number of noise realizations
        Returns
        -------
        2D np.ndarray
              E reconstructed mode. Convergence = E
        2D np.ndarray
              B reconstructed mode.
        2D np.ndarray
              E reconstructed mode of the sparse component. Convergence = E
        2D np.ndarray
              B reconstructed mode  of the sparse component.
        """
        (
            gamma1,
            gamma2,
            nx,
            ny,
            eta,
            Esn,
            mask,
            ind,
            tau,
            niter,
            Nsigma,
        ) = self._prepare_data(
            InshearData, msg="MCALens estimator", niter=niter, Nsigma=Nsigma
        )

        RMS_ShearMap = np.sqrt(InshearData.Ncov / 2.0)  # shape = (nx, ny)
        SigmaNoise = np.min(RMS_ShearMap)  # float
        Esn_Sparse = SigmaNoise / RMS_ShearMap  # shape = (nx, ny)
        Esn_Sparse[Esn_Sparse == np.inf] = 0

        # calculate the wiener filter coefficients
        Wfc = self._get_Wfc(PowSpecSignal, nx, ny, None, eta)  # shape = (nx, ny)

        # Detection of the significant wavelet coefficents.
        # to avoid border artefacts, we first make a rough very smooth estimate
        # using a smoothed KS, compute the residual and the residul to the
        # rough estimate. This avoid the detection of many wavelet coeff along
        # the border.

        if Inpaint:
            # TODO: check code: self.get_resi should be computed on shear maps, not convergence maps.
            # TODO: to be placed before or after "if PropagateNoise is True"? Inconsistent between methods.
            lmin = 0
            resi1, resi2 = self.get_resi(xg, gamma1, gamma2, Esn)
            lmax = self.get_lmax_dct_inpaint(resi1, resi2)

        ks = self.gamma_to_cf_kappa(gamma1, gamma2)
        rec = ks.real
        ks = self.smooth(rec, sigma=15)
        resi1, resi2 = self.get_resi(ks, gamma1, gamma2, Esn_Sparse)

        resi1 += ks
        self.WT.transform(resi1)
        self.WT_ActiveCoef = self.get_active_wt_coef(
            InshearData,
            OnlyPos=OnlyPos,
            UseRea=False,
            SigmaNoise=SigmaNoise,
            Nsigma=Nsigma,
            ComputeWTCoef=False,
        )  # shape = ([nimgs], ns, nx, ny)
        self.WT_ActiveCoef[..., -1, :, :] = 0

        # Replace the shear measurements by noise realisations
        if PropagateNoise:
            # Nonlinear operator: uncertainty intervals depend on the input image
            # Each input image gets its own set of noise realizations
            if Nrea > 0:
                self.WT_ActiveCoef = self.WT_ActiveCoef[
                    ..., np.newaxis, :, :, :
                ]  # shape = ([nimgs], 1, ns, nx, ny)
            inpshape = gamma1.shape[:-2]  # typically, inpshape = (nimgs,)
            gamma1, gamma2 = self._noise_realizations(
                InshearData, mask, Nrea=Nrea, inpshape=inpshape
            )  # shape = ([nimgs], [Nrea], nx, ny)

        # shape = ([nimgs], [Nrea], nx, ny)
        # TODO: complex or complex128? Same question for real-valued arrays
        xg = np.zeros_like(gamma1, dtype=complex)  # Gaussian + sparse components
        xs = np.zeros_like(gamma1, dtype=complex)  # sparse component
        xw = np.zeros_like(gamma1, dtype=complex)  # Gaussian component

        for n in range(niter):
            resi1, resi2 = self.get_resi(
                xg, gamma1, gamma2, Esn_Sparse
            )  # shape = ([nimgs], [Nrea], nx, ny)

            # sparse component
            xt = resi1 + 1j * resi2  # xg + H^T(eta / Sn * (y- H * xg))
            self.WT.transform(xt.real)
            self.WT.coef *= self.WT_ActiveCoef  # shape = ([nimgs], [Nrea], ns, nx, ny)
            signif_resi = self.WT.recons()  # shape = ([nimgs], [Nrea], nx, ny)
            rec = xs.real + signif_resi  # shape = ([nimgs], [Nrea], nx, ny)

            if PropagateNoise is False:
                if OnlyPos:
                    rec[rec < 0] = 0
            if Bmode:
                self.WT.transform(xs.imag)
                self.WT.threshold(
                    SigmaNoise=SigmaNoise,
                    Nsigma=Nsigma,
                    ThresCoarse=True,
                    hard=True,
                    FirstDetectScale=FirstDetectScale,
                    Verbose=False,
                )
                reci = self.WT.recons()  # shape = ([nimgs], [Nrea], nx, ny)
            else:
                reci = 0
            xs = rec + 1j * reci  # shape = ([nimgs], [Nrea], nx, ny)
            xg = xw + xs

            InpMethod1 = 1
            if InpMethod1:
                nw = 1
            else:
                nw = 1
            for i in range(nw):
                if Inpaint and InpMethod1:
                    # TODO: why not using self.step_dct_inpaint as in the other methods?
                    t1, t2 = self.get_resi(
                        xg, gamma1, gamma2, Esn, mask=mask, niter=InpNiter, dctmax=lmax
                    )
                else:
                    t1, t2 = self.get_resi(xg, gamma1, gamma2, Esn)
                xt = t1 + 1j * t2  # shape = ([nimgs], [Nrea], nx, ny)
                xt += xw  # xg + H^T(eta / Sn * (y- H * xg))
                xw = self.mult_wiener(
                    xt, Wfc
                )  # wiener filtering in fourier space; shape = ([nimgs], [Nrea], nx, ny)
                xg = xw + xs  # shape = ([nimgs], [Nrea], nx, ny)

            ZeroMeanCst = False
            if ZeroMeanCst is True:
                xgr = xg.real
                xgi = xg.imag
                mxg = np.mean(
                    xgr, axis=(-2, -1), keepdims=True
                )  # shape = ([nimgs], [Nrea], 1, 1)
                # mxgm = np.mean(xgr[mask == 1], axis=(-1), keepdims=True) # shape = ([nimgs], [Nrea], 1)
                xgr = xgr - mxg  # shape = ([nimgs], [Nrea], nx, ny)
                xg = xgr + 1j * xgi  # shape = ([nimgs], [Nrea], nx, ny)

            if self.Verbose:
                ind = np.where(mask == 1)
                if ktr is not None:
                    print(
                        "   Sparse rec Iter: ",
                        n + 1,
                        ", Err = %5.4f"
                        % (
                            np.std(
                                (
                                    xg.real[..., ind[0], ind[1]]
                                    - ktr[..., ind[0], ind[1]]
                                )
                            )
                            / np.std(ktr[..., ind[0], ind[1]])
                            * 100.0
                        ),
                        ", Resi ke (x100) =  %5.4f"
                        % (np.std(resi1[..., ind[0], ind[1]]) * 100.0),
                        ", Resi kb (x100) =  %5.4f"
                        % (np.std(resi2[..., ind[0], ind[1]]) * 100.0),
                    )
                else:
                    print(
                        "   Sparse rec Iter: ",
                        n + 1,
                        ", Resi ke =  %5.4f"
                        % (np.std(resi1[..., ind[0], ind[1]] / tau)),
                        ", Resi kb = %5.4f"
                        % (np.std(resi2[..., ind[0], ind[1]]) / tau),
                    )

        return xg.real, xg.imag, xs.real, xs.imag

    def step_wt_recons(self, xg):
        self.WT.transform(xg.real)
        self.WT.coef *= self.WT_ActiveCoef
        xg.real = self.WT.recons()
        return xg

    def sparse_recons(
        self,
        InshearData,
        UseNoiseRea=False,
        WT_Support=None,
        WT_Sigma=None,
        niter=None,
        Nsigma=None,
        ThresCoarse=False,
        Inpaint=False,
        ktr=None,
        FirstDetectScale=None,
        Nrea=None,
        hard=True,
        FirstGuess=None,
    ):
        """
        Reconstruction of the convergence field using sparsity. The detection levels

        convergence map with an absolute value large than Nsigma * NoiseStandardDeviation.
        It returns a cube  A[s,i,j] containing 0 or 1.
        If A[s,i,j] == 1 then we consider we have a detection at scale s and position (i,j).

        Parameters
        ----------
        InshearData : Class  shear_data
            Input Class describing the obervations.

        Parameters
        ----------
        InshearData : TYPE
            DESCRIPTION.
        WT_Support : 3D np.ndarray, optional
            This variable is not used anymore.
            We keep it for back compatibility,
            The default is None.
        WT_Sigma : 3D np.ndarray, optional
            Noise standard deviation on each wavelet coefficient. The default is None
            and it is calculated in the routine.
        niter : int, optional
            number of iterations. Default is DEF_niter
        Nsigma : float, optional
            Detection level on wavelet coefficients. The default is self.DEF_Nsigma.
        ThresCoarse : Bool, optional
            If True the coarsest wavelet scale is removed. The default is False.
        Inpaint : Bool, optional
            if true, inpainting the missing data. Default is False.
        ktr : 2D np.ndarray, optional
            true convergence map, known in case simulated data are used.
            Errors are calculated at each iteration.  Default is None.
        FirstDetectScale : int, optional
             No wavelet coefficient are detected in the finest wavelet scales.
             The default is DEF_FirstDetectScale.
        Nrea : int, optional
            Generate Nrea noise realizations to estimate the noise level on
            each wavelet coefficient. The default is DEF_Nrea.

        Returns
        -------
        2D np.ndarray
              E reconstructed mode. Convergence = E
        2D np.ndarray
              B reconstructed mode.
        """
        if niter is None:
            niter = self.DEF_niter
        if Nsigma is None:
            Nsig = self.DEF_Nsigma
        else:
            Nsig = Nsigma
        if not hard:
            Nsig = Nsig / 2.0  # for soft thresholding, the thresholding must
            # be small than for hard thresholding
        if FirstDetectScale is None:
            FirstDetectScale = self.DEF_FirstDetectScale
        if Nrea is None:
            Nrea = self.DEF_Nrea
        (nx, ny) = InshearData.g1.shape

        if UseNoiseRea:
            if WT_Sigma is None:
                self.WT_Sigma = self.get_wt_noise_level(InshearData, Nrea=Nrea)
            else:
                self.WT_Sigma = WT_Sigma
            WeightResi = InshearData.mask
        else:
            RMS_ShearMap = np.sqrt(InshearData.Ncov / 2.0)
            tau = np.min(RMS_ShearMap)
            self.SigmaNoise = tau

            ind = np.where(RMS_ShearMap > 0)
            WeightResi = np.zeros((nx, ny))
            WeightResi[ind] = tau / RMS_ShearMap[ind]

            self.WT_Sigma = tau
        gamma1 = InshearData.g1
        gamma2 = InshearData.g2
        mask = InshearData.mask
        self.nx = nx
        self.ny = ny

        Iz = 1j * np.zeros((nx, ny))
        xg = np.zeros((nx, ny)) + 1j * np.zeros((nx, ny))
        if FirstGuess is not None:
            if FirstGuess.dtype == "float64":
                xg = FirstGuess + Iz
            elif FirstGuess.dtype == "complex128":
                xg = FirstGuess
            else:
                print(
                    "Warning: first guess must be float64 ou complex128. It is not used."
                )

        if self.Verbose:
            if UseNoiseRea:
                print(
                    "Sparse rec. with noise rea: ",
                    nx,
                    ny,
                    ", Niter = ",
                    niter,
                    ", Nsigma = ",
                    Nsigma,
                    ", Nrea = ",
                    Nrea,
                    "Inpaint = ",
                    Inpaint,
                    "Hard = ",
                    hard,
                )

        # Initialisation for  inpainting
        WT_Inpaint = 0
        DCT_Inpaint = 0
        if Inpaint:
            DCT_Inpaint = 1  # we will use on WT inpainting,
            # DCT inpating is availabble if we
            # replace WT_Inpaint = 1  by DCT_Inpaint = 1
        if DCT_Inpaint:
            lmin = 0
            lmax = self.get_lmax_dct_inpaint(gamma1, gamma2)
        if WT_Inpaint:
            lmin = 0
            resi1, resi2 = self.get_resi(xg, gamma1, gamma2, WeightResi)
            self.WT.transform(resi1)
            lmax = np.max(np.abs(self.WT.coef[0 : self.WT.ns - 1, :, :]))

        # Main iteration
        Verbose = self.Verbose
        rec = np.zeros((nx, ny))
        reci = np.zeros((nx, ny))
        for n in range(niter):
            xn = np.copy(xg)
            # print("XG=", xg.shape)
            resi1, resi2 = self.get_resi(xg, gamma1, gamma2, WeightResi)
            xg += resi1 + 1j * resi2

            self.WT.transform(xg.real)
            self.WT.threshold(
                SigmaNoise=self.WT_Sigma,
                Nsigma=Nsig,
                ThresCoarse=False,
                hard=hard,
                FirstDetectScale=FirstDetectScale,
                Verbose=False,
            )
            rec[:, :] = self.WT.recons()

            if WT_Inpaint:
                lval = lmin + (lmax - lmin) * (1 - erf(2.8 * n / niter))  # exp decay
                self.WT.threshold(
                    SigmaNoise=1.0, Nsigma=lval, hard=False, KillCoarse=False
                )
                inp = self.WT.recons()
                rec[:, :] = (1 - mask) * inp + mask * xg.real
            xg[:, :] = rec + Iz
            # xg = self.step_wt_recons(t)

            if DCT_Inpaint:
                rec[:, :] = self.step_dct_inpaint(xg, xn, mask, n, niter, lmin, lmax)
                rec[:, :] = (1 - mask) * rec + mask * xg.real
                xg[:, :] = rec + Iz

            if Verbose:
                if ktr is not None:
                    print(
                        "   Sparse rec Iter: ",
                        n + 1,
                        ", Err = %5.4f"
                        % (
                            LA.norm((xg.real - ktr) * mask)
                            / LA.norm(ktr * mask)
                            * 100.0
                        ),
                        ", Resi ke (x100) =  %5.4f" % (LA.norm(resi1 * mask) * 100.0),
                        ", Resi kb (x100) =  %5.4f" % (LA.norm(resi2 * mask) * 100.0),
                    )
                else:
                    print(
                        "   Sparse rec Iter: ",
                        n + 1,
                        ", Sol =  %5.4f" % LA.norm(xg.real),
                    )

        # debiaising sterp
        if self.niter_debias > 0:
            xn = np.copy(xg)
            resi1, resi2 = self.get_resi(xg, gamma1, gamma2, WeightResi)
            xn += resi1 + 1j * resi2
            self.WT.transform(xn.real)
            # determine the set of active coefficients
            if UseNoiseRea:
                self.WT_ActiveCoef = self.get_active_wt_coef(
                    InshearData,
                    WT_Sigma=self.WT_Sigma,
                    UseRea=UseNoiseRea,
                    SigmaNoise=1.0,
                    Nsigma=Nsigma,
                    ComputeWTCoef=False,
                )
            else:
                self.WT_ActiveCoef = self.get_active_wt_coef(
                    InshearData,
                    UseRea=False,
                    SigmaNoise=tau,
                    Nsigma=Nsigma,
                    ComputeWTCoef=False,
                )
            for n in range(self.niter_debias):
                resi1, resi2 = self.get_resi(xg, gamma1, gamma2, mask)
                xg += resi1 + 1j * resi2
                self.WT.coef *= self.WT_ActiveCoef
                rec[:, :] = self.WT.recons()
                xg[:, :] = rec + Iz
        return xg.real, xg.imag

    def sparse_recons_covmat(
        self,
        gamma1,
        gamma2,
        NcvIn,
        niter=None,
        Nsigma=None,
        ThresCoarse=False,
        Inpaint=False,
        ktr=None,
        FirstDetectScale=DEF_FirstDetectScale,
        Bmode=True,
        FirstGuess=None,
    ):
        """
        This routine should not be used anymore.
        It is equivalent to self.sparse_recons using the parameter UseNoiseRea=False
        """

        if niter is None:
            niter = self.DEF_niter
        if Nsigma is None:
            Nsigma = self.DEF_Nsigma

        RMS_ShearMap = np.sqrt(NcvIn / 2.0)
        (nx, ny) = gamma1.shape
        self.nx = nx
        self.ny = ny

        xg = np.zeros((nx, ny)) + 1j * np.zeros((nx, ny))
        Iz = 1j * np.zeros((nx, ny))
        if FirstGuess is not None:
            if FirstGuess.dtype == "float64":
                xg = FirstGuess + Iz
            elif FirstGuess.dtype == "complex128":
                xg = FirstGuess
            else:
                print(
                    "Warning: first guess must be float64 ou complex128. It is not used."
                )

        index = np.where(NcvIn < 1e2)
        mask = np.zeros((nx, ny))
        mask[index] = 1

        # find the minimum noise variance
        tau = np.min(RMS_ShearMap)
        SigmaNoise = tau
        # compute signal coefficient
        Esn = tau / RMS_ShearMap

        if self.Verbose:
            print(
                "Sparse (l_0) Rec. with covmat: ",
                nx,
                ny,
                ", Niter = ",
                niter,
                ", Nsigma = ",
                Nsigma,
                "Inpaint = ",
                Inpaint,
            )

        Verbose = self.Verbose
        if Inpaint:
            lmin = 0
            lmax = self.get_lmax_dct_inpaint(gamma1, gamma2)
        rec = np.zeros((nx, ny))
        reci = np.zeros((nx, ny))
        for n in range(niter):
            xn = np.copy(xg)
            resi1, resi2 = self.get_resi(xg, gamma1, gamma2, Esn)
            xg += resi1 + 1j * resi2
            self.WT.transform(resi1)
            self.WT.threshold(
                SigmaNoise=SigmaNoise,
                Nsigma=Nsigma,
                ThresCoarse=ThresCoarse,
                hard=True,
                FirstDetectScale=FirstDetectScale,
                Verbose=False,
            )
            rec[:, :] = self.WT.recons()
            if Bmode:
                self.WT.transform(xg.imag)
                self.WT.threshold(
                    SigmaNoise=SigmaNoise,
                    Nsigma=Nsigma,
                    ThresCoarse=ThresCoarse,
                    hard=True,
                    FirstDetectScale=FirstDetectScale,
                    Verbose=False,
                )
                reci[:, :] = self.WT.recons()
            else:
                reci[:, :] = 0
            xg[:, :] = rec + 1j * reci

            if Inpaint:
                xg[:, :] = self.step_dct_inpaint(
                    xg, xn, mask, n, niter, lmin, lmax, InpaintAlsoImag=Bmode
                )

            if Verbose:
                ind = np.where(mask == 1)
                if ktr is not None:
                    print(
                        "   Sparse rec Iter: ",
                        n + 1,
                        ", Err = %5.4f"
                        % (
                            np.std((xg.real[ind] - ktr[ind])) / np.std(ktr[ind]) * 100.0
                        ),
                        ", Resi ke (x100) =  %5.4f" % (np.std(resi1[ind]) * 100.0),
                        ", Resi kb (x100) =  %5.4f" % (np.std(resi2[ind]) * 100.0),
                    )
                else:
                    print(
                        "   Sparse rec Iter: ",
                        n + 1,
                        ", Resi ke =  %5.4f" % (np.std(resi1[ind] / tau)),
                        ", Resi kb = %5.4f" % (np.std(resi2[ind]) / tau),
                    )

        return xg.real, xg.imag

    def rea_sparse_wiener_filtering(
        self,
        InshearData,
        PowSpecSignal,
        WT_Support=None,
        niter=None,
        WT_Sigma=None,
        Nsigma=DEF_Nsigma,
        Inpaint=False,
        FirstDetectScale=DEF_FirstDetectScale,
        Nrea=DEF_Nrea,
        ktr=None,
    ):
        """
        This routine should not be used. It is a test routine
        """

        if niter is None:
            niter = self.DEF_niter
        if Nsigma is None:
            Nsigma = self.DEF_Nsigma

        if WT_Support is None:
            WT_Support = self.get_active_wt_coef(
                InshearData,
                UseRea=True,
                SigmaNoise=1.0,
                Nsigma=Nsigma,
                Nrea=Nrea,
                WT_Sigma=WT_Sigma,
            )
        if FirstDetectScale > 0:
            WT_Support[0:FirstDetectScale, :, :] = 0
        WT_Support[0:FirstDetectScale, :, :] = 0
        WT_Support[self.WT.ns - 1, :, :] = 0
        gamma1 = InshearData.g1
        gamma2 = InshearData.g2
        (nx, ny) = gamma1.shape
        mask = InshearData.mask
        self.nx = nx
        self.ny = ny
        xg = np.zeros((nx, ny)) + 1j * np.zeros((nx, ny))
        mask = InshearData.mask
        if FirstDetectScale > 0:
            WT_Support[0:FirstDetectScale] = 0

        Ncv = InshearData.Ncov / 2.0
        xg = np.zeros((nx, ny)) + 1j * np.zeros((nx, ny))
        xs = np.zeros((nx, ny)) + 1j * np.zeros((nx, ny))
        # find the minimum noise variance
        tau = np.min(Ncv)

        # set the step size
        # eta = 1.83 * tau
        eta = tau

        # compute signal coefficient
        Esn = eta / Ncv

        # calculate the wiener filter coefficients
        Px_map = get_ima_spectrum_map(PowSpecSignal, nx, ny)

        Wfc = Px_map / (Px_map + eta)
        if Inpaint:
            lmin = 0
            lmax = self.get_lmax_dct_inpaint(gamma1, gamma2)

        xs = xg
        nw = 1
        for n in range(niter):
            xn = np.copy(xg)
            resi1, resi2 = self.get_resi(xg, gamma1, gamma2, mask)
            # sparse component
            xs = resi1 + 1j * resi2  # xg + H^T(eta / Sn * (y- H * xg))
            self.WT.transform(xs.real)
            self.WT.coef *= WT_Support
            xs.real = self.WT.recons()
            xs.imag[:, :] = 0

            # Winer component
            # calculate the residual
            xg = xg - xs
            for i in 0, range(nw):
                t1, t2 = self.H_operator_eb2g(
                    xg.real + xs.real, xg.imag + xs.imag
                )  # H * xg
                t1, t2 = self.H_adjoint_g2eb(
                    Esn * (gamma1 - t1), Esn * (gamma2 - t2)
                )  # H^T(eta / Sn * (y- H * xg))
                t = xg + (t1 + 1j * t2)  # xg + H^T(eta / Sn * (y- H * xg))
                xg = self.mult_wiener(t, Wfc)  # wiener filtering in fourier space

            xg = xg + xs
            if Inpaint:
                xg = self.step_dct_inpaint(xg, xn, mask, n, niter, lmin, lmax)

            if self.Verbose:
                ind = np.where(mask == 1)
                if ktr is not None:
                    print(
                        "   Sparse rec Iter: ",
                        n + 1,
                        ", Err = %5.4f"
                        % (
                            np.std((xg.real[ind] - ktr[ind])) / np.std(ktr[ind]) * 100.0
                        ),
                        ", Resi ke (x100) =  %5.4f" % (np.std(resi1[ind]) * 100.0),
                        ", Resi kb (x100) =  %5.4f" % (np.std(resi2[ind]) * 100.0),
                    )
                else:
                    print(
                        "   Sparse rec Iter: ",
                        n + 1,
                        ", Resi ke =  %5.4f" % (np.std(resi1[ind] / tau)),
                        ", Resi kb = %5.4f" % (np.std(resi2[ind]) / tau),
                    )

        return xg.real, xs.real


############ END CLASS #######################


# from lenspack.utils
def bin2d(x, y, npix=10, v=None, w=None, extent=None, verbose=False):
    """Bin samples of a spatially varying quantity according to position.

    The (weighted) average is taken of values falling into the same bin. This
    function is relatively general, but it is mainly used within this package
    to produce maps of the two components of shear from a galaxy catalog.

    Parameters
    ----------
    x, y : array_like
        1D position arrays.
    npix : int or list or tuple as (nx, ny), optional
        Number of bins in the `x` and `y` directions. If an int N is given,
        use (N, N). Binning defaults to (10, 10) if not provided.
    v : array_like, optional
        Values at positions (`x`, `y`). This can be given as many arrays
        (v1, v2, ...) of len(`x`) to bin simultaneously. If None, the bin
        count in each pixel is returned.
    w : array_like, optional
        Weights for `v` during averaging. If provided, the same weights are
        applied to each input `v`.
    extent : array_like, optional
        Boundaries of the resulting grid, given as (xmin, xmax, ymin, ymax).
        If None, bin edges are set as the min/max coordinate values of the
        input position arrays.
    verbose : boolean, optional
        If True, print details of the binning.

    Returns
    -------
    ndarray or tuple of ndarray
        2D numpy arrays of values `v` binned into pixels. The number of
        outputs matches the number of input `v` arrays.

    Examples
    --------
    >>> # 100 values at random positions within the ranges -0.5 < x, y < 0.5
    >>> # and binned within -1 < x, y < 1 to a (5, 5) grid.
    >>> x = np.random.random(100) - 0.5
    >>> y = np.random.random(100) - 0.5
    >>> v = np.random.randn(100) * 5
    >>> bin2d(x, y, v=v, npix=5, extent=(-1, 1, -1, 1))
    array([[ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ],
           [ 0.        ,  4.43560619, -2.33308373,  0.48447844,  0.        ],
           [ 0.        ,  1.94903524, -0.29253335,  1.3694618 ,  0.        ],
           [ 0.        , -1.0202718 ,  0.37112266, -1.43062585,  0.        ],
           [ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ]])

    """
    # Regroup extent if necessary
    if extent is not None:
        assert len(extent) == 4
        extent = [extent[:2], extent[2:]]

    if v is None:
        # Return the simple bin count map
        bincount, xbins, ybins = np.histogram2d(x, y, bins=npix, range=extent)
        result = bincount.T
    else:
        # Prepare values to bin
        v = np.atleast_1d(v)
        if len(v.shape) == 1:
            v = v.reshape(1, len(v))

        # Prepare weights
        if w is not None:
            w = np.atleast_1d(w)
            has_weights = True
        else:
            w = np.ones_like(x)
            has_weights = False

        # Compute weighted bin count map
        wmap, xbins, ybins = np.histogram2d(x, y, bins=npix, range=extent, weights=w)
        # Handle division by zero (i.e., empty pixels)
        wmap[wmap == 0] = np.inf
        # Compute mean values per pixel
        result = tuple(
            (
                np.histogram2d(x, y, bins=npix, range=extent, weights=(vv * w))[0]
                / wmap
            ).T
            for vv in v
        )

        # Clean up
        if len(result) == 1:
            result = result[0]

    if verbose:
        if v is not None:
            print(
                "Binning {} array{} with{} weights.".format(
                    len(v), ["", "s"][(len(v) > 1)], ["out", ""][has_weights]
                )
            )
        else:
            print("Returning bin count map.")
        print("npix : {}".format(npix))
        print("extent : {}".format([xbins[0], xbins[-1], ybins[0], ybins[-1]]))
        print("(dx, dy) : ({}, {})".format(xbins[1] - xbins[0], ybins[1] - ybins[0]))

    return result


# from lenspack.geometry.projections.gnom
def radec2xy(ra0, dec0, ra, dec):
    """Project spherical sky coordinates to a tangent plane.

    Parameters
    ----------
    ra0 : float
        Right ascension of the projection origin.
    dec0 : float
        Declination of the projection origin.
    ra : float or array_like
        Right ascension of point(s) to project.
    dec : float or array_like
        Declination of point(s) to project.

    Notes
    -----
    All input units are assumed to be degrees.

    Returns
    -------
    x, y : tuple of floats or numpy arrays
        Projected coordinate(s) in the tangent plane relative to (0, 0), i.e.
        the origin in the projected space.

    Raises
    ------
    Exception
        For input arrays of different sizes.

    Examples
    --------
    ...

    """
    # Standardize inputs
    ra = np.atleast_1d(ra)
    dec = np.atleast_1d(dec)

    if len(ra) != len(dec):
        raise Exception("Input ra and dec must have the same length.")

    # Convert input coordinates to radians
    alpha0 = np.deg2rad(ra0)
    delta0 = np.deg2rad(dec0)
    alpha = np.deg2rad(ra)
    delta = np.deg2rad(dec)

    # Project points
    denom = np.cos(delta0) * np.cos(delta) * np.cos(alpha - alpha0) + np.sin(
        delta0
    ) * np.sin(delta)
    x = np.cos(delta) * np.sin(alpha - alpha0) / denom
    y = (
        np.cos(delta0) * np.sin(delta)
        - np.sin(delta0) * np.cos(delta) * np.cos(alpha - alpha0)
    ) / denom

    # Potentially remove unnecessary array layers
    if len(x) == 1:
        x, y = x[0], y[0]

    return x, y
