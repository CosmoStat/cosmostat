#! /usr/bin/env Python
"""
Created on Nov 7, 2023

@authors: Andreas Tersenov & Jean-Luc Starck
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
from pycs.astro.wl.mass_mapping import *


def get_wt_noiselevel(W, NoiseSigmaMap, Mask=None):
    """
    Compute the noise standard deviation relative to each wavelet coefficient, when
    the input noise is Gaussian, but not necessary stationary.
    NoiseStdMap is the input noise standard deviation map.

    Parameters
    ----------
    W : Class obtained with the command WT = starlet2d()
        The class has to be initialized.
        Ex: WT = starlet2d(gen2=False,l2norm=False, verb=False)
            WT.init_starlet(nx, ny, nscale=5)
       Starlet transform call (see starlet.py in CosmoStat package).
    NoiseSigmaMap : 2D array.
        Noise standard deviation map.
        It must have the same size as the size given in the starlet class W.
    Mask: 2D array
        mask on the data. Mask[i,j] = 1 if the pixel is observed, 0 otherwise
    Returns
    -------
    3D array
        x = get_wt_noiselevel(W, NoiseStdMap):
        x[i,:,:] is the noise standard deviation map at the scale i.
    """
    if Mask is None:
        NoiseStdMap = np.copy(NoiseSigmaMap)
    else:
        # we will consider in the non observed area (mask == 0) that
        # the noise is equivalent to the maximum noise level in the map.
        # This gives an over-estimation of the noise at these locations after
        # inpaiting.
        NoiseStdMap = np.copy(NoiseSigmaMap)
        ind = np.where(Mask == 1)
        MaxNoise = np.max(NoiseStdMap[ind])
        ind = np.where(Mask == 0)
        NoiseStdMap[ind] = MaxNoise

    nx = W.nx
    ny = W.ny
    im = np.zeros((nx, ny))
    im = im.astype("float64")
    im[int(nx / 2), int(ny / 2)] = np.float64(1.0)
    wt = W.transform(im)
    VarianceCoeff = W.coef * 0.0
    VarianceMap = NoiseStdMap**2
    for i in np.arange(W.ns):
        coef = W.coef[i, :, :]
        coef = coef**2
        VarianceCoeff[i, :, :] = conv(VarianceMap, coef)
        VarianceCoeff[VarianceCoeff < 0] = 0
        # info(VarianceCoeff[i, :, :])
    return np.sqrt(VarianceCoeff)


def get_snr_noiselevel(W, Map, Wnoise, KeepSign=False):
    """
    Take the wavelet transform of an image and compute the Signal-to-Noise Ratio (SNR)
    for each wavelet coefficient.

    Parameters
    ----------
    W : Starlet Call
        Wavelet transform class.
    Map : 2D array
        Image to analyse.
    Wnoise : 3D array
        Wnoise[j,x,y] is the noise level for a wavelet coefficient a scale j and position (x,y)
    KeepSign : Boolean, optional
        If KeepSign is True, the sign of the wavelet coefficients are kept, and the
        SNR out cube contains both positive and negative values. The default is False.

    Returns
    -------
    WSnr : 3D array
        WSnr[j,x,y] is the SNR for a wavelet coefficient a scale j and position (x,y).

    """
    W.transform(Map)
    WSnr = np.copy(Wnoise)
    WSnr[:, :, :] = 0.0
    ind = np.where(Wnoise != 0)
    WSnr[ind] = W.coef[ind] / Wnoise[ind]
    if KeepSign is False:
        WSnr = np.abs(WSnr)
    return WSnr


def get_peaks(image, threshold=None, ordered=True, mask=None, include_border=False):
    """Identify peaks in an image (2D array) above a given threshold.

    A peak, or local maximum, is defined as a pixel of larger value than its
    eight neighbors. A mask may be provided to exclude certain regions from
    the search. The border is excluded by default.

    Parameters
    ----------
    image : array_like
        Two-dimensional input image.
    threshold : float, optional
        Minimum pixel amplitude to be considered as a peak. If not provided,
        the default value is set to the minimum of `image`.
    ordered : bool, optional
        If True, return peaks in decreasing order according to height.
    mask : array_like (same shape as `image`), optional
        Boolean array identifying which pixels of `image` to consider/exclude
        in finding peaks. A numerical array will be converted to binary, where
        only zero values are considered masked.
    include_border : bool, optional
        If True, include peaks found on the border of the image. Default is
        False.

    Returns
    -------
    X, Y, heights : tuple of 1D numpy arrays
        Pixel indices of peak positions and their associated heights.

    Notes
    -----
    From LensPack, added here to avoid circular dependencies

    Examples
    --------
    ...

    """
    image = np.atleast_2d(image)

    # Deal with the mask first
    if mask is not None:
        mask = np.atleast_2d(mask)
        if mask.shape != image.shape:
            print("Warning: mask not compatible with image -> ignoring.")
            mask = np.ones(image.shape)
        else:
            # Make sure mask is binary, i.e. turn nonzero values into ones
            mask = mask.astype(bool).astype(float)
    else:
        mask = np.ones(image.shape)

    # Add 1 pixel padding if including border peaks
    if include_border:
        image = np.pad(image, pad_width=1, mode="constant", constant_values=image.min())
        mask = np.pad(mask, pad_width=1, mode="constant", constant_values=1)

    # Determine threshold level
    if threshold is None:
        # threshold = image[mask.astype('bool')].min()
        threshold = image.min()
    else:
        threshold = max(threshold, image.min())

    # Shift everything to be positive to properly handle negative peaks
    offset = image.min()
    threshold = threshold - offset
    image = image - offset

    # Extract the center map
    map0 = image[1:-1, 1:-1]

    # Extract shifted maps
    map1 = image[0:-2, 0:-2]
    map2 = image[1:-1, 0:-2]
    map3 = image[2:, 0:-2]
    map4 = image[0:-2, 1:-1]
    map5 = image[2:, 1:-1]
    map6 = image[0:-2, 2:]
    map7 = image[1:-1, 2:]
    map8 = image[2:, 2:]

    # Compare center map with shifted maps
    merge = (
        (map0 > map1)
        & (map0 > map2)
        & (map0 > map3)
        & (map0 > map4)
        & (map0 > map5)
        & (map0 > map6)
        & (map0 > map7)
        & (map0 > map8)
    )

    bordered = np.lib.pad(merge, (1, 1), "constant", constant_values=(0, 0))
    peaksmap = image * bordered * mask
    X, Y = np.nonzero(peaksmap > threshold)

    # Extract peak heights
    heights = image[X, Y] + offset

    # Compensate for border padding
    if include_border:
        X = X - 1
        Y = Y - 1

    # Sort peaks according to height
    if ordered:
        inds = np.argsort(heights)[::-1]
        return X[inds], Y[inds], heights[inds]

    return X, Y, heights


# Class to compute multiscale peaks and wavelet l1 norm for an image
class HOS_starlet_l1norm_peaks:
    NBins = 0  # int: Number of bins to use
    TabBins = 0  # array indicated the bins
    WT = 0  # WT transform class
    Ns = 0  # int, number of scales
    DataCoef = 0  # 3D array, containing the wavelet transform of the data
    # 3D array, containing the noise standard deviation for each wavelet
    # coefficient
    NoiseCoef = 0
    WSNR = (
        0  # 3D array, containing the signal to noise ratio for each wavelet coefficient
    )
    PIX_ARCMIN = 1  # pixel resolution
    scales_arcmin = 0  # 1D array, resolution in arcmin of each scale
    TabBinsCenter = 0  # 1D array, center of each bin
    Peaks_PosX = 0  # list of X peak positions at each scale
    Peaks_PosY = 0  # list of Y peak positions at each scale
    Peaks_Height = 0  # list of  peak amplitudes at each scale
    Peaks_Count = 0  # list of  peak counts at each scale
    Mono_Peaks_PosX = 0  # list of X mono-scale peak positions
    Mono_Peaks_PosY = 0  # list of Y mono-scale peak positions
    Mono_Peaks_Height = 0  # list of  mono-scale peak amplitudes
    Mono_Peaks_Count = 0  # list of  mono-scale peak counts
    l1bins = 0  # list of l1 bins at each scale
    l1norm = 0  # list of l1 norm at each scale

    def __init__(self, WTrans, PixResolArcMin=1):  # __init__ is the constructor
        self.WT = WTrans  # Starlet wavelet Class defined in starlet.py
        self.Ns = WTrans.ns  # number of scales
        self.PIX_ARCMIN = PixResolArcMin  # resolution
        self.scales_arcmin = [2 ** (i + 1) * PixResolArcMin for i in range(self.Ns)]

    def set_data(self, Map, SigmaMap, Mask=None):
        """
        Calculate the wavelet transform of the data, and estimation the noise standard deviation
        for every wavelet coefficient.
        If some pixels don't have a value, a mask has to be given.
        Parameters
        ----------
        Map : 2D array
            IMage to analyze.
        SigmaMap : 2D array
            Noise standard deviation map.
        Mask : 2D array, optional
            Mask indicating where we have observations. Mask[x,y]=1 means we have an
            observation at location (x,y), Mask[x,y]=0 otherwise. The default is None.

        Returns
        -------
        None.

        """
        self.NoiseCoef = get_wt_noiselevel(self.WT, SigmaMap, Mask=Mask)
        self.WSNR = get_snr_noiselevel(self.WT, Map, self.NoiseCoef, KeepSign=True)
        self.DataCoef = self.WT.coef

    def set_bins(self, Min=-2, Max=6, nbins=31):
        """
        Set the bins for the histogram and the l1-norm.

        Parameters
        ----------
        Min : int, optional
            Minimm value for the histogram. The default is -2.
        Max : int, optional
            Maximum value for the histogram. The default is 6.
        nbins : int, optional
            Number of bins. The default is 31.

        Returns
        -------
        None.

        """
        self.NBins = nbins
        self.TabBins = np.linspace(Min, Max, self.NBins)
        self.TabBinsCenter = 0.5 * (self.TabBins[:-1] + self.TabBins[1:])

    def tvscalepeaks(self, ScaleNumber, maxcoef=None, lut="inferno", OnlyPeaks=None):
        if maxcoef is None:
            maxcoef = 2.0
        i = ScaleNumber
        Scale = self.DataCoef[i, :, :]

        X = self.Peaks_PosX[i]
        Y = self.Peaks_PosY[i]
        if OnlyPeaks is not None:
            Scale[:, :] = 0
            Scale[X, Y] = 1
        else:
            Max = Scale.max()
            Scale[X, Y] = maxcoef * Max
        tvimap(Scale, title="Scale" + str(i + 1), lut=lut)

    def tvpeaks(self, maxcoef=2, lut="inferno", OnlyPeaks=None):
        """
        Plot the scales, chaning the values at the position of the peaks
        to see them. The values at pixel positions are set to twice the maximum.

        Parameters
        ----------
        maxcoef : float, optional
            Multiplicative factor applied on the maximum value. The default is 2.
        lut : string, optional
            selected LUT. The default is 'inferno'.
        OnlyPeaks : Boolean, optional
            If set, only peaks are shown, not the wavelet images. The default is None.

        Returns
        -------
        None.

        """
        for i in np.arange(self.WT.ns):
            self.tvscalepeaks(i, maxcoef=maxcoef, lut=lut, OnlyPeaks=OnlyPeaks)

    def get_mono_scale_peaks(self, image, sigma_noise, smoothing_sigma=2, mask=None):
        """
        Calculate mono-scale peak counts with Gaussian smoothing.

        Parameters
        ----------
        image : 2D array
            Input image.
        sigma_noise : float
            Standard deviation of the noise.
        smoothing_sigma : float, optional
            Standard deviation for Gaussian smoothing. Default is 2.
        mask : 2D array, optional
            Mask indicating where observations are present.

        Returns
        -------
        None.

        """
        # create a Gaussian kernel
        im = np.zeros_like(image)
        im[im.shape[0] // 2, im.shape[1] // 2] = 1.0
        gaussian_kernel = ndimage.gaussian_filter(im, smoothing_sigma)

        # calculate the noise level for the smoothed pixels
        variance_map = sigma_noise**2
        square_gaussian_kernel = gaussian_kernel**2
        smoothed_noise_sigma = np.sqrt(conv(variance_map, square_gaussian_kernel))

        # calculate the SNR image
        image_smoothed = ndimage.gaussian_filter(image, smoothing_sigma)
        snr_image = np.zeros_like(image_smoothed)
        ind = np.where(smoothed_noise_sigma != 0)
        snr_image[ind] = image_smoothed[ind] / smoothed_noise_sigma[ind]

        # Get mono-scale peaks
        X, Y, heights = get_peaks(snr_image, mask=mask)
        self.Mono_Peaks_PosX = X
        self.Mono_Peaks_PosY = Y
        self.Mono_Peaks_Height = heights

        # Calculate histogram of peak counts
        counts, bin_edges = np.histogram(heights, bins=self.TabBins)
        self.Mono_Peaks_Count = counts

    def get_wtpeaks(self, Mask=None, verbose=False):
        """
        Calculate the histogram of of peak counts at all scales

        Parameters
        ----------
        Mask : 2D array, optional
            DESCRIPTION. The default is None.
            Caculate the histogram of of peak counts  only in the mask area.
            The default is None, the whole image is used at every scale.
        verbose : bool, optional
            DESCRIPTION. If True, print the minimum and maximum values at each scale. The default is False.
        Returns
        -------
        None.

        """
        self.Peaks_Count = []
        self.Peaks_PosX = []
        self.Peaks_PosY = []
        self.Peaks_Height = []
        for i in np.arange(self.WT.ns):
            Scale = self.WSNR[i, :, :]
            X, Y, heights = get_peaks(
                Scale, threshold=None, ordered=True, mask=Mask, include_border=False
            )
            self.Peaks_PosX.append(X)
            self.Peaks_PosY.append(Y)
            self.Peaks_Height.append(heights)
            counts, bin_edges = np.histogram(heights, bins=self.TabBins)
            self.Peaks_Count.append(counts)

            if verbose:
                Min = Scale.min()
                Max = Scale.max()
                print("Scale ", i + 1, ": Min = ", Min, ", Max = ", Max)

    def plot_mono_peaks_histogram(self):
        """
        Plot histogram of mono-scale peak counts.

        Returns
        -------
        None.
        """
        plt.figure()
        plt.plot(self.TabBinsCenter, self.Mono_Peaks_Count)
        plt.title("Mono-Scale Peak Counts Histogram")
        plt.xlabel("SNR")
        plt.ylabel("Counts")
        plt.yscale("log")
        plt.grid()
        plt.show()

    def plot_peaks_histo(
        self,
        Scale=None,
        title="Starlet Peak Counts Histogram",
        xlabel="SNR",
        ylabel="Peak Counts",
        log_scale=False,
    ):
        """
        PLot either the histogram of peak counts for a given scale or for all scales.
        Parameters:
        ----------
            title (str): Title of the plot.
            xlabel (str): Label for the x-axis.
            ylabel (str): Label for the y-axis.
            log_scale (bool, optional): Whether to use a logarithmic scale for the y-axis. Default is False.
        """
        plt.figure()
        if Scale is not None:
            pc = self.Peaks_Count[Scale]
            plt.plot(self.TabBinsCenter, pc)
        else:
            if isinstance(self.Peaks_Count, list):  # Multiscale case
                for scale, pc in enumerate(self.Peaks_Count):
                    plt.plot(self.TabBinsCenter, pc, label=f"Scale {scale + 1}")

        plt.legend()
        if log_scale:
            plt.yscale("log")
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.grid(True)
        plt.title(title)
        plt.show()

    def plot_l1norm(
        self,
        Scale=None,
        xlim=None,
        title="Starlet l1-norm",
        xlabel="SNR",
        ylabel="l1-norm",
        log_scale=False,
    ):
        """
        PLot either the l1 norm per bin for a given scale or for all scales.
        Parameters
        ----------
        Scale : int, optional
            Plot the l1 norm per bin at the specified scale. Scale=0 for the first scale.
            The default is None, all scales are plotted.
        xlim : int, optional
            PLot only values for x > xlim. The default is None.
        title : string, optional
            Title of the plot. The default is 'Starlet l1-norm'.
        xlabel : string, optional
            x-label. The default is 'SNR'.
        ylabel : string, optional
            y label. The default is 'l1-norm'.
        log_scale : Boolean, optional
            Use log scale for the y axis . The default is False.

        Returns
        -------
        None.

        """

        plt.figure()

        if Scale is not None:
            l1 = self.l1norm[Scale]
            plt.plot(self.l1bins[Scale], l1, label=f"Scale {Scale}")
        else:
            for scale, l1 in enumerate(self.l1norm):
                plt.plot(self.l1bins[scale], l1, label=f"Scale {scale}")
        plt.legend()
        plt.xticks()
        plt.yticks()
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.grid(True)
        plt.title(title)
        if xlim:
            plt.xlim(xlim)
        if log_scale:
            plt.yscale("log")
        plt.show()

    def get_wtl1(self, nbins=None, Mask=None, min_snr=None, max_snr=None):
        """
        Calculate the wavelet l1 norm per bin.

        Parameters
        ----------
        nbins : int, optional
            Number of bins. The default is NBins, defined using the set_bins function.
        Mask : 2D array, optional
            Caculate the l1 norm only in the mask area.
            The default is None, the whole image is used at every scale.
        min_snr : float, optional
            Minimum SNR value for calculating l1 norms. The default is None, which will use the minimum SNR found in the data.
        max_snr : float, optional
            Maximum SNR value for calculating l1 norms. The default is None, which will use the maximum SNR found in the data.


        Returns
        -------
        None.

        """
        if nbins is None:
            nbins = self.NBins

        l1_coll = []
        bins_coll = []

        for i in np.arange(self.WT.ns):
            ScaleSNR = self.WSNR[i, :, :]
            if Mask is not None:
                ind = np.where(Mask == 0)
                ScaleSNR = ScaleSNR[ind]

            min_snr_value = min_snr if min_snr is not None else np.min(ScaleSNR)
            max_snr_value = max_snr if max_snr is not None else np.max(ScaleSNR)

            thresholds_snr = np.linspace(min_snr_value, max_snr_value, nbins + 1)
            bins_snr = 0.5 * (thresholds_snr[:-1] + thresholds_snr[1:])
            digitized = np.digitize(ScaleSNR, thresholds_snr)
            bin_l1_norm = [
                np.sum(np.abs(ScaleSNR[digitized == j]))
                for j in range(1, len(thresholds_snr))
            ]
            bins_coll.append(bins_snr)
            l1_coll.append(bin_l1_norm)

        self.l1bins = np.array(bins_coll)
        self.l1norm = np.array(l1_coll)


#############  TESTS routine ############


def get_rms_error(Res, TrueSol, Mask, sigma=0):
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


def test_hos_test1(lut="inferno"):
    # For space constraints, we provide this experiment starting from a pixelized
    # shear map, and not from the catalog.

    # Read shear data and covariance matrix
    # DIR='/Users/starck/git/cosmostat/cosmostat/examples/mcalens_paper/'
    DIR = "/Users/atersenov/Software/cosmostat_hos/cosmostat/examples/mcalens_paper/"
    g1 = readfits(DIR + "exp_wiener_miceDSV_g1.fits")
    g2 = readfits(DIR + "exp_wiener_miceDSV_g2.fits")
    Ncov = readfits(DIR + "exp_wiener_miceDSV_covmat.fits")

    # Read the true convergence map and its theoretical power spectrum
    ktr = readfits(DIR + "exp_wiener_miceDSV_true_convergence_map.fits")
    ps1d = readfits(DIR + "exp_wiener_miceDSV_signal_powspec.fits")

    # Mass mapping class initialization
    (nx, ny) = Ncov.shape
    M = massmap2d(name="mass")
    M.init_massmap(nx, ny)
    M.DEF_niter = 200
    Inpaint = True
    M.niter_debias = 30
    M.Verbose = True
    g1t, g2t = M.k2g(ktr)

    # Read the noise power spectrum. It has been derived from
    # from the covariance matrix using 1000 realization with the command:
    #    Pn = M.get_noise_powspec(Ncov ,mask=mask,nsimu=1000, inpaint=True)
    pn = readfits(DIR + "exp_wiener_miceDSV_noise_powspec.fits")

    # derive the mask from the covariance matrix
    index = np.where(Ncov < 1e2)
    mask = np.zeros((nx, ny))
    mask[index] = 1
    ind = np.where(mask != 1)
    mask[ind] = 0

    # Make simu. Use a convergence image of 512x512,
    # get gamma (512x512)
    # extract the center of the gamma image ==> get 256x256 shear image
    # reconstruct the center of kappa, and as it is reconstructed from
    # the truncated gamma map, we can see the effect of the border.
    # We use mice noise covariance to have a realistic sky coverage and noise variance
    # We have a parameter to improve the SNR to better investigate the noise
    # behavior
    AnotherSIMU = False
    if AnotherSIMU is True:
        # DIR='/Users/starck/tex/PPT/Y23/EuclidTutorialSchool2023/'
        DIR = "/Users/atersenov/Software/test_cosmostat/data/"
        M512 = massmap2d(name="mass")
        M512.init_massmap(nx, ny)
        ktr512 = readfits(DIR + "WLconv_z0.50_3316r.fits")
        ktr = ktr512[128 : 128 + 256, 128 : 128 + 256]
        g1t512, g2t512 = M512.k2g(ktr512)
        g1t = g1t512[128 : 128 + 256, 128 : 128 + 256]
        g2t = g2t512[128 : 128 + 256, 128 : 128 + 256]
        krec = M.g2k(g1t, g2t)
        # DIR='/Users/starck/git/cosmostat/cosmostat/examples/mcalens_paper/'
        DIR = (
            "/Users/atersenov/Software/cosmostat_hos/cosmostat/examples/mcalens_paper/"
        )
        Ncov = readfits(DIR + "exp_wiener_miceDSV_covmat.fits")
        ImproveSNR = 10.0
        Ncov[index] = Ncov[index] / ImproveSNR
        n1 = np.random.normal(loc=0.0, scale=np.sqrt(Ncov / 2.0)) * mask
        n2 = np.random.normal(loc=0.0, scale=np.sqrt(Ncov / 2.0)) * mask
        g1 = (g1t + n1) * mask
        g2 = (g2t + n2) * mask

    # Shear data class initialization
    d = shear_data()
    d.g1 = g1
    d.g2 = g2
    d.Ncov = Ncov
    d.mask = mask
    d.ps1d = ps1d

    # make Fig. 8 of MCAlens paper
    lut = "inferno"
    tvilut(
        mask,
        title="DES-MICE Mask",
        lut=lut,
        filename="fig_mice_mask.png",
        vmin=0,
        vmax=1,
    )
    tvilut(
        smooth2d(ktr, 1.0) * mask,
        title="Convergence map",
        lut=lut,
        filename="fig_mice_kappa.png",
        vmin=-0.03,
        vmax=0.03,
    )
    tvilut(g1, title="g1", lut=lut, filename="fig_mice_g1.png", vmin=-0.3, vmax=0.3)
    tvilut(g2, title="g2", lut=lut, filename="fig_mice_g2.png", vmin=-0.3, vmax=0.3)

    # Proximal iterative filtering.
    ke_inp_pwiener, kb_winp = M.prox_wiener_filtering(
        g1, g2, ps1d, Ncov, Pn=pn, Inpaint=True
    )  # , Pn=Pn) # ,ktr=InShearData.ktr)
    tvilut(
        ke_inp_pwiener,
        title="Inpainted Wiener",
        lut=lut,
        filename="fig_mice_inp_wiener.png",
        vmin=-0.004,
        vmax=0.004,
    )

    # KS
    ks = M.gamma_to_cf_kappa(g1, g2)
    ksb = ks.imag
    ks = ks.real
    tvilut(ks, title="Kaiser-Squire", lut=lut)

    # KS + Inpainting
    ksi = M.iks(g1, g2, mask)
    ksi = ksi.real
    tvilut(ksi, title="iKS with Inpainting", lut=lut)

    # Prox MSE + Inpainting

    M = massmap2d(name="mass")
    M.init_massmap(nx, ny)
    M.DEF_niter = 200
    Inpaint = True
    M.niter_debias = 30
    M.Verbose = True

    sig = 3
    M.Verbose = True
    ProxMSE, ProxMSE_B = M.prox_mse(g1, g2, Ncov, ktr=ktr, Inpaint=True, sigma=0)
    tvilut(ProxMSE, title="Prox MSE with Inpainting", lut=lut)
    Err_proxmse = get_rms_error(ProxMSE, ktr, mask, sigma=sig)
    print("   Prox MSE Error: ", Err_proxmse)

    # Example of error calculation at a give scale
    sig = 3
    tvilut(smooth2d((ksi - ktr) * mask, sig))
    tvilut(smooth2d((ProxMSE - ktr) * mask, sig))

    tvilut(smooth2d((ks) * mask, sig))
    tvilut(smooth2d((ProxMSE) * mask, sig))

    tvilut(((ks) * mask))
    tvilut(((ProxMSE) * mask))

    sig = 4  # sigma allows us to see the error at a given scale
    print("== ERROR Sigma = ", sig)
    Err_ks = get_rms_error(ks, ktr, mask, sigma=sig)
    print("   Kaiser Squires Error: ", Err_ks)
    Err_ksi = get_rms_error(ksi, ktr, mask, sigma=sig)
    print("   Iterative Kaiser Squires (with inpainting) Error: ", Err_ksi)
    Err_kiw = get_rms_error(ke_inp_pwiener, ktr, mask, sigma=sig)
    print("   Iterative Wiener Error: ", Err_kiw)
    Err_proxmse = get_rms_error(ProxMSE, ktr, mask, sigma=sig)
    print("   Prox MSE Error: ", Err_proxmse)

    # Prox MSE Error:  623.308541569542  sans inpaint
    # Prox MSE Error:  32.79202612245059 sans inpaint sig=5

    # MCAlens
    vm = 0.2
    vmin = -0.1
    M.Verbose = True
    M.DEF_niter = 100
    UseMCAlens = False
    if UseMCAlens is True:
        k1r5, k1i, k2r5, k2i = M.sparse_wiener_filtering(
            d, d.ps1d, Nsigma=5, Inpaint=True, Bmode=True, InpNiter=20
        )
        tvilut(k1r5, title="MCAlens (lambda=5)", fs=5, vmin=vmin, vmax=vm, lut=lut)
        tvilut(
            (k1r5 - k2r5),
            title="MCAlens-Gaussian Component (lambda=5)",
            fs=5,
            vmin=vmin,
            vmax=vm,
            lut=lut,
        )
        tvilut(
            k2r5,
            title="MCAlens-Sparse Component (lambda=5)",
            fs=5,
            vmin=vmin,
            vmax=vm,
            lut=lut,
        )
        Err_mcalens = get_rms_error(k1r5, ktr, mask, sigma=sig)
        print("Iterative MCAlens: ", Err_mcalens)

    # Statistics
    MinSNR = -2
    MaxSNR = 6
    Nb = 31
    nx = g1.shape[0]
    nx = g2.shape[1]
    ns = 5

    WT = starlet2d(gen2=False, l2norm=False, verb=False)
    WT.init_starlet(nx, ny, nscale=ns)
    H = HOS_starlet_l1norm_peaks(WT)

    H.set_bins(Min=MinSNR, Max=MaxSNR, nbins=Nb)
    SigmaMap = np.sqrt(d.Ncov / 2)
    H.set_data(ksi, SigmaMap, Mask=mask)

    info(H.DataCoef, name="DATA")
    info(H.NoiseCoef, name="NoiseCoef")
    info(H.WSNR, name="WSNR")

    H.get_wtpeaks(Mask=mask)
    pc = H.Peaks_Count
    H.plot_peaks_histo(log_scale=True)

    H.tvpeaks(maxcoef=3, lut="rainbow", OnlyPeaks=True)
    H.get_wtl1(Nb * 2)
    H.plot_l1norm()


def test_hos_cfis():
    # Exemple for the CFIS experiment
    # DIR='/Users/starck/tex/PPT/Y23/EuclidTutorialSchool2023/'
    DIR = "/Users/atersenov/Software/test_cosmostat/data/"
    kappa_map = readfits(DIR + "WLconv_z0.50_3316r.fits")
    SHAPE_NOISE = 0.44
    PIX_ARCMIN = 0.4
    N_GAL = 7
    NSCALES = 5
    NBINS = 40
    KAPPA_SNR = np.linspace(-2, 6, 31)

    # Make a noise
    sigma_noise_CFIS = SHAPE_NOISE / (np.sqrt(2 * N_GAL * PIX_ARCMIN**2))
    noise_map_CFIS_z05 = sigma_noise_CFIS * np.random.randn(
        kappa_map.shape[0], kappa_map.shape[1]
    )  # generate noise map
    kappa_map_noisy = kappa_map + noise_map_CFIS_z05  # Add noise to the mass map

    # tvilut(k2r5*mask, title='CFIS kappa',fs=5, vmin=vmin, vmax=vm,lut=lut, filename='fig_ramses_mcalens_sparse_sig5_masked.png')

    tvilut(kappa_map_noisy, title="CFIS kappa")
    k = kappa_map_noisy
    m = kappa_map * 0 + 1
    nx = 512
    ny = 512
    ns = 5
    m[100:200, :] = 0
    tvilut(m * kappa_map_noisy, title="mask")

    WT = starlet2d(gen2=False, l2norm=False, verb=False)
    WT.init_starlet(nx, ny, nscale=ns)
    H = HOS_starlet_l1norm_peaks(WT)
    SigmaMap = kappa_map * 0 + sigma_noise_CFIS

    MinSNR = -2
    MaxSNR = 6
    Nb = 31
    H.set_bins(Min=MinSNR, Max=MaxSNR, nbins=Nb)
    H.set_data(kappa_map_noisy, SigmaMap, Mask=m)
    H.get_wtpeaks(Mask=m)
    pc = H.Peaks_Count
    # plot_peak_count_histograms(H.TabBinsCenter, pc, 'Peak Counts Histogram for Gaussian SNR Map', 'SNR smooth', 'Peak Counts', log_scale=True)
    H.plot_peaks_histo(log_scale=True)
    H.plot_peaks_histo(Scale=0, log_scale=True, title="Scale 1")
    H.plot_peaks_histo(Scale=1, log_scale=True)

    H.tvpeaks(maxcoef=3, lut="rainbow", OnlyPeaks=True)
    H.get_wtl1(NBINS * 2)
    H.plot_l1norm()
    H.plot_l1norm(Scale=0)

    # Plot l1-norm histograms for different scales
    # plot_l1norm_histograms(H.l1bins, H.l1norm, 'L1-norm Histograms for Different Scales', 'L1-norm', 'Frequency')

    H.get_wtl1(NBINS * 2, Mask=m)
    # plot_l1norm_histograms(H.l1bins, H.l1norm, 'L1-norm Histograms for Different Scales', 'L1-norm', 'Frequency')
    H.plot_l1norm()


############ END CLASS #######################
