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
        It must have the same size as the size given in the starlet classs W.
    Mask: 2D array
        mask on the data. Mask[i,j] = 1 if the pixel is observed, 0 otherwise
    Returns
    -------
    3D array
        x = get_wt_noiselevel(W, NoiseStdMap):
        x[i,:,:] is the noise standard deviation map at the scale i.
    """
    if Mask is None:
        NoiseStdMap = NoiseSigmaMap
    else:
        # we will consider in the non observed area (mask == 0) that 
        # the noise is equivalent to the maximum noise level in the map.
        # This gives an over-estimation of the noise at these locations after 
        # inpaiting. 
        NoiseStdMap = NoiseSigmaMap
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
    VarianceCoeff = W.coef * 0.
    VarianceMap = NoiseStdMap**2
    for i in np.arange(WT.ns):
        coef = W.coef[i, :, :]
        coef = coef**2
        VarianceCoeff[i, :, :] = conv(VarianceMap,coef)
    return np.sqrt(VarianceCoeff)


def get_snr_noiselevel(W, Map, Wnoise, KeepSign=False):
    WT.transform(Map)
    WSnr = Wnoise * 0.
    ind = np.where(Wnoise != 0)
    WSnr[ind] = W.coef[ind] / Wnoise[ind]
    if KeepSign is False:
        WSnr = np.abs(WSnr)
    return WSnr


def get_peaks(image, threshold=None, ordered=True, mask=None,
                 include_border=False):
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
    From LensPack, added here to avoid circular dependancies

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
        image = np.pad(image, pad_width=1, mode='constant',
                       constant_values=image.min())
        mask = np.pad(mask, pad_width=1, mode='constant', constant_values=1)

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
    merge = ((map0 > map1) & (map0 > map2) & (map0 > map3) & (map0 > map4) &
             (map0 > map5) & (map0 > map6) & (map0 > map7) & (map0 > map8))

    bordered = np.lib.pad(merge, (1, 1), 'constant', constant_values=(0, 0))
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


class HOS_starlet_l1norm_peaks:
    NBins = 0
    TabBins = 0
    WT = 0
    Ns = 0
    DataCoef = 0
    NoiseCoef = 0
    WSNR = 0
    PIX_ARCMIN = 1
    scales_arcmin = 0
    TabBinsCenter = 0
    Peaks_PosX = 0
    Peaks_PosY = 0
    Peaks_Height = 0
    Peaks_Count = 0
    DataMask = 0
    l1bins = 0
    l1_coll = 0
    
    def __init__(self, WTrans, PixResolArcMin=1):  # __init__ is the constructor
        self.WT = WTrans  # Starlet wavelet Class defined in starlet.py
        self.Ns = WTrans.ns
        self.PIX_ARCMIN = PixResolArcMin
        self.scales_arcmin = [2**(i+1) * PIX_ARCMIN for i in range(ns)]

    def set_wtnoise(self, SigmaMap, Mask=None):
        NoiseCoef = get_wt_noiselevel(WT, SigmaMap, Mask=Mask)

    def set_bins(self, Min=-2, Max=6, nbins=31):
        self.NBins = nbins
        self.TabBins = np.linspace(Min, Max, self.NBins)
        self.TabBinsCenter = 0.5 * (self.TabBins[:-1] + self.TabBins[1:])
    
    def tvscalepeaks(self, ScaleNumber, maxcoef=None, lut='inferno', OnlyPeaks=None):
        if maxcoef is None:
            maxcoef = 2.
        i = ScaleNumber
        Scale = self.DataCoef[i,:,:]

        X = self.Peaks_PosX[i]
        Y = self.Peaks_PosY[i]
        if OnlyPeaks is not None:
            Scale[:,:] = 0
            Scale[X,Y] = 1
        else:
         Max = Scale.max()
         Scale[X,Y] = maxcoef*Max
        tvimap(Scale,title='Scale'+str(i+1),lut=lut)
        
    def tvpeaks(self, maxcoef=None, lut='inferno', OnlyPeaks=None):
        if maxcoef is None:
            maxcoef = 2.
        for i in np.arange(self.WT.ns):
            self.tvscalepeaks(i,maxcoef=maxcoef, lut=lut, OnlyPeaks=OnlyPeaks)
            
    def get_wtpeaks(self, Map, SigmaMap, Mask=None):
        self.DataMask = Mask
        self.NoiseCoef = get_wt_noiselevel(self.WT, SigmaMap, Mask=Mask)
        self.WSNR = get_snr_noiselevel(self.WT, Map, self.NoiseCoef, KeepSign=True)
        self.DataCoef = self.WT.coef
        self.Peaks_Count = []
        self.Peaks_PosX = []
        self.Peaks_PosY = []
        self.Peaks_Height = []
        for i in np.arange(self.WT.ns):
            Scale = self.WSNR[i,:,:]
            Min = Scale.min()
            Max = Scale.max()
            print("Scale ",i+1, ": Min = ", Min, ", Max = ", Max)
            X, Y, heights = get_peaks(Scale, threshold=None, ordered=True, mask=Mask,
                         include_border=False)
            self.Peaks_PosX.append(X)
            self.Peaks_PosY.append(Y)
            self.Peaks_Height.append(heights)
            counts, bin_edges = np.histogram(heights, bins=self.TabBins)
            self.Peaks_Count.append(counts)
        
    def plot_peaks_histo(self, Scale=None, title='Starlet Peak Counts Histogram', xlabel='SNR', ylabel='Peak Counts', log_scale=False):
        """
        Plot histograms of peak counts.

        Parameters:
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
                    plt.plot(self.TabBinsCenter, pc, label=f'Scale {scale + 1}')
         
        plt.legend()
        if log_scale:
            plt.yscale('log')
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.grid(True)
        plt.title(title)
        plt.show()

    def plot_l1norm(self, Scale=None, xlim=None, title='Starlet l1-norm', xlabel='SNR', ylabel='l1-norm', log_scale=False):
        plt.figure()

        if Scale is not None:
            l1 = self.l1norm[Scale]
            plt.plot(self.l1bins[Scale], l1, label=f'Scale {Scale}')
        else:
            for scale, l1 in enumerate(self.l1norm):
                plt.plot(self.l1bins[scale], l1, label=f'Scale {scale}')
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
            plt.yscale('log')
        plt.show()
        
    def get_wtl1(self, nbins=None, Mask=None):
        if nbins is None:
            nbins = NBins
        l1_coll = []
        bins_coll = []
        for i in np.arange(self.WT.ns):
            ScaleSNR = self.WSNR[i,:,:]
            if Mask is not None:
                ind = np.where(Mask == 0)
                ScaleSNR = ScaleSNR[ind]
            thresholds_snr = np.linspace(np.min(ScaleSNR), np.max(ScaleSNR), nbins + 1)
            bins_snr = 0.5 * (thresholds_snr[:-1] + thresholds_snr[1:])
            digitized = np.digitize(ScaleSNR, thresholds_snr)
            bin_l1_norm = [np.sum(np.abs(ScaleSNR[digitized == j]))
                           for j in range(1, len(thresholds_snr))]
            bins_coll.append(bins_snr)
            l1_coll.append(bin_l1_norm)
        self.l1bins = np.array(bins_coll)
        self.l1norm = np.array(l1_coll)

#     """High order statistics class for wavelet peaks and l1 norm
#     """

#     UseMask = False
#     Mask = 0
#     PIX_ARCMIN = 1
    
#     def __init__(self, name="HOS", PixArmin=None):  # __init__ is the constructor
#         if PixArmin is not None:
#             self.PIX_ARCMIN = PixArmin

#     def set_linear_bins(Min, Max, Nbins=10):
#         # Define bin ed

def tvimap(map_data, title='', lut='inferno', vmin=None, vmax=None,filename=None):
    """
    Plot a 2D map using a colormap.

    Parameters:
        map_data (numpy.ndarray): The 2D map data.
        title (str): Title of the plot.
        lut (str): Colormap name ('rainbow','inferno', 'gist_stern', etc)
        vmin (float): Minimum value for colormap scaling.
        vmax (float): Maximum value for colormap scaling.
    """
    plt.figure()
    img = plt.imshow(map_data, cmap=lut, vmin=vmin, vmax=vmax, origin='lower')
    plt.title(title)
    plt.colorbar(img)
    if filename is not None:
        plt.savefig(filename)
    plt.show()  

############ END CLASS #######################

# if __name__ == '__main__':
#     print ( "Main :)")

# Exemple for the CFIS experiment
kappa_map = readfits("/Users/starck/tex/PPT/Y23/EuclidTutorialSchool2023/WLconv_z0.50_3316r.fits")
SHAPE_NOISE = 0.44
PIX_ARCMIN = 0.4
N_GAL = 7
NSCALES = 5
NBINS = 40 
KAPPA_SNR = np.linspace(-2, 6, 31)

# Make a noise 
sigma_noise_CFIS = SHAPE_NOISE / (np.sqrt(2 * N_GAL * PIX_ARCMIN**2))
noise_map_CFIS_z05 = sigma_noise_CFIS * np.random.randn(kappa_map.shape[0], kappa_map.shape[1]) # generate noise map
kappa_map_noisy = kappa_map + noise_map_CFIS_z05 # Add noise to the mass map



; tvilut(k2r5*mask, title='CFIS kappa',fs=5, vmin=vmin, vmax=vm,lut=lut, filename='fig_ramses_mcalens_sparse_sig5_masked.png')

tvilut(kappa_map_noisy, title='CFIS kappa')
k = kappa_map_noisy
m = kappa_map * 0 + 1
nx=512
ny=512
ns=5
m [100:200,:] = 0
tvilut(m*kappa_map_noisy, title='mask')

WT = starlet2d(gen2=False,l2norm=False, verb=False)
WT.init_starlet(nx, ny, nscale=ns)
H = HOS_starlet_l1norm_peaks(WT)
SigmaMap = kappa_map * 0 + sigma_noise_CFIS


MinSNR = -2
MaxSNR = 6
Nb=31
H.set_bins(Min=MinSNR, Max=MaxSNR, nbins=Nb)
H.get_wtpeaks(kappa_map_noisy, SigmaMap, Mask=m)
pc = H.Peaks_Count
# plot_peak_count_histograms(H.TabBinsCenter, pc, 'Peak Counts Histogram for Gaussian SNR Map', 'SNR smooth', 'Peak Counts', log_scale=True)
H.plot_peaks_histo(log_scale=True)
H.plot_peaks_histo(Scale=0, log_scale=True, title='Scale 1')
H.plot_peaks_histo(Scale=1, log_scale=True)

H.tvpeaks(maxcoef=3, lut='rainbow',OnlyPeaks=True)
H.get_wtl1(NBINS*2)
H.plot_l1norm()
H.plot_l1norm(Scale=0)
              

# Plot l1-norm histograms for different scales
# plot_l1norm_histograms(H.l1bins, H.l1norm, 'L1-norm Histograms for Different Scales', 'L1-norm', 'Frequency')

H.get_wtl1(NBINS*2,Mask=m)
# plot_l1norm_histograms(H.l1bins, H.l1norm, 'L1-norm Histograms for Different Scales', 'L1-norm', 'Frequency')
H.plot_l1norm()

# for i in np.arange(WT.ns):
#     print(i)
#     histo = pc[i]
#     print(histo)
#     plot_peak_count_histograms(H.TabBinsCenter, histo, 'Peak Counts Histogram for Gaussian SNR Map', 'SNR smooth', 'Peak Counts', log_scale=True)
    



# print("Result")
# NoiseMap = k * 0. + 10.
# Wtest = get_wt_noiselevel(WT, NoiseMap)
# for i in np.arange(WT.ns):
#     print("Scale ", i)
#     info(Wtest[i, :, :])


# for i in np.arange(WT.ns):
#     # wt[i, :, :] *= WT.Starlet_Gen1TabNorm[i]
#     print(WT.Starlet_Gen1TabNorm[i])


# g1 = np.zeros((64,64))
# g2 = g1
# nx=64
# ny=64
# g1[32,32]= 1
# g2[32,32]= 1
# m = g1 * 0 + 1
# g1 = np.random.normal(loc=0.0, scale=np.sqrt(1. / 1.0)*m)
# g2 = np.random.normal(loc=0.0, scale=np.sqrt(1. / 1.0)*m)

# nx=256
# ny=256
# M = massmap2d(name='mass')
# M.init_massmap(nx,ny)    
# k = M.gamma_to_kappa(g1, g2)

# Use_Nrea=20
# M.get_wt_noise_level(InshearData, Nrea=Use_Nrea)
# info(k)
# tvilut(k)







