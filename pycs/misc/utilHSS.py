"""
Created on Nov 4, 2015

@author: mjiang
"""

import numpy as np
import scipy.fftpack as pfft


def fftshift2d1d(cubefft):
    return np.fft.fftshift(cubefft, axes=(1, 2))


def ifftshift2d1d(cubefftSh):
    return np.fft.ifftshift(cubefftSh, axes=(1, 2))


def fft2d1d(cube):
    return np.fft.fft2(cube)


def ifft2d1d(cubefft):
    return np.fft.ifft2(cubefft)


def fftshiftNd1d(imagefft, N):
    assert len(imagefft.shape) == (N + 1)
    if N == 1:
        axes = 1
    if N == 2:
        axes = (1, 2)
    return np.fft.fftshift(imagefft, axes=axes)


def ifftshiftNd1d(imagefftSh, N):
    assert len(imagefftSh.shape) == (N + 1)
    if N == 1:
        axes = 1
    if N == 2:
        axes = (1, 2)
    return np.fft.ifftshift(imagefftSh, axes=axes)


def fftNd1d(image, N):
    assert len(image.shape) == (N + 1)
    if N == 1:
        out = np.fft.fft(image)
    if N == 2:
        out = np.fft.fft2(image)
    return out


def ifftNd1d(imagefft, N):
    assert len(imagefft.shape) == (N + 1)
    if N == 1:
        out = np.fft.ifft(imagefft)
    if N == 2:
        out = np.fft.ifft2(imagefft)
    return out


def mad(alpha):
    dim = np.size(np.shape(alpha))
    if dim == 1:
        alpha = alpha[np.newaxis, :]
    axes = tuple(range(1, dim))
    sigma = np.median(np.abs(alpha - np.median(alpha, axis=axes)), axis=axes) / 0.6745
    return sigma


def softTh(alpha, thTab, weights=None, reweighted=False):
    # TODO: remove for loop
    nz = np.size(thTab)
    #     print (weights.shape)
    #     print (thTab.shape)
    for i in np.arange(nz):
        if not reweighted:
            (alpha[i])[abs(alpha[i]) <= thTab[i]] = 0
            (alpha[i])[alpha[i] > 0] -= thTab[i]
            (alpha[i])[alpha[i] < 0] += thTab[i]
        else:
            (alpha[i])[np.abs(alpha[i]) <= thTab[i] * weights[i]] = 0
            (alpha[i])[alpha[i] > 0] -= (thTab[i] * weights[i])[alpha[i] > 0]
            (alpha[i])[alpha[i] < 0] += (thTab[i] * weights[i])[alpha[i] < 0]


def hardTh(alpha, thTab, weights=None, reweighted=False):
    # TODO: remove for loop
    nz = np.size(thTab)
    #     print (weights.shape)
    #     print (thTab.shape)
    if np.ndim(alpha) == 1:
        alpha = alpha[np.newaxis, :]
    for i in np.arange(nz):
        if not reweighted:
            (alpha[i])[abs(alpha[i]) <= thTab[i]] = 0
        else:
            (alpha[i])[np.abs(alpha[i]) <= thTab[i] * weights[i]] = 0
    alpha = alpha.squeeze()


def filter_Hi(sig, Ndim, fc, fc2=1.0 / 8):
    """
    The function is a high-pass filter applied on signals with fc as the cut-off frequency

    @param sig: 1D or 2D signal as entry, the number of sources is n(n>=1)

    @param Ndim: The dimension of the signal, Ndim=1 means 1D signal, Ndim=2 means 2D image

    @param fc: The cut-off frequency is given by normalized numerical frequency

    @return: The high frequency part of the signal
    """
    if Ndim == 1:
        (nz, ny) = np.shape(sig)
        sig_Hi = np.zeros_like(sig)
        sig_Hi[:, int(fc * ny) : -int(fc * ny)] = sig[:, int(fc * ny) : -int(fc * ny)]
    elif Ndim == 2:
        (nz, nx, ny) = np.shape(sig)
        sig_Hi = np.zeros_like(sig)
        sig_Hi[:, int(fc * nx) : -int(fc * nx), int(fc * ny) : -int(fc * ny)] = sig[
            :, int(fc * nx) : -int(fc * nx), int(fc * ny) : -int(fc * ny)
        ]
    return sig_Hi


def div0(a, b):
    """ignore / 0, div0( [-1, 0, 1], 0 ) -> [0, 0, 0]"""
    with np.errstate(divide="ignore", invalid="ignore"):
        c = np.true_divide(a, b)
        c[~np.isfinite(c)] = 0  # -inf inf NaN
    return c
