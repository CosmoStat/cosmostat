"""
Created on Mar 30, 2015

@author: Ming  Jiang and Jean-Luc Starck

Routines for GMCA evaluation
"""

import numpy as np
from os import remove
from subprocess import check_call
from subprocess import call
from datetime import datetime
from astropy.io import fits
import shlex
from pycs.misc.cosmostat_init import *
from pycs.misc.cosmostat_init import writefits
from skimage import data, color
from skimage.transform import resize

def load_source_images(shape=(128, 128)):
    # Load grayscale source images and resize
    img1 = resize(data.camera(), shape, anti_aliasing=True)
    img2 = color.rgb2gray(resize(data.astronaut(), shape, anti_aliasing=True))
    
    # Normalize images
    img1 = (img1 - np.mean(img1)) / np.std(img1)
    img2 = (img2 - np.mean(img2)) / np.std(img2)

    # Stack into a single array: (2, H, W)
    sources = np.stack([img1, img2], axis=0)
    return sources

def mix_sources_images(sources):
    n_sources, H, W = sources.shape
    n_pixels = H * W
    
    # Flatten sources: (2, H*W)
    S = sources.reshape(n_sources, n_pixels)

    # Create random mixing matrix: (3, 2)
    A = np.random.randn(3, 2)
    
    # Mix: (3, H*W)
    mixed = A @ S

    # Reshape to image form: (3, H, W)
    mixed_images = mixed.reshape(3, H, W)
    return mixed_images, A


def mix_sources_images_noise(sources, noise_level=0.05):
    n_sources, H, W = sources.shape
    n_pixels = H * W

    # Flatten source images
    S = sources.reshape(n_sources, n_pixels)

    # Generate random mixing matrix
    A = np.random.randn(3, 2)

    # Mix sources
    mixed = A @ S  # Shape: (3, n_pixels)

    # Add Gaussian noise
    noise = np.random.normal(scale=noise_level, size=mixed.shape)
    mixed += noise
    
    # Reshape back to images
    mixed_images = mixed.reshape(3, H, W)
    return mixed_images, A


def reorder_and_fix_sign(true_sources, estimated_sources):
    """
    Reorder and apply sign correction to estimated_sources so they match true_sources.

    Parameters:
        true_sources: np.ndarray of shape (n_sources, H, W)
        estimated_sources: np.ndarray of shape (n_sources, H, W)

    Returns:
        corrected_sources: np.ndarray of shape (n_sources, H, W)
    """
    n_sources, H, W = true_sources.shape
    S_true = true_sources.reshape(n_sources, -1)
    S_est = estimated_sources.reshape(n_sources, -1)

    # Normalize
    S_true = (S_true - S_true.mean(axis=1, keepdims=True))
    S_true /= np.linalg.norm(S_true, axis=1, keepdims=True)
    S_est = (S_est - S_est.mean(axis=1, keepdims=True))
    S_est /= np.linalg.norm(S_est, axis=1, keepdims=True)

    # Correlation matrix
    corr = S_true @ S_est.T  # (n_true, n_est)

    # Reorder and sign-correct
    used = set()
    corrected_sources = np.zeros_like(true_sources)

    for i in range(n_sources):
        idx = np.argmax(np.abs(corr[i]))
        while idx in used:
            corr[i, idx] = 0
            idx = np.argmax(np.abs(corr[i]))
        used.add(idx)
        sign = np.sign(corr[i, idx])
        corrected_sources[i] = sign * estimated_sources[idx]

    return corrected_sources

def compute_sdr(true_sources, estimated_sources):
    """
    Compute SDR for each pair of true and estimated sources.

    Parameters:
        true_sources: np.ndarray of shape (n_sources, H, W)
        estimated_sources: np.ndarray of shape (n_sources, H, W)

    Returns:
        sdr_values: list of SDR values for each source
    """
    sdr_values = []
    for i in range(true_sources.shape[0]):
        s_true = true_sources[i].flatten()
        s_est = estimated_sources[i].flatten()
        noise = s_true - s_est
        sdr = 10 * np.log10(np.sum(s_true ** 2) / np.sum(noise ** 2))
        sdr_values.append(sdr)
    return sdr_values


def amari_error(A_true, A_est):
    """
    Compute the Amari error between the true and estimated mixing matrices.

    Parameters:
        A_true: np.ndarray (n_obs, n_sources) — ground truth mixing matrix
        A_est: np.ndarray (n_obs, n_sources) — estimated mixing matrix

    Returns:
        Amari error (float)
    """
    try:
        # Estimate the unmixing matrix
        W_est = np.linalg.pinv(A_est)  # shape: (n_sources, n_obs)
        G = W_est @ A_true             # shape: (n_sources, n_sources)
    except np.linalg.LinAlgError:
        return np.inf

    # DEBUG
    if G.shape[0] != G.shape[1]:
        raise ValueError(f"G should be square, but got shape {G.shape}. Check input shapes.")

    G = np.abs(G)
    row_sums = np.sum(G, axis=1, keepdims=True)
    col_sums = np.sum(G, axis=0, keepdims=True)

    row_error = np.sum(np.sum(G / row_sums, axis=1) - 1)
    col_error = np.sum(np.sum(G / col_sums, axis=0) - 1)

    return (row_error + col_error) / (2 * G.shape[0])

# Metrics

def evaluate(A0, S0, A, S, corrPerm=False):
    """Computes the NMSE and the CA.

    Parameters
    ----------
    A0: np.ndarray
        (m,n) float array, ground truth mixing matrix
    S0: np.ndarray
        (n,p) float array, ground truth sources
    A: np.ndarray
        (m,n) float array, estimated mixing matrix
    S: np.ndarray
        (n,p) float array, estimated sources
    corrPerm: bool
        correct permutation of A and S (in-place updates)
    perScale: bool
        calculate NMSE per wavelet scale
    nscales: int
        number of wavelet detail scales
    S0wt: np.ndarray
        (m,n,nscales+1) float array, wavelet transform of S0, optional (to accelerate)

    Returns
    -------
    (float,float) or (float,float,np.ndarray)
        CA,
        NMSE,
        NMSE per scale if perScale ((nscales+1,) float array)
    """

    if not corrPerm:
        A = A.copy()
        S = S.copy()

    n = np.shape(A0)[1]

    corr_perm(A0, S0, A, S, inplace=True)

    # CA = -10 * np.log10(np.mean(np.abs(np.dot(np.linalg.pinv(A), A0) - np.eye(n))))
    CA = (np.mean(np.abs(np.dot(np.linalg.pinv(A), A0) - np.eye(n))))
    # NMSE = -10 * np.log10(np.sum((S0-S)**2)/np.sum(S0**2))
    NMSE =  (np.sum((S0-S)**2)/np.sum(S0**2))

    return CA, NMSE



def corr_perm(A0, S0, A, S, inplace=False, optInd=False):
    """Correct the permutation of the solution.

    Parameters
    ----------
    A0: np.ndarray
        (m,n) float array, ground truth mixing matrix
    S0: np.ndarray
        (n,p) float array, ground truth sources
    A: np.ndarray
        (m,n) float array, estimated mixing matrix
    S: np.ndarray
        (n,p) float array, estimated sources
    inplace: bool
        in-place update of A and S
    optInd: bool
        return permutation

    Returns
    -------
    None or np.ndarray or (np.ndarray,np.ndarray) or (np.ndarray,np.ndarray,np.ndarray)
        A (if not inplace),
        S (if not inplace),
        ind (if optInd)
    """

    A0 = A0.copy()
    S0 = S0.copy()
    if not inplace:
        A = A.copy()
        S = S.copy()

    n = np.shape(A0)[1]

    for i in range(0, n):
        S[i, :] *= (1e-24 + np.linalg.norm(A[:, i]))
        A[:, i] /= (1e-24 + np.linalg.norm(A[:, i]))
        S0[i, :] *= (1e-24 + np.linalg.norm(A0[:, i]))
        A0[:, i] /= (1e-24 + np.linalg.norm(A0[:, i]))

    try:
        diff = abs(np.dot(np.linalg.inv(np.dot(A0.T, A0)), np.dot(A0.T, A)))
    except np.linalg.LinAlgError:
        diff = abs(np.dot(np.linalg.pinv(A0), A))
        print('Warning! Pseudo-inverse used.')

    ind = np.arange(0, n)

    for i in range(0, n):
        ind[i] = np.where(diff[i, :] == max(diff[i, :]))[0][0]

    A[:] = A[:, ind.astype(int)]
    S[:] = S[ind.astype(int), :]

    for i in range(0, n):
        p = np.sum(S[i, :] * S0[i, :])
        if p < 0:
            S[i, :] = -S[i, :]
            A[:, i] = -A[:, i]

    if inplace and not optInd:
        return None
    elif inplace and optInd:
        return ind
    elif not optInd:
        return A, S
    else:
        return A, S, ind


def nmse(S0, S):
    """Compute the normalized mean square error (NMSE) in dB.

    Parameters
    ----------
    S0: np.ndarray
        (n,p) float array, ground truth sources
    S: np.ndarray
        (n,p) float array, estimated sources

    Returns
    -------
    float
        NMSE (dB)
    """
    return -10 * np.log10(np.sum((S0-S)**2)/np.sum(S0**2))


def ca(A0, A):
    """Compute the criterion on A (CA) in dB.

    Parameters
    ----------
    A0: np.ndarray
        (m,n) float array, ground truth mixing matrix
    A: np.ndarray
        (m,n) float array, estimated mixing matrix

    Returns
    -------
    float
        CA (dB)
    """
    return -10 * np.log10(np.mean(np.abs(np.dot(np.linalg.pinv(A), A0) - np.eye(np.shape(A0)[1]))))




