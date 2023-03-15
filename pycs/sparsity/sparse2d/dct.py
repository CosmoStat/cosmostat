# -*- coding: utf-8 -*-

import os
import tempfile
import numpy as np
from shutil import which
from subprocess import call
from astropy.io import fits
from scipy.ndimage import convolve1d
from scipy.fftpack import dct, idct


def dct2d(image, norm="ortho"):
    """Compute the discrete cosine transform (type 2) of an image.

    Parameters
    ----------
    image : array_like, 2D
        Input image.
    norm : {None, 'ortho', 'isap'}, optional
        Normalization option. See scipy.fftpack.dct documentation (Type II)
        for a description of the None and 'ortho' options. The 'isap' option
        is available to match the output from the im_dct iSAP binary, which
        involves an additional scaling of the zero-frequency elements.
        Default is 'ortho'.

    Returns
    -------
    2D numpy array
        Type 2 DCT.

    Notes
    -----
    Using no normalization (i.e. norm=None) will not automatically
    recover the original image after performing the inverse transformation.
    Each transform brings an overall scaling factor of 2N.

    See Also
    --------
    idct2d
        Inverse 2D DCT.

    Examples
    --------
    ...

    """
    # Check inputs
    image = np.array(image)
    assert len(image.shape) == 2, "Input image must be 2D."
    assert norm in (None, "ortho", "isap"), "Invalid norm."

    # Compute DCT along each axis
    if norm == "isap":
        result = dct(dct(image, norm="ortho", axis=0), norm="ortho", axis=1)
        result[:, 0] *= np.sqrt(2)
        result[0, :] *= np.sqrt(2)
    else:
        result = dct(dct(image, norm=norm, axis=0), norm=norm, axis=1)

    return result


def idct2d(image, norm="ortho"):
    """Compute the inverse discrete cosine transform (type 2) of an image.

    Parameters
    ----------
    image : array_like (2D)
        Input image.
    norm : {None, 'ortho', 'isap'}, optional
        Normalization option. Default is 'ortho'.

    Returns
    -------
    2D numpy array
        Inverse type 2 DCT.

    See Also
    --------
    dct2d
        Forward 2D DCT.

    Examples
    --------
    ...

    """
    # Check inputs
    image = np.array(image)
    assert len(image.shape) == 2, "Input image must be 2D."
    assert norm in (None, "ortho", "isap"), "Invalid norm."

    # Compute inverse DCT along each axis
    if norm == "isap":
        image[:, 0] /= np.sqrt(2)
        image[0, :] /= np.sqrt(2)
        result = idct(idct(image, norm="ortho", axis=0), norm="ortho", axis=1)
    else:
        result = idct(idct(image, norm=norm, axis=0), norm=norm, axis=1)

    return result


def blockdct2d(image, norm="ortho", blocksize=None, overlap=False):
    """Compute a block (local) discrete cosine transform of an image.

    This is an extension of dct2d to perform the transform on sub-blocks
    of the image.

    Parameters
    ----------
    image : array_like, 2D
        Input image.
    norm : {None, 'ortho', 'isap'}, optional
        Normalization option. See scipy.fftpack.dct documentation (Type II)
        for a description of the None and 'ortho' options. The 'isap' option
        is available to match the output from the im_dct iSAP binary, which
        involves an additional scaling of the zero-frequency elements.
        Default is 'ortho'.
    blocksize : int, optional
        Size of sub-blocks for a local DCT.
    overlap : bool, optional
        Whether to overlap sub-blocks.

    Returns
    -------
    2D numpy array
        Local type 2 DCT.

    See Also
    --------
    iblockdct2d
        Inverse local 2D DCT.

    Examples
    --------
    ...

    TODO
    -----
    This needs MORE TESTING before deployment !

    """
    # Check inputs
    image = np.array(image)
    assert len(image.shape) == 2, "Input image must be 2D."
    assert image.shape[0] == image.shape[1], "Input image must be square."
    assert norm in (None, "ortho", "isap"), "Invalid norm."

    # Determine output shape based on blocksize
    n = image.shape[0]
    if blocksize is not None:
        if blocksize == n:
            result = np.zeros_like(image)
        elif blocksize not in [n / 2, n / 4, n / 8]:
            print("Warning: invalid blocksize --> using {}".format(n))
            blocksize = n
            result = np.zeros_like(image)
        else:
            if overlap:
                size = 2 * n - blocksize
                result = np.zeros((size, size))
            else:
                result = np.zeros_like(image)
    else:
        blocksize = n
        result = np.zeros_like(image)

    print(blocksize)

    # Compute DCT on sub blocks
    if overlap:
        for ii in range(2 * n / blocksize - 1):
            for jj in range(2 * n / blocksize - 1):
                i1 = ii * blocksize
                i2 = i1 + blocksize
                j1 = jj * blocksize
                j2 = j1 + blocksize
                imsub = image[i1 / 2 : i1 / 2 + blocksize, j1 / 2 : j1 / 2 + blocksize]
                result[i1:i2, j1:j2] = dct2d(imsub, norm=norm)
    else:
        for ii in range(0, n, blocksize):
            for jj in range(0, n, blocksize):
                i1 = ii
                i2 = ii + blocksize
                j1 = jj
                j2 = jj + blocksize
                imsub = image[i1:i2, j1:j2]
                result[i1:i2, j1:j2] = dct2d(imsub, norm=norm)

    return result


def iblockdct2d(image, norm="ortho", blocksize=None, overlap=False):
    """Compute the inverse block (local) discrete cosine transform of an image.

    This is an extension of idct2d to perform the transform on sub-blocks
    of the image.

    Parameters
    ----------
    image : array_like, 2D
        Input image.
    norm : {None, 'ortho', 'isap'}, optional
        Normalization option. See scipy.fftpack.dct documentation (Type II)
        for a description of the None and 'ortho' options. The 'isap' option
        is available to match the output from the im_dct iSAP binary, which
        involves an additional scaling of the zero-frequency elements.
        Default is 'ortho'.
    blocksize : int, optional
        Size of sub-blocks for a local inverse DCT.
    overlap : bool, optional
        Whether to overlap sub-blocks.

    Returns
    -------
    2D numpy array
        Local type 2 inverse DCT.

    Examples
    --------
    ...

    TODO
    -----
        This needs MORE TESTING before deployment !

    """
    if norm not in [None, "ortho", "isap"]:
        print("Warning: invalid norm --> using isap")
        norm = "isap"

    # Determine output shape
    n = image.shape[0]
    if blocksize is not None:
        if blocksize == n:
            result = np.zeros_like(image)
        else:
            if overlap:
                size = (n + blocksize) / 2
                result = np.zeros((size, size))
            else:
                result = np.zeros_like(image)
    else:
        blocksize = n
        result = np.zeros_like(image)

    # Compute inverse DCT on sub blocks
    if overlap:
        for ii in range(n / blocksize):
            for jj in range(n / blocksize):
                i1 = ii * blocksize
                i2 = i1 + blocksize
                j1 = jj * blocksize
                j2 = j1 + blocksize
                i1r = i1 / 2
                i2r = i1r + blocksize
                j1r = j1 / 2
                j2r = j1r + blocksize
                imsub = image[i1:i2, j1:j2]
                result[i1r:i2r, j1r:j2r] += idct(imsub, norm=norm)

        # Take averages
        step = blocksize / 2
        counts = np.ones_like(result)
        counts[step:-step, :step] = 2
        counts[step:-step, -step:] = 2
        counts[:step, step:-step] = 2
        counts[-step:, step:-step] = 2
        counts[step:-step, step:-step] = 4
        result /= counts
    else:
        for ii in range(0, n, blocksize):
            for jj in range(0, n, blocksize):
                i1 = ii
                i2 = ii + blocksize
                j1 = jj
                j2 = jj + blocksize
                imsub = image[i1:i2, j1:j2]
                result[i1:i2, j1:j2] = idct(imsub, norm=norm)

    return result
