# -*- coding: utf-8 -*-

import os
import tempfile
import numpy as np
from shutil import which
from subprocess import call
from astropy.io import fits
from scipy.ndimage import convolve1d
from scipy.fftpack import dct, idct


def mr_transform(image, nscales=4, type=2, verbose=False):
    """Compute the multi-resolution wavelet transform of an image.

    Parameters
    ----------
    image : array_like, 2D
        Input image.
    nscales : int, optional
        Number of wavelet scales to compute. Default is 4.
    type : int, optional
        Type of the multiresolution transform. See the original mr_transform
        documentation for details. Default is 2, which corresponds to the
        'bspline wavelet transform: a trous algorithm', i.e. the starlet.
    verbose : bool, optional
        If True, print details of the temporary file I/O process.

    Returns
    -------
    3D numpy array
        Result of the wavelet transform.

    Notes
    -----
    This function is a wrapper for the mr_transform C++ binary of the iSAP
    code package (see References). The astropy package is necessary to write
    out `image` as a temporary fits file on which mr_transform can act.

    References
    ----------
    * Starck, Murtagh, Fadili, 'Sparse Image and Signal Processing: Wavelets
      and Related Geometric Multiscale Analysis', Cambridge University Press,
      Cambridge (GB), 2016.
    * http://www.cosmostat.org/software/isap

    Examples
    --------
    ...

    """
    # Verify that mr_transform is installed
    assert which("mr_transform"), "Cannot find mr_transform. Is it installed?"

    # Create a temporary directory to hold the image and its transform
    tmpdir = tempfile.mkdtemp()
    saved_umask = os.umask(0o077)
    image_path = os.path.join(tmpdir, "image.fits")
    mr_path = os.path.join(tmpdir, "image.mr")
    if verbose:
        print("\nCreating {}".format(image_path))
        print("         {}".format(mr_path))

    # Call mr_transform on the saved image
    try:
        if verbose:
            print("Writing image to fits.")
        fits.writeto(image_path, image)
        callstr = [
            "mr_transform",
            "-t",
            str(type),
            "-n",
            str(nscales + 1),
            image_path,
            mr_path,
        ]
        if verbose:
            print("Executing " + " ".join(callstr))
        call(callstr)
        mr = fits.getdata(mr_path)
    except IOError as e:
        print("Something went wrong... trying again.")
        # print("Removing {}".format(image_path))
        os.remove(image_path)
        # print("Unmasking {}".format(tmpdir))
        os.umask(saved_umask)
        # print("Removing {}".format(tmpdir))
        os.rmdir(tmpdir)
        # print("Calling mr_transform.")
        return mr_transform(image, nscales=nscales, type=type, verbose=verbose)
    else:
        # If successful, remove file paths
        os.remove(image_path)
        os.remove(mr_path)
        if verbose:
            print("Success.")
            print("Removed {}".format(image_path))
            print("        {}".format(mr_path))
        # Remove temporary directory
        os.umask(saved_umask)
        os.rmdir(tmpdir)
        if verbose:
            print("        {}".format(tmpdir))

    if os.path.exists(tmpdir) or os.path.exists(image_path) or os.path.exists(mr_path):
        print("Warning : not all files or directories were removed in")
        print(mr_path)

    return mr
