##########################################################################
# XXX - Copyright (C) XXX, 2017
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# System import
import numpy


def fftconvolve(array, kernel):
    """Function that convolves an image with a kernal using FFT.

    Parameters
    ----------
    image: array
        a 2D image data.
    kernel: array
        a 2D kernel.

    Returns
    -------
    out: array
        the convolved image.
    """
    x = numpy.fft.fftshift(numpy.fft.fftn(image))
    y = numpy.fft.fftshift(numpy.fft.fftn(kernel))

    return numpy.real(numpynp.fft.fftshift(numpy.fft.ifftn(numpy.fft.ifftshift(x * y))))


def fftdeconvolve(image, kernel):
    """Function that deconvolves and image with a kernal using FFT.

    Parameters
    ----------
    image: array
        a 2D image data.
    kernel: array
        a 2D kernel.

    Returns
    -------
    out: array
        the deconvolved image.
    """
    x = numpy.fft.fftshift(numpy.fft.fftn(image))
    y = numpy.fft.fftshift(numpy.fft.fftn(kernel))

    return numpy.real(numpy.fft.fftshift(numpy.fft.ifftn(numpy.fft.ifftshift(x / y))))
