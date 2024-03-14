#  @cosmmostat_init.py
#
#  Some convenient ROUTINES for interactive use
#
#  Functions which makes the life easier for iSAP users.
#
#  @author Jean-Luc Starck
#  @version 1.0
#  @date 2020
#

import numpy as np
from os import remove
from subprocess import check_call
from datetime import datetime
import matplotlib.pyplot as plt
import astropy
from astropy.io import fits
from subprocess import check_call
import pylab
import readline
from mpl_toolkits.axes_grid1 import make_axes_locatable

# import pyqtgraph
# from pyqtgraph.Qt import QtGui
import numpy
from matplotlib.colors import PowerNorm
import seaborn as sns
import scipy
from scipy import ndimage
from scipy import signal

################################################

# def pwd():
# 	check_call(['pwd'])

################################################

# def ls():
# 	check_call(['ls'])

################################################
# similiar to the rebin function


def smooth2d(map, sigma):
    return ndimage.gaussian_filter(map, sigma=sigma)


def rebin2d(a, shape):
    sh = shape[0], a.shape[0] // shape[0], shape[1], a.shape[1] // shape[1]
    return a.reshape(sh).mean(-1).mean(1)


def rebin1d(a, shape):
    sh = shape[0], a.shape[0] // shape[0]
    return a.reshape(sh).mean(-1)


def pad_width2d(size):
    pad_width = ((size, size), (size, size))
    return pad_width


def unpad(x, pad_width):
    slices = []
    for c in pad_width:
        e = None if c[1] == 0 else -c[1]
        slices.append(slice(c[0], e))
    return x[tuple(slices)]


def unpad2d(x, size):
    return unpad(x, pad_width2d(size))


def gauss2d(size, fwhm=20, center=None):
    """Make a square gaussian kernel.

    size is the length of a side of the square
    fwhm is full-width-half-maximum, which
    can be thought of as an effective radius.
    """

    x = np.arange(0, size, 1, float)
    y = x[:, np.newaxis]

    if center is None:
        x0 = y0 = size // 2
    else:
        x0 = center[0]
        y0 = center[1]

    return np.exp(-4 * np.log(2) * ((x - x0) ** 2 + (y - y0) ** 2) / fwhm**2)


# Test
# pad_width = ((0, 0), (1, 0), (3, 4))
# a = np.random.rand(10, 10, 10)
# b = np.pad(a, pad_width, mode='constant')
# c = unpad(b, pad_width)
# np.testing.assert_allclose(a, c)

################################################


def hthres(alpha, Thres):
    Res = np.copy(alpha)
    Res[np.abs(Res) <= Thres] = 0
    return Res


def hard_thresholding(alpha, Thres):
    # alpha: shape = ([nimgs], [Nrea], ns, nx, ny)
    # Thres: shape = (ns, [nx, ny])
    if Thres.ndim == 1:
        Thres = Thres[..., np.newaxis, np.newaxis]
    alpha[np.abs(alpha) <= Thres] = 0


################################################


def soft_thresholding(alpha, Thres):
    # alpha: shape = ([nimgs], [Nrea], ns, nx, ny)
    # Thres: shape = (ns, [nx, ny])
    if Thres.ndim == 1:
        Thres = Thres[..., np.newaxis, np.newaxis]
    Res = np.abs(alpha) - Thres
    Res[Res < 0] = 0
    alpha[:, :] = np.sign(alpha) * Res


def sthres(alpha, Thres):
    Res = np.copy(alpha)
    Res = np.abs(alpha) - Thres
    Res[Res < 0] = 0
    Res = np.sign(alpha) * Res
    return Res


################################################


def info(Data, name=0, mask=None):
    if mask is None:
        if name:
            print(
                name,
                ": Size = ",
                Data.shape,
                "Type = ",
                Data.dtype,
                ", Mean = ",
                np.mean(Data),
                " Sigma = ",
                np.std(Data),
                "Min = ",
                np.min(Data),
                " Max = ",
                np.max(Data),
            )
        else:
            print(
                "Size = ",
                Data.shape,
                "Type = ",
                Data.dtype,
                ", Mean = ",
                np.mean(Data),
                " Sigma = ",
                np.std(Data),
                "Min = ",
                np.min(Data),
                " Max = ",
                np.max(Data),
            )
    else:
        ind = np.where(mask != 0)
        if name:
            print(
                name,
                ": Size = ",
                Data.shape,
                "Type = ",
                Data.dtype,
                ", Mean = ",
                np.mean(Data[ind]),
                " Sigma = ",
                np.std(Data[ind]),
                "Min = ",
                np.min(Data[ind]),
                " Max = ",
                np.max(Data[ind]),
            )
        else:
            print(
                "Size = ",
                Data.shape,
                "Type = ",
                Data.dtype,
                ", Mean = ",
                np.mean(Data[ind]),
                " Sigma = ",
                np.std(Data[ind]),
                "Min = ",
                np.min(Data[ind]),
                " Max = ",
                np.max(Data[ind]),
            )


################################################


def indgen(N):
    return np.arange(N)


def reverse(a):
    return np.flipud(a)


def dump(obj):
    for attr in dir(obj):
        print("obj.%s = %r" % (attr, getattr(obj, attr)))


def vsize(Data):
    # isinstance(P, (list, tuple, np.ndarray))
    if np.isscalar(Data):
        vs = np.zeros([1], dtype=np.int64)
    else:
        d = np.asarray(Data, dtype=np.int64)
        dim = d.ndim
        vs = np.zeros([dim + 1], dtype=np.int64)
        vs[0] = dim
        vs[1:] = np.array(d.shape)
    return vs


################################################


def readfits(FileName, plot=False, verbose=False):
    #    print("READ FITS2: "+FileName)
    fitsfile = fits.open(FileName, ignore_missing_end=True)
    if verbose:
        fitsfile.info()
    ima = fitsfile[0]
    data = ima.data
    if verbose:
        print(ima.header["NAXIS"])
        print("==>Size = ", vsize(data))
    Header = ima.header
    if plot:
        plt.imshow(data, origin="lower")
        plt.show()
    #        plt.close
    #    print("READ FITS1 END: "+FileName)
    return data


def writefits(FileName, Data):
    #    hdu = fits.PrimaryHDU(Data)
    #    hdul = fits.HDUList([hdu])
    fits.writeto(FileName, data=Data, overwrite=True)


################################################


def dft2d(ima):
    return np.fft.fftshift(np.fft.fft2(ima))


def idft2d(ima):
    return np.fft.ifft2(np.fft.ifftshift((ima)))


def idft2dr(ima):
    z = np.fft.ifft2(ima)
    return z.real


def conv(ima1, ima2):
    return scipy.signal.fftconvolve(ima1, ima2, mode="same")


#    x = np.fft.fft2(ima1) * np.fft.fft2(ima2)
#    z = np.fft.ifft2(x)
#    return z.real
#     return idft2dr(dft2d(ima1) * dft2d(ima2))


def dft2dnorm(ima):
    z = np.fft.fftshift(np.fft.fft2(ima))
    z = z * z.conj()
    return z.real


# fig, ax = plt.subplots(1, 1, facecolor='w', figsize=(FigSize,FigSize))
# fig, [ax] = plt.subplots(nrows=1, ncols=1, facecolor='w', figsize=(10, 10))

################################################


def tvilut(
    d,
    title=None,
    xtitle=None,
    ytitle=None,
    vmin=None,
    vmax=None,
    lut=None,
    fs=None,
    filename=None,
):
    if fs is None:
        FigSize = 10
    else:
        FigSize = fs
    if lut is None:
        lut = "rainbow"  # 'inferno'   'gist_stern'
    fig, ax = plt.subplots(1, 1, facecolor="w", figsize=(FigSize, FigSize))
    cm = plt.cm.get_cmap(lut)
    img = ax.imshow(d, cmap=cm, origin="lower", vmin=vmin, vmax=vmax)
    if title:
        ax.set_title(title)
    if xtitle:
        ax.set_xlabel(xtitle)
    if ytitle:
        ax.set_ylabel(ytitle)
        # ax.set_xlabel(r"$\Sigma$")
        # ax.set_ylabel("DEC")
        # ax.set_xlim(10, 40)
    # fig.colorbar(img)

    # cbar = fig.colorbar(cax, ticks=[-1, 0, 1])
    # cbar = fig.colorbar(cax, ticks=[d.min(), 0, d.max()])
    # cbar = fig.colorbar(cax)
    # cbar.ax.set_yticklabels(['< -1', '0', '> 1'])  # vertically oriented
    # colorbar

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    fig.colorbar(img, cax=cax)
    if filename is not None:
        plt.savefig(filename)
    fig.show()


#  	plt.imshow(d, origin='lower' )
# plt.colorbar()
# plt.show()


################################################


def tvimap(map_data, title="", lut="inferno", vmin=None, vmax=None, filename=None):
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
    img = plt.imshow(map_data, cmap=lut, vmin=vmin, vmax=vmax, origin="lower")
    plt.title(title)
    plt.colorbar(img)
    if filename is not None:
        plt.savefig(filename)
    plt.show()


################################################


def history():
    print(
        "\n".join(
            [
                str(readline.get_history_item(i))
                for i in range(readline.get_current_history_length())
            ]
        )
    )


################################################
def h():
    history()


################################################
def journal(FileName):
    readline.write_history_file(FileName)


################################################


def tvima(im, vmax=0, gamma=0.5, cmap="gist_stern"):
    if vmax == 0:
        vmax = im.max()
    # , norm=PowerNorm(gamma=gamma))
    plt.imshow(np.rot90(im, 2), cmap=cmap, vmax=vmax)
    plt.colorbar()
    plt.show()


################################################


def tvima2(ima1, ima2, tit1="data1", tit2="data2", vmin=None, vmax=None):
    """
    This routine plot 2 images

    Parameters
    ----------
    ima1, ima2 : np.ndarray
            Input data array

    Returns
    -------

    Notes
    -----
    plot the 2 images
    """
    if vmin is None:
        vmin = -0.1
    if vmax is None:
        vmax = 0.1
    plt.rcParams["figure.figsize"] = (15, 5)
    plt.rcParams.update({"font.size": 13})
    fig1, [ax1, ax2] = plt.subplots(nrows=1, ncols=2)
    cax1 = ax1.imshow(
        ima1, interpolation="none", cmap="inferno", origin="lower", vmin=vmin, vmax=vmax
    )
    ax1.set_title(tit1), plt.colorbar(cax1, ax=ax1)
    cax2 = ax2.imshow(
        ima2, interpolation="none", cmap="inferno", origin="lower", vmin=vmin, vmax=vmax
    )
    ax2.set_title(tit2), plt.colorbar(cax2, ax=ax2)
    plt.show()


################################################


def tvimacont(im, levels, vmin=0, vmax=0, gamma=0.5, cmap="gist_stern"):
    if vmin == 0:
        vmin = im.min()
    if vmax == 0:
        vmax = im.max()
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(111)
    ax.imshow(im, vmin=vmin, vmax=vmax, origin="lower")
    #    plt.colorbar()
    ax.contour(im, levels=levels, colors="white", alpha=0.5)
    plt.show()


################################################


def tv_frames(TabFrame, scales=None, multiview=False):
    """Display the different bands on the requested scales.

    Parameters
    ----------
    transform: WaveletTransformBase derived instance
        a wavelet decomposition.
    scales: list of int, default None
        the desired scales, if None compute at all scales.
    multiview: bool, default False
        if True use a slider to select a specific band.
    """
    # Set default scales
    Dim = np.shape(TabFrame)
    nb_scale = Dim[0]
    scales = scales or range(nb_scale)

    # Create application and tab widget
    app = pyqtgraph.mkQApp()
    tabs = QtGui.QTabWidget()
    tabs.setWindowTitle("Wavelet Transform")

    # Go through each scale
    pen = pyqtgraph.intColor(2)
    for scale in scales:
        # Create the plots for this scales with scrolling possibilities
        scroller = QtGui.QScrollArea()
        tabs.addTab(scroller, "Scale {0}".format(scale))

        # Go through each band of the current scale
        # > using multiview
        # TODO: update this code
        if multiview:
            raise NotImplementedError("Multiview transform view not yet implemented.")
            window = pyqtgraph.image(numpy.asarray(transform[scale]))
            scroller.setWidget(window)
        # > using mosaic
        else:
            window = pyqtgraph.GraphicsWindow()
            scroller.setWidget(window)

            scale_data = np.rot90(TabFrame[scale], 1)
            if not isinstance(scale_data, list):
                scale_data = [scale_data]

            for index, subband_data in enumerate(scale_data):
                # Deal with complex data
                if numpy.iscomplex(subband_data).any():
                    subband_data = numpy.abs(subband_data)
                subband_data = numpy.lib.pad(
                    subband_data, 1, "constant", constant_values=subband_data.max()
                )

                # Select viewer
                if subband_data.ndim == 1:
                    ax = window.addPlot()
                    ax.plot(subband_data, pen=pen)
                elif subband_data.ndim == 2:
                    row = index // 2
                    col = index % 2
                    box = window.addViewBox(
                        row=row,
                        col=col,
                        border="00ff00",
                        lockAspect=True,
                        enableMouse=False,
                    )
                    image = pyqtgraph.ImageItem(subband_data)
                    box.addItem(image)
                else:
                    raise ValueError(
                        "This function currently support only " "1D or 2D data."
                    )
                window.nextRow()

    # Display the tab
    tabs.show()

    # Run application
    app.exec_()


def tvnima(images, ncol=1, nline=None, title=None, scale=1, vmin=None, vmax=None):
    sns.set_style("ticks")
    n_images = len(images)
    if nline is None:
        nline = 1
    if ncol is None:
        ncol = n_images / nline + 1

    if title is None:
        title = ["Image (%d)" % i for i in range(1, n_images + 1)]
    fig = plt.figure(facecolor="white")  # ,figsize=(5,5))

    for n, (image, tit) in enumerate(zip(images, title)):
        # a = fig.add_subplot(nline+1,ncol,n+1)  #
        a = fig.add_subplot(ncol, np.ceil(n_images / float(ncol)) + 1, n + 1)
        cm = plt.cm.get_cmap("jet")
        plt.imshow(image, cmap=cm, vmin=vmin, vmax=vmax)
        a.set_title(tit)  # , plt.colorbar(cax2, ax=a)

    plt.colorbar()
    fig.set_size_inches(np.array(fig.get_size_inches()) * n_images * scale)
    plt.show()


def show_images(images, cols=1, titles=None, scale=1):
    sns.set_style("ticks")
    n_images = len(images)
    if titles is None:
        titles = ["Image (%d)" % i for i in range(1, n_images + 1)]
    fig = plt.figure(facecolor="white")
    for n, (image, title) in enumerate(zip(images, titles)):
        a = fig.add_subplot(cols, np.ceil(n_images / float(cols)) + 1, n + 1)
        plt.imshow(image)
        a.set_title(title)
    # plt.colorbar()
    fig.set_size_inches(np.array(fig.get_size_inches()) * n_images * scale)
    plt.show()
