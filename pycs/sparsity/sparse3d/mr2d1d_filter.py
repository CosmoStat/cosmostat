#!/usr/bin/env python
# coding: utf-8

# In[1]:


import matplotlib.pyplot as plt

get_ipython().run_line_magic("matplotlib", "inline")

import numpy as np
import pysparse
from pycs.misc.cosmostat_init import *
from pycs.misc.stats import *

from scipy.ndimage.filters import gaussian_filter


# ## MR2D1D convenience class
#
# The Wavelet2D1DTransform class below is a wrapper for MR2D1D which includes convenience methods for decomposing / reconstructing a data cube. It also has a method for extracting the meta information from the output transform that specifies the sub-band coefficient indicies (locations in the flattened output array) and their 3D shape parameters.

# In[2]:


class Wavelet2D1DTransform(object):
    """Wavelet decomposition of a 3D data cube."""

    # NOT CURRENTLY USED
    # TODO Compute values more accurately by excluding the borders of each sub-band cube
    NOISE_TAB = np.array(
        [
            [0.9610, 0.9261, 0.9611, 0.9810, 0.9933],
            [0.2368, 0.2282, 0.2369, 0.2417, 0.2449],
            [0.1049, 0.1011, 0.1050, 0.1071, 0.1082],
            [0.0527, 0.0507, 0.0528, 0.0539, 0.0543],
            [0.0283, 0.0272, 0.0283, 0.0288, 0.0295],
        ]
    )

    def __init__(self, transform_type=2):
        """Wrapper for pysparse's 2D-1D wavelet transform.

        Parameters
        ----------
        transform_type : int
            Type of wavelet transform to perform with `decompose`. See pysap
            documentation for all options. Default is 2, which takes a starlet
            transform in the 2D domain (undecimated) and a 7/9 filter bank
            in the 1D domain (decimated).

        """
        self.transform_type = transform_type

    def decompose(self, cube, num_scales_2d, num_scales_1d):
        """Forward 2D1D wavelet transform.

        Parameters
        ----------
        cube : 3D array
            Data to transform.
        num_scales_2d : int
            Number of wavelet scales in the 2D domain (i.e. the last two axes).
        num_scales_1d : int
            Number of wavelet scales in the 1D domain (i.e. first axis).

        Returns
        -------
        inds : nested list of tuples
            Pairs of index values arranged such that inds[i][j] gives the
            (start, end) indices for the band (2d_scale, 1d_scale) = (i, j).
            See `_extract_index_ranges` docstring.
        coeffs : 1D array
            A flattened array object containing all the wavelet coefficients in
            increasing order of scale. Also included are the sizes of the
            transformed cube bands along each axis as three ints just before
            the coefficients themselves. Unpacking this array to access
            specific coefficients for a given (i, j) band requires `inds`.

        """
        # Compute the transform
        self._mr2d1d = pysparse.MR2D1D(
            type_of_transform=self.transform_type,
            normalize=False,
            verbose=True,
            NbrScale2d=num_scales_2d,
            Nbr_Plan=num_scales_1d,
        )
        coeffs = self._mr2d1d.transform(cube)

        # Determine the starting/ending index values of the 1d array that
        # correspond to the wavelet coefficients of band (scale_2d, scale_1d)
        inds, shapes = self._extract_metadata(coeffs)

        return inds, shapes, coeffs

    def reconstruct(self, coeffs):
        """Inverse 2D1D wavelet transform.

        Parameters
        ----------
        coeffs : 1D array
            Wavelet coefficients plus index markers packaged as a 1D array,
            i.e. the output of decompose().

        Returns
        -------
        cube : 3D array
            Reconstructed 3D data cube.

        """
        assert hasattr(self, "_mr2d1d"), "Need to call decompose first."

        return self._mr2d1d.reconstruct(coeffs)

    def energy_per_scale(self, num_scales_2d, num_scales_1d):
        return self.NOISE_TAB[:num_scales_2d, :num_scales_1d]

    @property
    def num_precomputed(self):
        return self.NOISE_TAB.shape

    def _extract_metadata(self, coeffs):
        """Get metadata of transformed coefficients for all pairs of scales.

        Parameters
        ----------
        coeffs : 1D array
            Wavelet coefficients plus index markers packaged as a 1D array,
            i.e. the output of decompose().

        Returns
        -------
        inds : nested list of tuples
            Pairs of index values arranged such that inds[i][j] gives the
            (start, end) indices for the band (2d_scale, 1d_scale) = (i, j).
            The actual coefficients of the (i, j) band can therefore be accessed
            as coeffs[start:end].
        shapes : nested list of tuples
            Triplets of sub-band cube shapes arranged such that shapes[i][j]
            gives the (nx, ny, nz) shape of the band (2d_scale, 1d_scale) = (i, j).

        """
        n_scales_2d = int(coeffs[0])
        n_scales_1d = int(coeffs[1])

        inds = [[() for _ in range(n_scales_1d)] for _ in range(n_scales_2d)]
        shapes = [[() for _ in range(n_scales_1d)] for _ in range(n_scales_2d)]

        # Starting index
        start = end = 2

        # Traverse the flattened array to pull out ranges for each index pair
        for ii in range(n_scales_2d):
            for jj in range(n_scales_1d):
                # Starting index for this band
                start = end + 3
                # Extract band sizes
                nx, ny, nz = map(int, coeffs[start - 3 : start])
                shapes[ii][jj] = (nx, ny, nz)
                # Total number of coefficients in this band
                ncoeff = nx * ny * nz
                # Ending index for this band
                end = start + ncoeff
                inds[ii][jj] = (start, end)

        return inds, shapes


# ## 3D Wavelet denoising
#
# The Denoiser2D1D class below is modeled after Aymeric's StarletDenoiser that does 2D denoising on each frequency slice of the input data. The overall structure is kept so that it will be easier to combine the two denoising options (2D or 3D) later for the user's convenience.
#
# NOTE : the denoise() method here accepts an optional noise cube realisation as input. The noise is transformed and its standard deviation is computed in each wavelet band to establish the noise level there. If not provided, the noise in each band is estimated from the data.

# In[3]:


class Denoiser2D1D(object):
    """Denoise a data cube using a 2D1D wavelet decomposition.

    The model is simply Y = X + N, where X is the noiseless signal and N is
    the noise.

    """

    def __init__(self, threshold_type="soft", verbose=True):
        """Initialise the denoiser.

        Parameters
        ----------
        threshold_type : string
            Supported values are 'soft' or 'hard'. Default: 'soft'
        verbose : bool
            If True, prints some more info, useful for debugging. Default: True

        """
        self.mr2d1d = Wavelet2D1DTransform()
        self._threshold_type = threshold_type
        self._verbose = verbose

    def __call__(self, *args, **kwargs):
        """Alias for self.denoise()"""
        return self.denoise(*args, **kwargs)

    def denoise(
        self,
        x,
        method="simple",
        threshold_level=3,
        threshold_increment_high_freq=2,
        num_scales_2d=None,
        num_scales_1d=None,
        noise_cube=None,
        **kwargs_method,
    ):
        """Denoise a data cube according to the chosen method.

        Parameters
        ----------
        x : array_like (3D)
            Input data cube. The frequency axis is assumed to be first.
        method : string
            Denoising method, either 'simple' or 'iterative'. The iterative method
            should give better results but takes longer to compute.
            Default: 'simple'.
        threshold_level : int
            Threshold level, as a detection signicance, in noise units (generally
            between 3 and 5 for '3-5 sigmas' detection).
            Default: 3
        threshold_increment_high_freq : int
            Increment of the above threshold_level for the highest frequencies (usually
            associated with pure noise).
            Default: 2
        num_scales_2d : int
            Number of starlet decomposition scales for the 2D images. Maximal value is
            int(np.log2(input_image_.shape[-1])).
            Default: None (max. value).
        num_scales_1d : int
            Number of wavelet scales for the 1D axis. Maximal value is
            int(np.log2(input_image_.shape[0])).
            Default: None (max. value).
        noise_cube : array_like, same shape as input `x`
            An estimate of the noise (e.g. by simulation). If not provided, the noise level
            is estimated automatically in each wavelet sub-band.
            Default: None

        kwargs_method : dict
            [See docstring of each method]

        Returns
        -------
        array_like
            Denoised array

        """
        # Set the number of 2D decomposition scales
        num_scales_2d_max = int(np.log2(x.shape[1]))
        if (
            num_scales_2d is None
            or num_scales_2d < 2
            or num_scales_2d > num_scales_2d_max
        ):
            # choose the maximum allowed number of scales
            num_scales_2d = num_scales_2d_max
            if self._verbose is True:
                print(
                    f"Number of 2D wavelet scales set to {num_scales_2d} "
                    "(maximum value allowed by input image)"
                )

        # Set the number of 1D decomposition scales
        num_scales_1d_max = int(np.log2(x.shape[0]))
        if (
            num_scales_1d is None
            or num_scales_1d < 2
            or num_scales_1d > num_scales_1d_max
        ):
            # choose the maximum allowed number of scales
            num_scales_1d = num_scales_1d_max
            if self._verbose is True:
                print(
                    f"Number of 1D wavelet scales set to {num_scales_1d} "
                    "(maximum value allowed by input image)"
                )

        # Check that the pre-computed noise scaling exists for the requested scales
        # if (num_scales_2d - 1 > self.mr2d1d.num_precomputed[0] or
        #     num_scales_1d - 1 > self.mr2d1d.num_precomputed[1]):
        #     raise NotImplementedError(f"Pre-computed noise in wavelet space has been implemented"
        #                               f" for up to {self.mr2d1d.NOISE_TAB.shape} scales "
        #                               f"[({num_scales_2d}, {num_scales_1d}) required)]")

        # Check that the noise realisation has the same shape as the input
        if noise_cube is not None:
            assert x.shape == noise_cube.shape, "Invalid noise estimate shape"

        # Initialise settings for the denoiser
        self._data = x
        self._num_bands = self._data.shape[0]
        self._num_pixels = self._data.shape[1] * self._data.shape[2]
        self._num_scales_2d = num_scales_2d
        self._num_scales_1d = num_scales_1d
        self._thresh_min = float(threshold_level)
        self._thresh_increm = float(threshold_increment_high_freq)
        self._noise = noise_cube

        # Select and run the denoiser
        if method == "simple":
            result = self._denoise_simple()
        elif method == "iterative":
            print("ITERATIVE DENOISING IS A WORK IN PROGRESS")
            result = self._denoise_iterative(**kwargs_method)
        else:
            raise ValueError(f"Denoising method '{method}' is not supported")

        return result

    def _denoise_simple(self):
        """Denoise the data using a one-step thresholding in 3D wavelet space.

        This is effectively equivalent to self._denoise_iterative() with num_iter=1
        and progressive_threshold=False. A positivity constraint is enforced.

        Returns
        -------
        array_like
            Denoised array

        """
        # Denoise considering the full 3D wavelet deconstruction
        inds, shapes, w_data = self.mr2d1d.decompose(
            self._data, self._num_scales_2d, self._num_scales_1d
        )

        # Forward transform the noise cube if provided
        if self._noise is not None:
            _, _, w_noise = self.mr2d1d.decompose(
                self._noise, self._num_scales_2d, self._num_scales_1d
            )

        # Threshold each wavelet sub-band
        for scale2d in range(self._num_scales_2d - 1):
            for scale1d in range(self._num_scales_1d - 1):
                # Extract the coefficients for this sub-band
                start, end = inds[scale2d][scale1d]
                c_data = w_data[start:end]

                # Noise standard deviation of this sub-band
                if self._noise is not None:
                    noise_level = w_noise[start:end].std()
                else:
                    # Estimate the noise std. dev. using the MAD
                    noise_level = self._estimate_noise(c_data)
                    print(scale2d, " ", scale1d, ": noise = ", noise_level)

                # Filter the data with one pass
                # TODO : COULD APPLY THE THRESHOLD INCREMENT HERE SOMEHOW
                thresh = self._thresh_min
                # if scale2d == 0:
                #     thresh += self._thresh_increm

                # Update band with thresholded coefficients
                w_data[start:end] = self._prox_sparsity_constraint(
                    c_data, thresh, noise_level
                )

        # Bring back to direct space by inverse transform
        result = self.mr2d1d.reconstruct(w_data)

        # Apply the positivity constraint
        result = self._prox_positivity_constraint(result)

        return result

    def _denoise_iterative(
        self, num_iter=20, num_reweight=2, progressive_threshold=True
    ):
        """Denoise the data using an iterative thresholding in wavelet space.

        Parameters
        ----------
        num_iter : int
            Number of iterations. Default: 20
        num_reweight : int
            Number of times l1-reweighting is applied. Default: 2
        progressive_threshold : bool
            If True, the threshold is exponentially decreased from an initial values
            estimated from the data, to the chosen final value. Default: True

        Returns
        -------
        array_like
            Denoised array

        """
        if self._threshold_type == "hard":
            num_reweight = 1

        num_iter_min_threshold = 5

        # Initialise the model
        model = np.zeros_like(self._data)

        # Gradient of the loss function
        grad = lambda x: self.gradient(x, self._data)

        # Gradient descent step size
        step_size = 1.0

        # No weights for the first pass
        weights = None

        for _ in range(num_reweight):
            # Forward transform the data
            inds, shapes, w_data = self.mr2d1d.decompose(
                self._data, self._num_scales_2d, self._num_scales_1d
            )

            # Forward transform the noise cube if provided
            if self._noise is not None:
                _, _, w_noise = self.mr2d1d.decompose(
                    self._noise, self._num_scales_2d, self._num_scales_1d
                )

            # Get initial threshold value (in noise units)
            if progressive_threshold and num_iter > num_iter_min_threshold:
                thresh = 0.0
                # Loop through all sub-bands to find the max usable threshold
                for scale2d in range(self._num_scales_2d - 1):
                    for scale1d in range(self._num_scales_1d - 1):
                        start, end = inds[scale2d][scale1d]
                        c_data = w_data[start:end]
                        if self._noise is not None:
                            noise_level = w_noise[start:end].std()
                        else:
                            noise_level = self._estimate_noise(c_data)
                        print("BUG: ", c_data.shape)
                        thresh = np.max(thresh, 0.9 * c_data / noise_level)
            else:
                thresh = self._thresh_min

            # for _ in range(num_iter):
            #    def prox(x, y):
            #
            #        return
            #    prox = lambda x, y: self._prox_sparsity_constraint(x, thresh, noise_level)
            #
            #    model_next =

        return model

    @staticmethod
    def step_forward_backward(x, grad, prox, step_size):
        """One step of the Forward Backward Algorithm.

        Parameters
        ----------
        x : array_like
            Variable being optimised
        grad : function
            Function that takes as a unique argument the variable x
        prox : function
            Function that takes two positional arguments: the variable being optimised, and the gradient step size
        step_size : float
            Step size for the gradient descent step

        Returns
        -------
        array_like
            Updated variable
        """
        x_next = prox(x - step_size * grad(x), step_size)
        return x_next

    @staticmethod
    def gradient(model, data):
        """Gradient of the data-fidelity term in the lost function, with respect to the main variable.
        The model is Y = X + N, where X in the denoised signal and N is the noise.
        The data-fidelity term is the mean-squared error ||Y' - Y||, where Y' represents the data.

        Parameters
        ----------
        model : array_like
            Any array
        data : array_like
            Any array

        Returns
        -------
        array_like
            Gradient of the data-fidelity term

        """
        res = data - model
        grad = -res
        return grad

    def _prox_sparsity_constraint(self, coeffs, thresh, noise_level, weights=None):
        """Proximal operator of the l0- or l1-sparsity constraint.

        Parameters
        ----------
        coeffs : array_like
            1D array of this sub-band
        thresh : float
            Threshold value in noise units
        noise_level : float
            Noise level at the scale represented by `array`
        weights : array_like
            Weights per pixel per scale. If None, no weights are applied

        Returns
        -------
        array_like
            Array on which the sparsity constraint has been applied

        """
        if weights is None:
            weights = np.ones_like(coeffs)

        # Threshold coefficients
        lmbda = thresh * noise_level * weights
        return self.threshold(coeffs, lmbda, self._threshold_type)

    @staticmethod
    def _prox_positivity_constraint(array):
        """Proximal operator of the positivity constraint

        Parameters
        ----------
        array : array_like
            Any array that supports index slicing

        Returns
        -------
        array_like
            Array with all negative entries set to zero

        """
        return np.maximum(0, array)

    @staticmethod
    def threshold(array, threshold_value, threshold_type="soft"):
        """Translate the noise estimation from direct space to starlet space

        Parameters
        ----------
        array : array_like
            1D or 2D array
        threshold_value : float
            Threshold in same units as array
        threshold_type : string
            Supported values are 'soft' or 'hard'

        Returns
        -------
        array_like
            Thresholded array

        Raises
        ------
        ValueError
            If the input array is not 1D or 2D
        ValueError
            If threshold_type is not supported
        """
        if len(array.shape) > 2:
            raise ValueError(f"Soft thresholding only supported for 1D or 2D arrays")
        if threshold_type == "soft":
            array_th = np.sign(array) * np.maximum(np.abs(array) - threshold_value, 0.0)
        elif threshold_type == "hard":
            array_th = np.copy(array)
            array_th[np.abs(array) <= threshold_value] = 0.0
        else:
            raise ValueError(f"Threshold type '{threshold_type}' is not supported")
        return array_th

    #     def _propagate_noise(self, noise, scale_2d, scale_1d):
    #         """Translate the noise estimation from direct space to starlet space

    #         Parameters
    #         ----------
    #         noise : float
    #             Noise value in direct space

    #         Returns
    #         -------
    #         float
    #             Noise per starlet scale

    #         """
    #         # Scale the noise to each wavelet sub-band
    #         norm = self.mr2d1d.energy_per_scale(self._num_scales_2d, self._num_scales_1d)
    #         return noise * norm[scale_2d, scale_1d]

    def _compute_weights(self, array, noise_levels, threshold):
        """Compute weight (per pixel per starlet scale) for l1-reweighting scheme
        (used for 'soft' threshold_type) based on the starlet decomposition of the input array

        Parameters
        ----------
        array : array_like
            Image on which the noise is estimated
        noise_levels : array_like
            1D array or list containing the noise level per starlet scale
        threshold : float
            Threshold in noise units

        Returns
        -------
        float
            Cube containing a weight value for each pixel and each scale
        """
        coeffs = self.starlet.decompose(array, self._num_scales)
        # construct cube with with constant threshold level per decomposition scale
        thresh_ = np.stack([nl * np.ones_like(array) for nl in noise_levels])
        thresh_[0, :, :] *= threshold + self._thresh_increm
        thresh_[1:, :, :] *= threshold
        # compute weights
        weights = 1.0 / (1 + np.exp(10 * (coeffs - thresh_)))
        return weights

    def _estimate_noise(self, array):
        """Estimate noise standard deviation from the median absolute deviation (MAD)

        Parameters
        ----------
        array : array_like
            Values on which the noise is estimated

        Returns
        -------
        float
            Noise standard deviation

        """
        mad = np.median(np.abs(array - np.median(array)))
        return 1.48 * mad

    def _estimate_threshold(self, array, noise_levels, fraction=0.9):
        """
        estimate maximum threshold, in units of noise, used for thresholding wavelets
        coefficients during optimization

        Parameters
        ----------
        data : array_like
            Imaging data.
        fraction : float, optional
            From 0 to 1, fraction of the maximum value of the image in transformed space,
            normalized by noise, that is returned as a threshold.

        Returns
        -------
        float
            Threshold level.
        """
        coeffs = self.starlet.decompose(array, self._num_scales)
        coeffs_norm = coeffs / noise_levels[:, None, None]
        coeffs_norm[noise_levels == 0] = 0
        # returns a fraction of max value, so only the highest coeffs is able to enter the solution
        return fraction * np.max(coeffs_norm[:-1])  # ignore last scale

    @staticmethod
    def exponential_decrease(
        curr_value, init_value, min_value, num_iter, num_iter_at_min_value
    ):
        """Computes a exponentially decreasing value, for a given loop index, starting at a specified value.

        Parameters
        ----------
        curr_value : float
            Current value to be updated
        init_value : float
            Value at iteration 0.
        min_value : float
            Minimum value, reached at iteration num_iter - num_iter_at_min_value - 1.
        num_iter : int
            Total number of iterations.
        num_iter_at_min_value : int
            Number of iteration for which the returned value equals `min_value`.

        Returns
        -------
        float
            Exponentially decreased value.

        Raises
        ------
        ValueError
            If num_iter - num_iter_at_min_value < 1, cannot compute the value.
        """
        num_iter_eff = num_iter - num_iter_at_min_value
        if num_iter_eff < 1:
            raise ValueError(
                f"Too low number of iterations ({num_iter}) to decrease threshold"
            )
        exp_factor = np.exp(np.log(min_value / init_value) / num_iter_eff)
        new_value = curr_value * exp_factor
        return max(new_value, min_value)


# ## Tests on mock data

# In[5]:


def run_test1():
    # Cube size
    nx = ny = 64
    nz = 80

    dirac2d = np.zeros((nx, ny))
    dirac2d[nx // 2, ny // 2] = 100.0
    dirac2d[7, 7] = 120.0
    dirac2d[-10, -20] = 100.0
    dirac2d[22, -22] = 80.0

    # Gaussian blobs shrinking along the z axis
    signal3d = [gaussian_filter(dirac2d, sigma) for sigma in np.linspace(3, 1, nz)]
    signal3d = np.array(signal3d)

    # Set noise sigma to be 20% of the max. signal in each z slice
    sigmas = np.max(signal3d, axis=(1, 2)) / 5
    noise3d = np.array([np.random.normal(0, sigma, (nx, ny)) for sigma in sigmas])
    data3d = signal3d + noise3d

    # Plot a z-slice
    iz = 10
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    ax1.imshow(signal3d[iz])
    ax1.set_title(f"Clean z-slice ({iz})", fontsize=16)
    ax2.imshow(data3d[iz])
    ax2.set_title(f"Noisy z-slice ({iz})", fontsize=16)

    # In[6]:

    # Create the denoiser  hard or soft
    denoiser3d = Denoiser2D1D(threshold_type="hard")
    # Denoise    simple or iterative
    denoised3d = denoiser3d(
        data3d,
        threshold_level=3,
        method="simple",
        num_scales_2d=None,
        num_scales_1d=None,
        noise_cube=None,
    )

    # In[7]:

    # Plot a z-slice
    iz = 10
    style = {"vmin": data3d[iz].min(), "vmax": data3d[iz].max()}

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    (ax1, ax2), (ax3, ax4) = axes
    img1 = ax1.imshow(data3d[iz])
    ax1.set_title("Data")
    img2 = ax2.imshow(denoised3d[iz], **style)
    ax2.set_title("Denoised 3D")
    img3 = ax3.imshow(
        (data3d[iz] - denoised3d[iz]) / noise3d[iz].std(), vmin=-5, vmax=5
    )
    ax3.set_title("(Normalised) Residual")
    img4 = ax4.imshow(signal3d[iz] - denoised3d[iz], **style)
    ax4.set_title("Absolute Error")

    # Add colorbars
    imgs = (img1, img2, img3, img4)
    for ax, img in zip(np.ravel(axes), imgs):
        fig.colorbar(img, ax=ax)

    # In[8]:

    # Plot an interesting y-slice
    iy = ny // 2 + 11
    style = {"vmin": data3d[:, :, iy].min(), "vmax": data3d[:, :, iy].max()}

    fig, axes = plt.subplots(2, 2, figsize=(10, 10))
    (ax1, ax2), (ax3, ax4) = axes
    img1 = ax1.imshow(data3d[:, :, iy])
    ax1.set_title("Data")
    img2 = ax2.imshow(denoised3d[:, :, iy], **style)
    ax2.set_title("Denoised 3D")
    img3 = ax3.imshow((data3d[:, :, iy] - denoised3d[:, :, iy]), vmin=-5, vmax=5)
    ax3.set_title("(Absolute) Residual")
    img4 = ax4.imshow(signal3d[:, :, iy] - denoised3d[:, :, iy], **style)
    ax4.set_title("Absolute Error")

    # Add colorbars
    imgs = (img1, img2, img3, img4)
    for ax, img in zip(np.ravel(axes), imgs):
        fig.colorbar(img, ax=ax)

    # In[ ]:


nx = ny = 64
nz = 80
d = readfits("test.fits")
d1 = d[0:nz, 0:nx, 0:ny]
d1 = np.random.normal(0, 10, (nz, nx, ny))

CDen = Denoiser2D1D(threshold_type="hard")
r = CDen(
    d,
    threshold_level=3,
    method="simple",
    num_scales_2d=5,
    num_scales_1d=None,
    noise_cube=None,
)
writefits("result.fits", r)

r5 = CDen(
    d,
    threshold_level=5,
    method="simple",
    num_scales_2d=5,
    num_scales_1d=None,
    noise_cube=None,
)
writefits("result5.fits", r5)
