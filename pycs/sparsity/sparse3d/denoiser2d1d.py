import matplotlib.pyplot as plt
import numpy as np

from pycs.sparsity.sparse3d.wavelet2d1d_transform import Wavelet2D1DTransform


class Denoiser2D1D(object):
    """
    Advanced 3D Spectral Cube Denoiser using 2D-1D Multi-Scale Wavelets

    This class implements sophisticated denoising algorithms for IFU spectral cubes
    using a hybrid 2D-1D wavelet decomposition. The approach combines spatial and
    spectral information to achieve superior noise reduction while preserving both
    morphological and kinematic features.

    The denoising framework assumes an additive noise model:
        Y = X + N
    where Y is the observed noisy cube, X is the true signal, and N is the noise.

    Available Denoising Methods:
    ---------------------------
    1. Iterative Hard: Multiple iterations with binary (hard) thresholding
    2. Iterative Soft: Multiple iterations and re-weighting with adaptive (soft) thresholding

    The iterative methods use advanced techniques:
    - Adaptive re-weighting to reduce bias from soft thresholding
    - Residual signal extraction to recover previously missed features
    - Plateau-based convergence criteria for robust stopping

    Key Advantages:
    --------------
    - Preserves spatial morphology through 2D starlet decomposition
    - Maintains spectral line profiles via 1D wavelet analysis
    - Adapts thresholds to local noise characteristics
    - Recovers faint emission through iterative refinement
    - Supports both synthetic and observational noise modeling

    Typical Applications:
    --------------------
    - ALMA/JWST IFU observations with low SNR
    - High-redshift galaxy emission line recovery
    - Continuum-subtracted line cube cleaning
    - Kinematic analysis preprocessing
    - Extended emission detection enhancement

    Attributes
    ----------
    mr2d1d : Wavelet2D1DTransform
        The wavelet transform object for decomposition/reconstruction
    _threshold_type : str
        Thresholding strategy ('soft' or 'hard')
    _verbose : bool
        Controls progress reporting and diagnostics
    _plot : bool
        Enables diagnostic plotting during processing

    Methods
    -------
    denoise(x, y, method='simple', **kwargs)
        Main denoising interface with multiple algorithm options
    __call__(*args, **kwargs)
        Convenience alias for denoise() method

    Private Methods
    ---------------

    _denoise_iterative_hard(**kwargs)
        Multi-iteration hard thresholding with L0 regularization
    _denoise_iterative_soft(**kwargs)
        Multi-iteration soft thresholding with adaptive reweighting
    _generate_hard_threshold_mask(coeffs, thresh, noise_level)
        Create binary masks for hard thresholding
    _residual_signal_extraction_l0/l1(...)
        Extract additional signal from residuals
    _apply_positivity(arr)
        Apply positivity constraint when enabled; return arr unchanged otherwise
    _estimate_noise(array)
        Robust noise estimation using median absolute deviation
    _compute_emission_rmse(model)
        Calculate reconstruction error in emission regions
    """

    def __init__(self, threshold_type='soft', verbose=True, plot=False):
        """
        Initialize the 3D denoising framework with specified parameters.

        Sets up the wavelet transform object and configures algorithm parameters
        for subsequent denoising operations. The choice of threshold type determines
        which iterative algorithms are available.

        Parameters
        ----------
        threshold_type : str, optional
            Thresholding strategy for coefficient processing, by default 'soft'.
            Options:
            - 'soft': Shrink coefficients by threshold amount (preserves gradients)
            - 'hard': Set coefficients below threshold to zero (creates sparsity)
            Soft thresholding generally produces smoother results but may over-smooth.
            Hard thresholding preserves sharp features but can introduce artifacts.
        verbose : bool, optional
            Enable detailed progress reporting and diagnostics, by default True.
            When True, prints iteration progress, convergence status, flux measurements,
            and algorithm-specific information during processing.
        plot : bool, optional
            Enable diagnostic plotting during denoising, by default False.
            When True, displays:
            - Coefficient distributions before/after thresholding
            - Iteration-by-iteration reconstruction progress
            - Residual analysis and signal extraction visualization
            - Final comparison plots (original/denoised/residual)
            Useful for algorithm development and result validation.

        Notes
        -----
        The wavelet transform object is created with default parameters and can
        handle cubes with dimensions up to the limits of available memory.
        Transform parameters are automatically configured based on input cube
        dimensions during the first call to denoise().

        For large cubes (>1GB), consider setting verbose=False to reduce I/O
        overhead during processing.
        """
        # Initialize the 2D-1D wavelet transform engine
        self.mr2d1d = Wavelet2D1DTransform()

        # Store algorithm configuration
        self._threshold_type = threshold_type
        self._verbose = verbose
        self._plot = plot

    def __call__(self, *args, **kwargs):
        """
        Convenience method to call denoise() directly on the object.

        Allows using the denoiser object as a function:
        denoiser(cube, signal, method='iterative')

        This is equivalent to denoiser.denoise(cube, signal, method='iterative')

        Returns
        -------
        Denoised cube or tuple of results depending on method used.
        """
        return self.denoise(*args, **kwargs)

    def denoise(self, x, y, method='iterative', threshold_level=3,
                threshold_increment_high_freq=2, num_scales_2d=None,
                num_scales_1d=None, noise_cube=None, emission_mask=None,
                positivity=False, positivity_final=False, **kwargs_method):
        """
        Denoise a 3D spectral cube using advanced 2D-1D wavelet techniques.

        This is the main interface for all denoising algorithms. The method applies
        multi-scale wavelet decomposition followed by intelligent thresholding to
        remove noise while preserving both spatial morphology and spectral features.

        The algorithm automatically determines optimal scale numbers if not specified,
        adapts thresholds to noise characteristics, and can incorporate ground truth
        for algorithm validation and noise modeling.

        Parameters
        ----------
        x : np.ndarray, shape (nz, ny, nx)
            Input noisy data cube. The spectral/frequency axis must be first (axis=0),
            followed by spatial dimensions. Each x[i,:,:] is a 2D image at wavelength i.
        y : np.ndarray, shape (nz, ny, nx)
            Clean signal cube (ground truth) with same shape as x.
            Used for algorithm validation, noise modeling, and performance assessment.
            In real observations, this would be unknown - here used for synthetic testing.
        method : str, optional
            Denoising algorithm to use, by default 'simple'.
            Options:
            - 'simple': Single-iteration thresholding (fast, good for high SNR)
            - 'iterative': Multi-iteration adaptive algorithm (slower, better for low SNR)
              For iterative method, the specific algorithm depends on threshold_type:
              * 'hard' → Iterative hard thresholding with L0 regularization
              * 'soft' → Iterative soft thresholding with adaptive reweighting
        threshold_level : float, optional
            Base threshold level in noise standard deviations, by default 3.
            Typical range: 2-5 for 2-5σ detection significance.
            Lower values preserve more signal but retain more noise.
            Higher values create cleaner results but may remove faint features.
        threshold_increment_high_freq : float, optional
            Additional threshold increment for highest frequency scales, by default 2.
            These scales typically contain pure noise and benefit from higher thresholds.
            Final threshold = threshold_level + threshold_increment_high_freq for finest scales.
        num_scales_2d : int, optional
            Number of 2D spatial starlet scales, by default None (auto-determined).
            If None, uses maximum: int(log2(min(ny, nx))).
            Range: 2 to maximum allowed by image dimensions.
            More scales capture finer spatial details but increase computation.
        num_scales_1d : int, optional
            Number of 1D spectral wavelet scales, by default None (auto-determined).
            If None, uses maximum: int(log2(nz)).
            Range: 2 to maximum allowed by spectral dimensions.
            More scales capture narrower spectral features but increase computation.
        noise_cube : np.ndarray, shape (nz, ny, nx), optional
            Independent noise realization with same shape as input, by default None.
            If provided, used to accurately estimate noise levels in each wavelet sub-band.
            If None, noise levels estimated from the data using robust statistics.
            Recommended for synthetic data where true noise is known.
        emission_mask : np.ndarray, shape (nz, ny, nx), optional
            Binary mask indicating emission regions (1=emission, 0=background), by default None.
            If None, creates mask of all ones (assumes entire cube contains signal).
            Used for targeted error calculation and emission-focused denoising.
            Particularly useful for line emission cubes with known source extent.
        positivity : bool, optional
            Whether to enforce non-negativity on the reconstructed model, by default True.
            Set to True for emission-only sources where negative flux is unphysical.
            Set to False when the data may contain absorption features or is mean-subtracted.
        **kwargs_method : dict
            Additional keyword arguments passed to specific denoising methods:

            For method='iterative' with threshold_type='hard':
            - num_iter : int, default 20
                Number of hard thresholding iterations

            For method='iterative' with threshold_type='soft':
            - num_iter_reweight : int, default 20
                Number of reweighting iterations
            - num_iter_debias : int, default 20
                Number of debiasing iterations
            - debias : bool, default True
                Whether to perform debiasing step

        Returns
        -------
        result : np.ndarray or tuple
            Denoised results. Return format depends on method:

            For method='simple':
                np.ndarray, shape (nz, ny, nx) : Denoised cube

            For method='iterative':
                tuple with multiple outputs depending on threshold_type:

                Hard thresholding returns:
                - best_model : Denoised cube at best iteration
                - deltas : Accumulated residual signals extracted
                - residual_stds : Standard deviation history per iteration
                - best_iteration : Iteration number of best result
                - noise_levels : Noise estimates per wavelet sub-band

                Soft thresholding returns:
                - best_model : Final denoised cube
                - model_1_step : Result after first reweighting phase
                - model_no_reweight : Result without reweighting
                - deltas : Accumulated residual signals from debiasing
                - residual_stds_reweight : Residual history during reweighting
                - residual_stds_debias : Residual history during debiasing
                - best_iteration : Iteration of best result
                - dists : Selected coefficient distributions for diagnostics
                - noise_levels : Noise level estimates per sub-band

        Raises
        ------
        ValueError
            If method is not 'simple' or 'iterative'
        AssertionError
            If noise_cube provided but shape doesn't match input x

        Notes
        -----
        Processing Time:
        - Simple method: ~10-30 seconds for 100³ cube
        - Iterative method: ~2-10 minutes depending on iterations and convergence

        Memory Requirements:
        - Peak usage ~4-5x input cube size during transform operations
        - Coefficient storage ~1.5x input cube size

        Algorithm Selection Guide:
        - High SNR (>10): Simple method often sufficient
        - Low SNR (<5): Iterative soft thresholding recommended
        - Preserving sharp features: Hard thresholding
        - Smooth reconstruction: Soft thresholding
        """

        # Validate and set default number of 2D decomposition scales
        # Use the smaller spatial dimension so non-square cubes don't exceed wavelet limits
        num_scales_2d_max = int(np.log2(min(x.shape[1], x.shape[2])))
        if num_scales_2d is None or num_scales_2d < 2 or num_scales_2d > num_scales_2d_max:
            num_scales_2d = num_scales_2d_max
            if self._verbose is True:
                print(f"Number of 2D wavelet scales set to {num_scales_2d} "
                      "(maximum value allowed by input image)")

        # Set the number of 1D decomposition scales
        num_scales_1d_max = int(np.log2(x.shape[0]))
        if num_scales_1d is None or num_scales_1d < 2 or num_scales_1d > num_scales_1d_max:
            num_scales_1d = num_scales_1d_max
            if self._verbose is True:
                print(f"Number of 1D wavelet scales set to {num_scales_1d} "
                      "(maximum value allowed by input image)")

        if noise_cube is not None:
            assert x.shape == noise_cube.shape, "Invalid noise estimate shape"

        # Initialise settings for the denoiser
        self._data = x
        self._signal = y
        self._num_bands = self._data.shape[0]
        self._num_pixels = self._data.shape[1] * self._data.shape[2]
        self._num_scales_2d = num_scales_2d
        self._num_scales_1d = num_scales_1d
        self._threshold_level = float(threshold_level)
        self._thresh_increm = float(threshold_increment_high_freq)
        self._noise = noise_cube
        self._positivity = positivity
        if emission_mask is None:
            emission_mask = np.ones_like(y)
        self._mask = emission_mask

        # Select and run the denoiser
        if method == 'simple':
            if self._verbose:   print('\n--- [ PERFORMING SIMPLE (ONE-STEP) DENOISING ] ---\n')
            result = self._denoise_simple()
        elif method == 'iterative':
            if self._verbose:   print('\n--- [ PERFORMING ITERATIVE DENOISING ] ---\n')
            if self._threshold_type == 'hard':
                result = self._denoise_iterative_hard(**kwargs_method)
            if self._threshold_type == 'soft':
                result = self._denoise_iterative_soft(**kwargs_method)
        else:
            raise ValueError(f"Denoising method '{method}' is not supported")

        if positivity_final:
            n_clipped = int(np.sum(result[0] < 0))
            result = (np.maximum(0, result[0]),) + result[1:]
            if self._verbose:
                print(f'(*) Positivity applied to final model: {n_clipped} negative voxels clipped to 0')

        if self._verbose:
            print(f'(*) Final aperture flux: {np.sum(result[0]):.3e}')

        return result

    def _apply_positivity(self, arr):
        """Apply positivity constraint when enabled; return arr unchanged otherwise.

        Parameters
        ----------
        arr : np.ndarray
            Input array (model or delta) to constrain.

        Returns
        -------
        np.ndarray
            Element-wise max(0, arr) when ``self._positivity`` is True,
            otherwise arr unchanged.
        """
        return np.maximum(0, arr) if self._positivity else arr

    def _decompose_and_estimate_noise(self):
        """Decompose data and compute per-sub-band noise standard deviations.

        Uses the independent noise cube (``self._noise``) when provided — the cube
        is transformed and per-band std is taken as the exact noise level, avoiding
        any bias from signal contamination. Falls back to MAD on the data's wavelet
        coefficients when no noise cube is available.

        Returns
        -------
        inds : list of list of tuples
            Sub-band index structure as returned by ``Wavelet2D1DTransform.decompose``.
        shapes : list of list of tuples
            Sub-band shape structure.
        w_data : np.ndarray, 1D
            Full wavelet coefficient array of the noisy data (read-only reference;
            callers must copy before modifying).
        noise_levels : list of float
            Noise standard deviation estimate for each non-coarse sub-band, in
            traversal order (scale2d outer, scale1d inner).
        """
        inds, shapes, w_data_buf = self.mr2d1d.decompose(
            self._data, self._num_scales_2d, self._num_scales_1d
        )
        # pysparse returns a reference to its internal C buffer; copy before any
        # subsequent decompose call can overwrite it.
        w_data = w_data_buf.copy()

        if self._noise is not None:
            # Transform an independent noise realization; std per band is unbiased
            _, _, w_noise = self.mr2d1d.decompose(
                self._noise, self._num_scales_2d, self._num_scales_1d
            )
            def _band_noise(start, end):
                return float(np.std(w_noise[start:end]))
        else:
            # Robust MAD fallback: works for non-white or non-Gaussian noise
            def _band_noise(start, end):
                return self._estimate_noise(w_data[start:end])

        noise_levels = []
        for scale2d in range(self._num_scales_2d):
            for scale1d in range(self._num_scales_1d):
                if scale2d == self._num_scales_2d - 1 and scale1d == self._num_scales_1d - 1:
                    continue  # coarse approximation band excluded
                start, end = inds[scale2d][scale1d]
                noise_levels.append(_band_noise(start, end))

        return inds, shapes, w_data, noise_levels

    def _denoise_simple(self):
        """
        Perform single-pass hard thresholding denoising (one IHT iteration).

        This is a fast, non-iterative denoiser equivalent to the first iteration of
        `_denoise_iterative_hard`. It decomposes the data into wavelet coefficients,
        estimates noise per sub-band via MAD, applies a binary hard threshold mask to
        discard insignificant coefficients, and reconstructs the denoised cube.

        The finest spatial scale receives a higher threshold
        (``threshold_level + threshold_increment_high_freq``) because it typically
        contains the largest proportion of pure noise.

        Returns
        -------
        model : np.ndarray, shape (nz, ny, nx)
            Denoised cube after a single thresholding pass with positivity applied.
        noise_levels : list of float
            Noise standard deviation estimates for each non-coarse wavelet sub-band,
            in traversal order (scale2d outer, scale1d inner).

        See Also
        --------
        _denoise_iterative_hard : Multi-iteration extension of this approach.
        _generate_hard_threshold_mask : Binary mask construction.
        """
        if self._verbose:
            print('(*) Decomposing data and estimating per-band noise levels')

        inds, shapes, w_data, noise_levels = self._decompose_and_estimate_noise()
        if self._verbose:
            src = 'noise cube (std)' if self._noise is not None else 'data (MAD)'
            print(f'    Noise estimated from {src} across {len(noise_levels)} sub-bands')
            print(f'(*) Applying hard thresholding (lambda = {self._threshold_level} sigma)')

        mask_coeff = np.zeros_like(w_data)
        i = 0
        for scale2d in range(self._num_scales_2d):
            for scale1d in range(self._num_scales_1d):
                start, end = inds[scale2d][scale1d]
                if scale2d == self._num_scales_2d - 1 and scale1d == self._num_scales_1d - 1:
                    continue  # leave coarse approximation band untouched
                thresh = self._threshold_level + (self._thresh_increm if scale2d == 0 else 0)
                mask_coeff[start:end] = self._generate_hard_threshold_mask(
                    w_data[start:end], thresh, noise_levels[i]
                )
                i += 1

        # Save coarse band before masking then restore it after — the mask loop leaves
        # it at 0.0 (via continue), which would make reconstruct() silently use
        # pysparse's internally cached coarse plane and give model ≈ data.
        start_c, end_c = inds[self._num_scales_2d - 1][self._num_scales_1d - 1]
        coarse_save = w_data[start_c:end_c].copy()

        w_data *= mask_coeff
        w_data[start_c:end_c] = coarse_save

        if self._verbose:
            print('(*) Reconstructing the denoised data from wavelet to real space')

        model = self._apply_positivity(self.mr2d1d.reconstruct(w_data))

        if self._verbose:
            residual_std = np.std(self._data - model)
            print(f'(*) Done. Aperture flux: {np.sum(model):.3e}, Residual STD: {residual_std:.3e}')

        return model, noise_levels

    def _residual_signal_extraction_l0(self, model, mask_coeff, iteration, noise_levels):
        """
        One iteration of residual signal extraction using L0 hard thresholding.

        Decomposes the current residual (data - model) into wavelet coefficients,
        applies hard thresholding restricted to the previously detected support mask,
        and reconstructs the extracted residual signal to update the model.

        Parameters
        ----------
        model : np.ndarray
            Current denoised cube estimate.
        mask_coeff : np.ndarray, dtype bool
            Support mask from the initial hard threshold pass — True where a
            coefficient was detected as significant in the data.
        iteration : int
            Current iteration index (controls first-iteration verbose messages).
        noise_levels : list of float
            Per-sub-band noise standard deviation estimates.

        Returns
        -------
        model : np.ndarray
            Updated model after adding the extracted residual signal.
        delta : np.ndarray
            Residual signal component added in this iteration.
        """
        max_voxel_index = np.argmax(self._signal)
        iz, _, _ = np.unravel_index(max_voxel_index, self._signal.shape)

        residual = self._data - model
        thresh = self._threshold_level

        if self._plot:
            fig, axs = plt.subplots(2, 3, figsize=(16, 13), constrained_layout=True)

            im1 = axs[0, 0].imshow(model[iz], vmin=np.min(self._signal[iz]), vmax=np.max(self._signal[iz]), cmap='RdBu_r')
            axs[0, 0].set_title('Previously Denoised (Iteration #{})'.format(iteration))

            im2 = axs[0, 1].imshow((self._signal - model)[iz], vmin=np.min(self._signal[iz]), vmax=np.max(self._signal[iz]), cmap='RdBu_r')
            axs[0, 1].set_title('SIGNAL Residual')

            im3 = axs[0, 2].imshow(self._signal[iz], vmin=np.min(self._signal[iz]), vmax=np.max(self._signal[iz]), cmap='RdBu_r')
            axs[0, 2].set_title('SIGNAL')

            for ax in axs[0]:
                ax.axis('off')

            cbar1 = fig.colorbar(im1, ax=axs[0, 0], orientation='horizontal', fraction=0.05, pad=0.02)
            cbar2 = fig.colorbar(im2, ax=axs[0, 1], orientation='horizontal', fraction=0.05, pad=0.02)
            cbar3 = fig.colorbar(im3, ax=axs[0, 2], orientation='horizontal', fraction=0.05, pad=0.02)
            cbar1.set_label('Flux')
            cbar2.set_label('Flux')
            cbar3.set_label('Flux')

        if iteration == 0:
            if self._verbose: print('(*) Decomposing residual into wavelet scales')

        inds, shapes, w_residual = self.mr2d1d.decompose(
            residual, self._num_scales_2d, self._num_scales_1d
        )

        if iteration == 0:
            if self._verbose: print('(*) Hard thresholding residual within detected support')

        i = 0
        for scale2d in range(self._num_scales_2d):
            for scale1d in range(self._num_scales_1d):
                start, end = inds[scale2d][scale1d]

                if scale2d == self._num_scales_2d - 1 and scale1d == self._num_scales_1d - 1:
                    continue

                c_data = w_residual[start:end]
                noise_level = noise_levels[i]
                mask = mask_coeff[start:end].astype(bool)

                # Keep residual coefficients that are both in the support AND above
                # the hard threshold; zero everything else.
                above_thresh = np.abs(c_data) > thresh * noise_level
                keep = mask & above_thresh
                w_residual[start:end] = np.where(keep, c_data, 0.0)

                i += 1

        if iteration == 0:
            if self._verbose: print('(*) Reconstructing the new signal coefficients into the real space')

        delta = self.mr2d1d.reconstruct(w_residual)

        if iteration == 0:
            if self._verbose: print('(*) Updating the model with the newly detected signal')

        model = model + delta

        if self._plot:
            im4 = axs[1, 0].imshow(residual[iz], cmap='RdBu_r')
            axs[1, 0].set_title('Residual')

            im5 = axs[1, 1].imshow(delta[iz], cmap='RdBu_r', vmin=self._signal[iz].min(), vmax=self._signal[iz].max())
            axs[1, 1].set_title('Residual Information')

            im6 = axs[1, 2].imshow(model[iz], cmap='RdBu_r', vmin=self._signal[iz].min(), vmax=self._signal[iz].max())
            axs[1, 2].set_title('Updated Model (Iteration #{})'.format(iteration + 1))

            for ax in axs[1]:
                ax.axis('off')

            cbar4 = fig.colorbar(im4, ax=axs[1, 0], orientation='horizontal', fraction=0.05, pad=0.02)
            cbar5 = fig.colorbar(im5, ax=axs[1, 1], orientation='horizontal', fraction=0.05, pad=0.02)
            cbar6 = fig.colorbar(im6, ax=axs[1, 2], orientation='horizontal', fraction=0.05, pad=0.02)
            cbar4.set_label('Flux')
            cbar5.set_label('Flux')
            cbar6.set_label('Flux')

            plt.subplots_adjust(hspace=1)
            plt.show()

        return model, delta

    def _denoise_iterative_hard(self, num_iter=20, patience=3):
        """
        Perform iterative hard thresholding denoising with L0 regularization.

        Applies a single hard threshold pass on the data to establish the signal
        support and initial model, then iteratively extracts residual signal within
        that support using L0 hard thresholding.

        Parameters
        ----------
        num_iter : int, optional
            Maximum number of debias iterations (default 20).

        Returns
        -------
        best_model : np.ndarray
            Final denoised cube at the iteration with lowest residual STD.
        deltas : np.ndarray
            Accumulated residual signal extracted during debias iterations.
        residual_stds : list of float
            Residual STD at each debias iteration.
        best_iteration : int
            Iteration number where the best model was achieved.
        noise_levels : list of float
            Estimated noise levels for each wavelet sub-band.
        """
        if self._verbose:
            print('----[ Denoising with ITERATIVE HARD THRESHOLDING ]----')

        if self._verbose:
            print('(*) Decomposing data and estimating per-band noise levels')
        inds, _, w_data_fixed, noise_levels = self._decompose_and_estimate_noise()
        if self._verbose:
            src = 'noise cube (std)' if self._noise is not None else 'data (MAD)'
            print(f'    Noise estimated from {src} across {len(noise_levels)} sub-bands')

        max_voxel_index = np.argmax(self._signal)
        iz = np.unravel_index(max_voxel_index, self._signal.shape)[0]

        # ---- Phase 1: single hard threshold pass --------------------------------
        if self._verbose:
            print(f'(*) Applying hard thresholding (lambda = {self._threshold_level} sigma)')

        thresh = self._threshold_level
        w_filtered = w_data_fixed.copy()
        mask_coeff = np.zeros_like(w_data_fixed, dtype=bool)
        i = 0
        for scale2d in range(self._num_scales_2d):
            for scale1d in range(self._num_scales_1d):
                start, end = inds[scale2d][scale1d]
                if scale2d == self._num_scales_2d - 1 and scale1d == self._num_scales_1d - 1:
                    continue
                t = thresh + (self._thresh_increm if scale2d == 0 else 0)
                band_mask = np.abs(w_data_fixed[start:end]) > t * noise_levels[i]
                mask_coeff[start:end] = band_mask
                w_filtered[start:end] = np.where(band_mask, w_data_fixed[start:end], 0.0)
                i += 1

        # Restore coarse band — always pass it unchanged so pysparse reconstruct()
        # uses the correct large-scale component from the data, not its internal cache.
        start_c, end_c = inds[self._num_scales_2d - 1][self._num_scales_1d - 1]
        w_filtered[start_c:end_c] = w_data_fixed[start_c:end_c]

        if self._verbose:
            print('(*) Reconstructing initial model from thresholded data')

        model_1_step = self._apply_positivity(
            self.mr2d1d.reconstruct(np.ascontiguousarray(w_filtered, dtype=np.float32))
        )

        if self._verbose:
            residual_std_init = np.std(self._data - model_1_step)
            print(f'(*) Aperture Flux: {np.sum(model_1_step):.3e}, Residual STD: {residual_std_init:.3e}')

        # ---- Phase 2: iterative residual signal extraction ----------------------
        if self._verbose:
            print('\n----[ DEBIAS ]----\n')
            print('(*) Iteratively extracting remaining signal from residual')

        p_init = patience
        epsilon = 5e-4
        converged = False
        global_min_residual_std = np.inf
        global_best_model = None
        global_best_iteration = 0

        for p in range(p_init, -1, -1):
            if self._verbose:
                print(f'\n[*] Trying with plateau condition: {p} consecutive stable residuals needed for convergence')

            plateau_counter = 0
            previous_residual_std = 1e-33
            deltas = np.zeros_like(self._data)
            residual_stds = []
            model = model_1_step.copy()

            for iteration in range(num_iter):

                if self._verbose:
                    print(f'\n\n--- [ DEBIAS ITERATION #{iteration + 1} ] ---\n')

                model, delta = self._residual_signal_extraction_l0(
                    model, mask_coeff, iteration, noise_levels
                )
                model = self._apply_positivity(model)
                deltas += delta

                aperture_flux = np.sum(model)
                residual_std = np.std(self._data - model)
                residual_stds.append(residual_std)

                if self._verbose:
                    print(f'(*) Aperture Flux: {aperture_flux:.3e}, Residual STD: {residual_std:.3e}')

                if residual_std < global_min_residual_std:
                    global_min_residual_std = residual_std
                    global_best_model = model.copy()
                    global_best_iteration = iteration + 1

                if p > 0:
                    if abs(residual_std - previous_residual_std) / previous_residual_std <= epsilon:
                        plateau_counter += 1
                    else:
                        plateau_counter = 0

                    if plateau_counter >= p:
                        if self._verbose:
                            print(f'\nflux: {aperture_flux:.3e}')
                            print(f'noise: {residual_std:.3e}')
                            print(f'Convergence achieved at iteration #{iteration + 1} with p = {p}')
                        converged = True
                        break

                previous_residual_std = residual_std

                if iteration == 0 and self._verbose:
                    print('(*) Repeating these steps until convergence')

            if converged:
                break

        best_model = global_best_model if global_best_model is not None else model
        best_iteration = global_best_iteration

        if not converged:
            if self._verbose:
                print(f'[Warning] Convergence not achieved for any p value from {p_init} to 0')
                print(f'Using best model at iteration #{best_iteration} with residual std = {global_min_residual_std:.3e}')

        if self._verbose:
            print(f'\nBest iteration : {best_iteration}')
            print(f'Final resid STD: {global_min_residual_std:.4e}')

        if self._plot:
            plt.figure(figsize=(28, 11))
            plt.subplot(121)
            plt.imshow(self._data[iz], cmap='RdBu_r')
            plt.title('Noisy Data')
            plt.colorbar()

            plt.subplot(122)
            plt.imshow(self._signal[iz], cmap='RdBu_r', vmin=self._signal[iz].min(), vmax=self._signal[iz].max())
            plt.title('Clean Signal')
            plt.colorbar()
            plt.show()

            plt.figure(figsize=(28, 9))
            plt.subplot(131)
            plt.imshow(model_1_step[iz], cmap='RdBu_r', vmin=self._signal[iz].min(), vmax=self._signal[iz].max())
            plt.title('One-Step Denoising')
            plt.colorbar()
            plt.axis('off')

            plt.subplot(132)
            plt.imshow(np.maximum(deltas, 0)[iz], cmap='RdBu_r', vmin=self._signal[iz].min(), vmax=self._signal[iz].max())
            plt.title('Residual Signal')
            plt.colorbar()
            plt.axis('off')

            plt.subplot(133)
            plt.imshow(best_model[iz], cmap='RdBu_r', vmin=self._signal[iz].min(), vmax=self._signal[iz].max())
            plt.title('Final Denoised')
            plt.colorbar()
            plt.axis('off')
            plt.show()

        return best_model, deltas, residual_stds, best_iteration, noise_levels

    def _residual_signal_extraction_l1(self, model, mask_coeff, all_weights, iteration, noise_levels):
        """
        Perform residual signal extraction using L1 soft-thresholding with adaptive weights.

        This method implements one iteration of weighted soft-thresholding for the debiasing
        step in iterative soft-thresholding algorithms. Unlike hard thresholding that uses
        binary masks, this approach applies adaptive weights that account for the bias
        introduced by soft-thresholding in previous iterations.

        Parameters
        ----------
        model : np.ndarray
            Current estimate of the denoised cube from previous iteration.
        mask_coeff : np.ndarray
            Boolean mask indicating coefficients previously identified as significant.
        all_weights : np.ndarray
            Adaptive weights for soft-thresholding computed from previous iterations.
        iteration : int
            Current iteration number for progress reporting and conditional processing.
        noise_levels : list of float
            Pre-computed noise level estimates for each wavelet sub-band.

        Returns
        -------
        model : np.ndarray
            Updated model after incorporating newly extracted residual signal.
        delta : np.ndarray
            Extracted signal component added to model in this iteration.
        """
        max_voxel_index = np.argmax(self._signal)
        iz, _, _ = np.unravel_index(max_voxel_index, self._signal.shape)

        residual = self._data - model
        thresh = self._threshold_level
        if self._plot:
            fig, axs = plt.subplots(2, 3, figsize=(16, 13), constrained_layout=True)

            im1 = axs[0,0].imshow(model[iz], vmin = np.min(self._signal[iz]), vmax = np.max(self._signal[iz]), cmap = 'RdBu_r')
            axs[0,0].set_title('Previously Denoised (Iteration #{})'.format(iteration))

            im2 = axs[0,1].imshow((self._signal - model)[iz], vmin = np.min(self._signal[iz]), vmax = np.max(self._signal[iz]), cmap = 'RdBu_r')
            axs[0,1].set_title('SIGNAL Residual')

            im3 = axs[0,2].imshow(self._signal[iz], vmin = np.min(self._signal[iz]), vmax = np.max(self._signal[iz]), cmap = 'RdBu_r')
            axs[0,2].set_title('SIGNAL')

            axs[0,0].axis('off')
            axs[0,1].axis('off')
            axs[0,2].axis('off')

            cbar1 = fig.colorbar(im1, ax=axs[0, 0], orientation='horizontal', fraction=0.05, pad=0.02)
            cbar2 = fig.colorbar(im2, ax=axs[0, 1], orientation='horizontal', fraction=0.05, pad=0.02)
            cbar3 = fig.colorbar(im3, ax=axs[0, 2], orientation='horizontal', fraction=0.05, pad=0.02)

            cbar1.ax.tick_params()
            cbar2.ax.tick_params()
            cbar3.ax.tick_params()

            cbar1.set_label('Flux')
            cbar2.set_label('Flux')
            cbar3.set_label('Flux')

        if iteration == 0:
            if self._verbose: print('(*) Decomposing residual into wavelet scales')

        inds, shapes, w_residual = self.mr2d1d.decompose(residual,
                                                            self._num_scales_2d,
                                                            self._num_scales_1d)

        if iteration == 0:
            if self._verbose: print('(*) Performing weighted de-biasing with previously calculated weights')

        i = 0
        for scale2d in range(self._num_scales_2d):
            for scale1d in range(self._num_scales_1d):

                start, end = inds[scale2d][scale1d]

                # Skip the coarse approximation band — leave it unchanged so low-frequency
                # structure from the residual is not lost or distorted.
                if scale2d == self._num_scales_2d - 1 and scale1d == self._num_scales_1d - 1:
                    continue

                # Use residual coefficients directly — no median centering, which would
                # suppress genuine DC components (e.g. diffuse emission plateaus).
                c_data = w_residual[start:end]

                noise_level = noise_levels[i]

                mask = mask_coeff[start:end].astype(bool)
                weights = all_weights[start:end]

                # Weighted soft-threshold on previously detected coefficients (mask=True).
                shrink = weights[mask] * thresh * noise_level
                w_residual[start:end][mask] = np.sign(c_data[mask]) * np.maximum(
                    np.abs(c_data[mask]) - shrink, 0.0
                )

                w_residual[start:end][~mask] = 0.0

                i += 1

        if iteration == 0:
            if self._verbose: print('(*) Reconstructing the new signal coefficients into the real space')

        # Inverse wavelet transform — delta may be positive or negative; do NOT clip here.
        # Clipping to positive then subtracting the mean (to undo artificial flux) is contradictory.
        # Instead, leave delta free-signed: the reconstruction of thresholded residual coefficients
        # is naturally near-zero-mean. Positivity is enforced on the full model at the call site.
        delta = self.mr2d1d.reconstruct(w_residual)

        if iteration == 0:
            if self._verbose: print('(*) Updating the model with the newly detected signal')

        model = model + delta
        if self._plot:
            im4 = axs[1,0].imshow(residual[iz], cmap = 'RdBu_r')
            axs[1,0].set_title('Residual')

            im5 = axs[1,1].imshow(delta[iz], cmap = 'RdBu_r', vmin = self._signal[iz].min(), vmax = self._signal[iz].max())
            axs[1,1].set_title('Residual Information')

            im6 = axs[1,2].imshow(model[iz], cmap = 'RdBu_r', vmin = self._signal[iz].min(), vmax = self._signal[iz].max())
            axs[1,2].set_title('Updated Model (Iteration #{})'.format(iteration+1))

            cbar4 = fig.colorbar(im4, ax=axs[1, 0], orientation='horizontal', fraction=0.05, pad=0.02)
            cbar5 = fig.colorbar(im5, ax=axs[1, 1], orientation='horizontal', fraction=0.05, pad=0.02)
            cbar6 = fig.colorbar(im6, ax=axs[1, 2], orientation='horizontal', fraction=0.05, pad=0.02)

            cbar4.ax.tick_params()
            cbar5.ax.tick_params()
            cbar6.ax.tick_params()

            cbar4.set_label('Flux')
            cbar5.set_label('Flux')
            cbar6.set_label('Flux')

            axs[1,0].axis('off')
            axs[1,1].axis('off')
            axs[1,2].axis('off')

            plt.subplots_adjust(hspace=1)
            plt.show()

        return model, delta

    def _denoise_iterative_soft(self, num_iter_reweight=20, num_iter_debias=20, debias=True, patience=3):
        """
        Perform iterative soft-thresholding denoising on the 3D data cube.

        This method applies multiple re-weighting iterations followed by an optional
        debiasing step to recover the underlying signal from noisy data. It uses a
        2D-1D multiscale wavelet decomposition for denoising, adaptive thresholding,
        and plateau-based convergence criteria.

        Parameters
        ----------
        num_iter_reweight : int, optional
            Number of iterations for the re-weighting denoising step (default is 20).
        num_iter_debias : int, optional
            Number of iterations for the debiasing step to extract residual signal (default is 20).
        debias : bool, optional
            Whether to perform the debiasing step (default is True).

        Returns
        -------
        best_model : np.ndarray
            Final denoised model after all iterations.
        model_1_step : np.ndarray
            Model after the first re-weighted denoising iteration.
        model_no_reweight : np.ndarray
            Model obtained without re-weighting in the first iteration.
        deltas : np.ndarray
            Accumulated residual signals extracted during debiasing.
        residual_stds_reweight : list of float
            Standard deviations of residuals during the re-weighting step.
        residual_stds_debias : list of float
            Standard deviations of residuals during the debiasing step.
        best_iteration : int
            Iteration number where the best model (lowest residual std) was achieved.
        dists : list of np.ndarray
            Selected wavelet sub-band distributions for diagnostic purposes.
        noise_levels : list of float
            Estimated noise levels for each wavelet sub-band.
        """
        if self._verbose:
            print('----[ Denoising with ITERATIVE SOFT THRESHOLDING ]----')

        if self._verbose:
            print('(*) Decomposing data and estimating per-band noise levels')
        inds, _, w_data_fixed, noise_levels = self._decompose_and_estimate_noise()
        if self._verbose:
            src = 'noise cube (std)' if self._noise is not None else 'data (MAD)'
            print(f'    Noise estimated from {src} across {len(noise_levels)} sub-bands')

        self.mean, self.std = np.mean(self._data), np.std(self._data)
        model = self._data.copy()
        thresh = self._threshold_level

        model_no_reweight = self._data.copy()

        p_init = patience

        max_voxel_index = np.argmax(self._signal)
        iz = np.unravel_index(max_voxel_index, self._signal.shape)[0]

        dists = []
        for scale2d in range(self._num_scales_2d):
            for scale1d in range(self._num_scales_1d):
                start, end = inds[scale2d][scale1d]
                if (scale2d == 5) and (scale1d == 0):
                    dists.append(w_data_fixed[start:end])

        converged_reweight = False
        global_min_residual_std_rw = np.inf
        global_best_model_rw = None
        global_best_iteration_rw = 0
        for p in range(p_init, -1, -1):
            if self._verbose:
                print(f'\n[*] Trying with plateau condition: {p} consecutive stable residuals needed for convergence')

            plateau_counter = 0
            previous_residual_std = 1e-33
            epsilon = 1e-3

            flux_history = []
            residual_stds_reweight = []

            for iteration in range(num_iter_reweight):

                if self._verbose:
                    print(f'\n\n--- [ DE-NOISING ITERATION #{iteration + 1} ] ---\n')

                # Decompose current model to compute per-band weights for adaptive thresholding.
                # The ISTA proximal step (mu=0.5 on ||y-x||^2) always maps back to the noisy
                # data, so w_data is w_data_fixed. We copy it to avoid modifying the reference.
                if iteration == 0:
                    if self._verbose: print('(*) Decomposing model for weight computation')
                _, _, w_data_weights = self.mr2d1d.decompose(
                    model, self._num_scales_2d, self._num_scales_1d
                )

                w_data = w_data_fixed.copy()

                if self._verbose and iteration == 1:
                    print('(*) Thresholding noisy data with updated adaptive weights (ISTA proximal step)')
                    print('(*) Calculating weights to account for soft-thresholding bias')

                w_data_copy = w_data.copy()
                mask_coeff = np.zeros_like(w_data, dtype=bool)
                all_weights = np.ones_like(w_data)

                i = 0
                for scale2d in range(self._num_scales_2d):
                    for scale1d in range(self._num_scales_1d):

                        start, end = inds[scale2d][scale1d]

                        if scale2d == self._num_scales_2d - 1 and scale1d == self._num_scales_1d - 1:
                            continue

                        c_data = w_data[start:end]

                        noise_level = noise_levels[i]

                        noise_level_weight = self._estimate_noise(w_data_weights[start:end])

                        # Apply the same high-frequency increment as IHT: finest 2D scale
                        # (scale2d==0) is the noisiest; a higher threshold here prevents
                        # spurious noise peaks from being reweighted to full amplitude.
                        t = thresh + (self._thresh_increm if scale2d == 0 else 0)
                        mask = np.abs(c_data) > t * noise_level

                        if iteration == 0:
                            weights = np.ones_like(c_data)
                        else:
                            # Epsilon scaled to noise prevents weight blow-up for near-zero coefficients
                            weights = thresh * noise_level_weight / (
                                np.abs(w_data_weights[start:end]) + noise_level_weight * 1e-6
                            )

                        all_weights[start:end] = weights
                        mask_coeff[start:end] = mask
                        shrink = weights[mask] * t * noise_level
                        w_data[start:end][mask] = np.sign(c_data[mask]) * np.maximum(
                            np.abs(c_data[mask]) - shrink, 0.0
                        )

                        # Coefficients below the detection threshold are always zero after soft-thresholding
                        # (max(|c| - thresh*noise, 0) == 0 by definition of the mask), so zero explicitly.
                        w_data[start:end][~mask] = 0.0

                        i += 1

                        if (scale2d == 5) and (scale1d == 0):
                            dists.append(w_data[start:end])

                        if self._plot:
                            if (scale2d==5) and (scale1d==0):

                                denosied_dist = w_data[start:end]
                                threshold_noise = thresh * noise_level

                                bins = np.linspace(w_data_copy[start:end].min(), w_data_copy[start:end].max(),100)

                                plt.figure(figsize = (11,7))
                                plt.hist(w_data[start:end], bins = bins, color = 'xkcd:blue', alpha = 0.5, label = 'Denoised')
                                plt.hist(w_data_copy[start:end], bins = bins, histtype='step', color = 'black', alpha = 1, label = 'Original')
                                plt.axvline(thresh * noise_level, color = 'black', linestyle = 'dashed', label = '{:.1f}'.format(self._threshold_level)+r'$\sigma$' + ' Threshold')
                                plt.axvline(-thresh * noise_level, color = 'black', linestyle = 'dashed')
                                plt.yscale('log')
                                plt.ylim(0,5e5)
                                plt.ylabel('$N_{C_{ij}}$')
                                plt.xlabel('$C_{ij}$')
                                plt.legend()
                                plt.grid(True)
                                plt.show()

                if iteration==0:
                        if self._verbose: print('(*) Reconstructing the new signal coefficients into the real space')

                model_denoised = self.mr2d1d.reconstruct(np.ascontiguousarray(w_data, dtype=np.float32))
                model_print = self.mr2d1d.reconstruct(np.ascontiguousarray(w_data_copy, dtype=np.float32))

                if iteration==0:
                        if self._verbose: print('(*) Applying the positivity constraint')
                model_denoised = self._apply_positivity(model_denoised)

                if self._plot:
                    plt.figure(figsize = (15,12))
                    plt.subplot(221)
                    plt.imshow(model_print[iz], cmap = 'RdBu_r',  vmin=self._data[iz].min(), vmax=self._data[iz].max())
                    plt.colorbar()
                    plt.axis('off')
                    plt.title('Input')

                    plt.subplot(222)
                    plt.imshow(model_denoised[iz], cmap='RdBu_r', vmin=self._signal[iz].min(), vmax=self._signal[iz].max())
                    plt.colorbar()
                    plt.axis('off')
                    plt.title('Denoised Iteration #{}'.format(iteration+1))

                    plt.figure(figsize = (15,12))
                    plt.subplot(223)
                    plt.imshow(self._signal[iz], cmap = 'RdBu_r',  vmin=self._signal[iz].min(), vmax=self._signal[iz].max())
                    plt.colorbar()
                    plt.axis('off')
                    plt.title('Signal')

                    plt.subplot(224)
                    plt.imshow((self._signal - model_denoised)[iz], cmap='RdBu_r', vmin=self._signal[iz].min(), vmax=self._signal[iz].max())
                    plt.colorbar()
                    plt.axis('off')
                    plt.title('SIGNAL Residual')

                    plt.show()

                if iteration==0:
                    if self._verbose: print('(*) Repeating these steps for subsequent iterations')

                model = model_denoised

                aperture_flux = np.sum(model)
                residual_std = np.std(self._data - model)

                residual_stds_reweight.append(residual_std)

                if self._verbose:
                    print(f"(*) Aperture Flux: {aperture_flux:.3e}, Clean Flux: {np.sum(self._signal):.3e}, Residual STD: {residual_std:.3e}")

                if residual_std < global_min_residual_std_rw:
                    global_min_residual_std_rw = residual_std
                    global_best_model_rw = model.copy()
                    global_best_iteration_rw = iteration + 1

                if p > 0:
                    if abs(residual_std - previous_residual_std) / previous_residual_std <= epsilon:
                        plateau_counter += 1
                    else:
                        plateau_counter = 0

                    if plateau_counter >= p:
                        if self._verbose:
                            print(f'\nflux: {aperture_flux}')
                            print(f'noise: {residual_std}')
                            print(f'Re-weight Convergence achieved at iteration #{iteration + 1} with p = {p}')
                        converged_reweight = True
                        break

                previous_residual_std = residual_std

                if iteration == 1 and self._verbose:
                    print(f'(*) Repeating these steps for subsequent {num_iter_reweight - 2} iterations')

                if iteration == 0:
                    model_no_reweight = model_denoised

            if converged_reweight:
                break

        best_model = global_best_model_rw if global_best_model_rw is not None else model
        best_iteration = global_best_iteration_rw

        if not converged_reweight:
            if self._verbose:
                print(f'[Warning] Re-weight convergence not achieved for any p value from {p_init} to 0')
                print(f'Using best model at iteration #{best_iteration} with residual std = {global_min_residual_std_rw:.3e}')

        if self._verbose: print('\n----[ DE-BIASING ]----\n')
        if self._verbose:   print('(*) Iteratively extracting remaining signal from residual')

        p_init_debias = patience
        epsilon_debias = 5e-4
        model_1_step = best_model.copy()

        converged_debias = False
        global_min_residual_std_db = np.inf
        global_best_model_db = None
        global_best_iteration_db = 0
        for p in range(p_init_debias, -1, -1):
            if self._verbose:
                print(f'\n[*] Trying with plateau condition: {p} consecutive stable residuals needed for convergence')

            plateau_counter = 0
            previous_residual_std = 1e-33
            deltas = np.zeros_like(model)
            residual_stds_debias = []

            model = model_1_step.copy()

            for iteration in range(num_iter_debias):

                model, delta = self._residual_signal_extraction_l1(
                    model, mask_coeff, all_weights, iteration, noise_levels
                )
                model = self._apply_positivity(model)
                deltas += delta

                aperture_flux = np.sum(model)
                residual_std = np.std(self._data - model)

                if self._verbose:
                    print(f"(*) Aperture Flux: {aperture_flux:.3e}, Residual STD: {residual_std:.3e}")

                flux_history.append(aperture_flux)
                residual_stds_debias.append(residual_std)

                if residual_std < global_min_residual_std_db:
                    global_min_residual_std_db = residual_std
                    global_best_model_db = model.copy()
                    global_best_iteration_db = iteration + 1

                if p > 0:
                    if abs(residual_std - previous_residual_std) / previous_residual_std <= epsilon_debias:
                        plateau_counter += 1
                    else:
                        plateau_counter = 0

                    if plateau_counter >= p:
                        if self._verbose:
                            print(f'\nflux: {aperture_flux:.3e}')
                            print(f'noise: {residual_std:.3e}')
                            print(f'Debias convergence achieved at iteration #{iteration + 1} with p = {p}')
                        converged_debias = True
                        break

                previous_residual_std = residual_std

                if iteration == 0 and self._verbose:
                    print('(*) Repeating these steps until convergence')

            if converged_debias:
                break

        best_model = global_best_model_db if global_best_model_db is not None else model
        best_iteration = global_best_iteration_db

        if not converged_debias:
            if self._verbose:
                print(f'[Warning] Debias convergence not achieved for any p value from {p_init_debias} to 0')
                print(f'Using best model at iteration #{best_iteration} with residual std = {global_min_residual_std_db:.3e}')

        if self._plot:
            plt.figure(figsize=(28, 11))
            plt.subplot(121)
            plt.imshow(self._data[iz], cmap='RdBu_r')
            plt.title('Noisy Data')
            plt.colorbar()

            plt.subplot(122)
            plt.imshow(self._signal[iz], cmap='RdBu_r', vmin=self._signal[iz].min(), vmax=self._signal[iz].max())
            plt.title('Clean Signal')
            plt.colorbar()
            plt.show()

            plt.figure(figsize=(28, 9))
            plt.subplot(131)
            plt.imshow(model_1_step[iz], cmap='RdBu_r', vmin=self._signal[iz].min(), vmax=self._signal[iz].max())
            plt.title('One-Step Denoising')
            plt.colorbar()
            plt.axis('off')

            plt.subplot(132)
            plt.imshow(np.maximum(deltas, 0)[iz], cmap='RdBu_r', vmin=self._signal[iz].min(), vmax=self._signal[iz].max())
            plt.title('Residual Signal')
            plt.colorbar()
            plt.axis('off')

            plt.subplot(133)
            plt.imshow(best_model[iz], cmap='RdBu_r', vmin=self._signal[iz].min(), vmax=self._signal[iz].max())
            plt.title('Final Denoised')
            plt.colorbar()
            plt.axis('off')
            plt.show()

        return best_model, model_1_step, model_no_reweight, deltas, residual_stds_reweight, residual_stds_debias, best_iteration, dists, noise_levels

    def _compute_emission_rmse(self, model):
        """
        Compute Root Mean Square Error in emission regions only.

        Parameters
        ----------
        model : np.ndarray
            Reconstructed/denoised data cube with same shape as original signal.

        Returns
        -------
        float
            Root mean square error between masked signal and model.
        """
        masked_diff = self._mask * self._signal - self._mask * model
        return np.sqrt(np.mean(masked_diff ** 2))

    def _estimate_noise(self, array):
        """
        Estimate noise standard deviation using robust Median Absolute Deviation (MAD).

        Parameters
        ----------
        array : np.ndarray
            Array of values (typically wavelet coefficients) for noise estimation.

        Returns
        -------
        float
            Estimated noise standard deviation (1.48 * MAD).
        """
        median_val = np.median(array)
        mad = np.median(np.abs(array - median_val))
        # Factor 1.48 ≈ 1/Φ^(-1)(3/4) where Φ^(-1) is inverse normal CDF
        return 1.48 * mad


def mock_noise_value(mock_cube, peak_snr):
    """
    Calculate noise level for synthetic data cube based on desired peak SNR.

    Parameters
    ----------
    mock_cube : np.ndarray
        Clean synthetic data cube without noise.
    peak_snr : float
        Desired signal-to-noise ratio at the peak of the cube.

    Returns
    -------
    float
        Noise standard deviation to achieve the desired peak SNR.
    """
    peak_signal = np.max(mock_cube)
    mock_cube_noise = peak_signal / peak_snr
    print(f'Max SNR: {peak_snr}')
    print(f'Mock noise level: {mock_cube_noise:.6e}')
    return mock_cube_noise
