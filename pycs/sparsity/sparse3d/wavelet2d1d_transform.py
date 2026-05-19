import numpy as np
import pysparse


class Wavelet2D1DTransform(object):
    """
    2D-1D Multi-Scale Wavelet Transform for 3D Spectral Cubes

    This class implements a hybrid wavelet decomposition that combines 2D spatial
    starlet transforms with 1D spectral wavelet analysis. The approach is specifically
    designed for IFU spectral cubes where spatial and spectral information should be
    treated differently to preserve both morphological and kinematic features.

    The decomposition treats the 3D cube as a stack of 2D images (spatial planes)
    along the spectral axis, applying:
    1. 2D starlet transform to capture spatial structure at multiple scales
    2. 1D wavelet transform along the spectral direction for each spatial pixel

    This preserves the multi-scale nature of astronomical sources while maintaining
    spectral line profiles and velocity structure.

    Attributes
    ----------
    NOISE_TAB : np.ndarray
        Pre-computed noise scaling factors for different wavelet sub-bands.
        Shape (5, 5) for up to 5 spatial and 5 spectral scales.
        Currently not used but available for future noise modeling.

    Methods
    -------
    decompose(cube, num_scales_2d, num_scales_1d)
        Perform forward 2D-1D wavelet transform
    reconstruct(coeffs)
        Perform inverse transform to recover spatial cube
    energy_per_scale(num_scales_2d, num_scales_1d)
        Get noise scaling for requested number of scales

    Notes
    -----
    - The transform preserves flux conservation
    - Coefficients are organized by scale with metadata for reconstruction
    - Transform type and normalization controlled by pysparse.MR2D1D parameters
    """

    # Pre-computed noise scaling factors for different wavelet sub-bands
    # TODO: Compute values more accurately by excluding borders of each sub-band cube
    NOISE_TAB = np.array([[0.9610, 0.9261, 0.9611, 0.9810, 0.9933],
                          [0.2368, 0.2282, 0.2369, 0.2417, 0.2449],
                          [0.1049, 0.1011, 0.1050, 0.1071, 0.1082],
                          [0.0527, 0.0507, 0.0528, 0.0539, 0.0543],
                          [0.0283, 0.0272, 0.0283, 0.0288, 0.0295]])

    def __init__(self, transform_type=2):
        """
        Initialize the 2D-1D wavelet transform with specified parameters.

        This is a wrapper for pysparse's MR2D1D transform, providing a simplified
        interface for 3D spectral cube analysis. The transform combines 2D spatial
        decomposition with 1D spectral analysis to preserve both morphological
        and kinematic information in astronomical data cubes.

        Parameters
        ----------
        transform_type : int, optional
            Type of wavelet transform to perform with `decompose`, by default 2.
            Available options from pysparse documentation:
            - 2: 2D starlet (undecimated) + 1D 7/9 filter bank (decimated) [default]

        Notes
        -----
        The transform object (self._mr2d1d) is created lazily when decompose()
        is first called. This allows automatic parameter detection based on
        input cube dimensions and avoids unnecessary memory allocation.

        The default choice (transform_type=2) provides:
        - Isotropic 2D starlets: Excellent for point sources and extended emission
        - 1D 7/9 filter bank: Good spectral localization for line profiles
        - Undecimated spatial: Preserves spatial resolution and avoids artifacts
        - Decimated spectral: Reduces computational cost while preserving features

        References
        ----------
        See pysparse documentation for detailed descriptions of available transforms
        and their mathematical properties.
        """
        self.transform_type = transform_type
        self._mr2d1d = None       # Initialized lazily in decompose method
        self._cached_scales = None  # (num_scales_2d, num_scales_1d) of cached object

    def decompose(self, cube, num_scales_2d, num_scales_1d):
        """
        Perform forward 2D-1D multi-scale wavelet decomposition of spectral cube.

        Decomposes the input 3D cube into wavelet coefficients using a hybrid approach:
        - 2D spatial decomposition preserves morphological structure at multiple scales
        - 1D spectral decomposition captures velocity/frequency features
        - Combined analysis exploits correlations between spatial and spectral dimensions

        The decomposition creates (num_scales_2d × num_scales_1d) sub-bands, each
        representing specific combinations of spatial and spectral scales. This allows
        for scale-dependent processing during denoising operations.

        Parameters
        ----------
        cube : np.ndarray, shape (nz, ny, nx)
            Input 3D spectral data cube. The spectral axis should be first (axis=0),
            followed by spatial dimensions. Typical astronomical convention where
            each cube[i,:,:] slice is a 2D image at a specific wavelength/velocity.
        num_scales_2d : int
            Number of scales for 2D spatial starlet decomposition.
            Must be >= 2 and <= int(log2(min(ny, nx))).
            Controls spatial resolution analysis - more scales capture finer details.
        num_scales_1d : int
            Number of scales for 1D spectral wavelet decomposition.
            Must be >= 2 and <= int(log2(nz)).
            Controls spectral feature resolution - more scales capture narrower lines.

        Returns
        -------
        inds : list of list of tuples
            Nested index structure for accessing coefficients by scale.
            inds[i][j] = (start, end) gives indices for sub-band (i,j) in coeffs array.
            Use as: coeffs[start:end] to extract coefficients for scale pair (i,j).
        shapes : list of list of tuples
            Nested shape information for each wavelet sub-band.
            shapes[i][j] = (nx, ny, nz) gives dimensions of sub-band (i,j).
        coeffs : np.ndarray, 1D
            Flattened array containing all wavelet coefficients plus metadata.
            Structure: [n_scales_2d, n_scales_1d, shape_info, coefficients...]
            Required for reconstruction via reconstruct() method.

        Notes
        -----
        Transform object is initialized on first call. Subsequent calls reuse
        the same transform for efficiency. Coefficients are organized in
        increasing order of scale (fine to coarse).

        See Also
        --------
        reconstruct : Inverse transform to recover original cube
        _extract_metadata : Internal method for parsing coefficient structure
        """
        # Reuse the cached MR2D1D object when scale parameters are unchanged.
        # Recreating it on every call causes repeated large buffer allocations
        # inside the C extension, which can exhaust memory or segfault in iterative loops.
        scales = (num_scales_2d, num_scales_1d)
        if self._mr2d1d is None or self._cached_scales != scales:
            self._mr2d1d = pysparse.MR2D1D(type_of_transform=self.transform_type,
                                           normalize=False,
                                           verbose=False,
                                           NbrScale2d=num_scales_2d,
                                           Nbr_Plan=num_scales_1d)
            self._cached_scales = scales

        # Ensure input is C-contiguous native float32 as required by the C extension
        cube = np.ascontiguousarray(cube, dtype=np.float32)

        # Perform the forward wavelet transform
        coeffs = self._mr2d1d.transform(cube)

        # Extract index mapping and shape information for coefficient organization
        # This metadata is essential for accessing specific sub-bands during processing
        inds, shapes = self._extract_metadata(coeffs)

        return inds, shapes, coeffs

    def reconstruct(self, coeffs):
        """
        Perform inverse 2D-1D wavelet transform to recover original cube.

        Reconstructs the 3D spectral cube from its multi-scale wavelet representation.
        This is the inverse operation of decompose(), combining information from all
        spatial and spectral scales to recover the original data structure.

        The reconstruction process:
        1. Parses the coefficient structure from decompose()
        2. Applies inverse 1D spectral wavelets
        3. Applies inverse 2D spatial starlets
        4. Combines all scales to form the final cube

        Parameters
        ----------
        coeffs : np.ndarray, 1D
            Wavelet coefficients as returned by decompose().
            Must include metadata (scales, shapes) followed by coefficient values.
            Can be modified coefficients (e.g., after thresholding for denoising).

        Returns
        -------
        cube : np.ndarray, shape (nz, ny, nx)
            Reconstructed 3D spectral data cube.
            Shape matches the original input to decompose().
            Preserves flux conservation (sum of pixels) unless coefficients modified.

        Notes
        -----
        The transform object must be initialized (via decompose()) before calling
        this method. The reconstruction uses the same transform parameters as
        the forward decomposition.

        Modified coefficients (e.g., thresholded for denoising) will produce
        a reconstructed cube with altered characteristics while preserving
        the overall structure encoded in the retained coefficients.
        """
        # assert hasattr(self, '_mr2d1d'), "Need to call decompose first."
        reconstructed = self._mr2d1d.reconstruct(coeffs)
        return reconstructed

    def energy_per_scale(self, num_scales_2d, num_scales_1d):
        """
        Get pre-computed noise energy scaling factors for wavelet sub-bands.

        Returns the noise scaling factors that account for how noise propagates
        through the wavelet transform at different scales. These values can be
        used to normalize thresholds for denoising operations.

        Parameters
        ----------
        num_scales_2d : int
            Number of 2D spatial scales to retrieve factors for.
        num_scales_1d : int
            Number of 1D spectral scales to retrieve factors for.

        Returns
        -------
        np.ndarray, shape (num_scales_2d, num_scales_1d)
            Noise scaling factors for each (2D scale, 1D scale) combination.
            Values represent the multiplicative factor for noise standard deviation
            in each wavelet sub-band relative to the input noise level.

        Notes
        -----
        Currently uses pre-computed values in NOISE_TAB. Future versions should
        compute these more accurately by excluding border effects in each sub-band.
        """
        return self.NOISE_TAB[:num_scales_2d, :num_scales_1d]

    @property
    def num_precomputed(self):
        """
        Get the maximum number of scales for which noise factors are pre-computed.

        Returns
        -------
        tuple of int
            (max_2d_scales, max_1d_scales) available in NOISE_TAB.
        """
        return self.NOISE_TAB.shape

    def _extract_metadata(self, coeffs):
        """
        Parse coefficient array structure to extract indexing and shape metadata.

        The coefficients array from pysparse contains metadata followed by the actual
        wavelet coefficients. This method parses the structure to create convenient
        indexing arrays for accessing specific sub-bands during processing.

        Coefficient Array Structure:
        - coeffs[0]: Number of 2D scales used
        - coeffs[1]: Number of 1D scales used
        - coeffs[2:5]: Shape (nx,ny,nz) for sub-band (0,0)
        - coeffs[5:5+nx*ny*nz]: Coefficients for sub-band (0,0)
        - coeffs[...]: Shape and coefficients for remaining sub-bands
        - Sub-bands ordered by: for scale2d in range(n_scales_2d):
                                   for scale1d in range(n_scales_1d): ...

        Parameters
        ----------
        coeffs : np.ndarray, 1D
            Flattened coefficient array from decompose() containing metadata
            and coefficients for all wavelet sub-bands.

        Returns
        -------
        inds : list of list of tuples
            Index structure where inds[i][j] = (start, end) gives the slice
            indices for accessing coefficients of sub-band (2d_scale=i, 1d_scale=j).
            Use as: coeffs[start:end] to get coefficients for that sub-band.
        shapes : list of list of tuples
            Shape structure where shapes[i][j] = (nx, ny, nz) gives the 3D
            dimensions of sub-band (2d_scale=i, 1d_scale=j). Essential for
            reshaping flattened coefficients back to 3D arrays.

        Notes
        -----
        The returned structures are essential for:
        - Accessing specific wavelet sub-bands for thresholding
        - Reshaping coefficients for 3D processing
        - Reconstruction operations
        """
        # Extract number of scales from the first two elements
        n_scales_2d = int(coeffs[0])
        n_scales_1d = int(coeffs[1])

        # Initialize nested lists for index ranges and shapes
        # Structure: inds[2d_scale][1d_scale] = (start_idx, end_idx)
        inds = [[() for _ in range(n_scales_1d)] for _ in range(n_scales_2d)]
        shapes = [[() for _ in range(n_scales_1d)] for _ in range(n_scales_2d)]

        # Parse through the coefficient array to extract metadata for each sub-band
        start = end = 2  # Skip the first two scale count elements

        # Traverse all scale combinations in order used by pysparse
        for ii in range(n_scales_2d):
            for jj in range(n_scales_1d):
                # Each sub-band starts with 3 shape values (nx, ny, nz)
                start = end + 3

                # Extract sub-band dimensions
                nx, ny, nz = map(int, coeffs[start-3 : start])
                shapes[ii][jj] = (nx, ny, nz)

                # Calculate coefficient count and ending index for this sub-band
                ncoeff = nx * ny * nz
                end = start + ncoeff

                # Store the index range for this sub-band
                inds[ii][jj] = (start, end)

        return inds, shapes
