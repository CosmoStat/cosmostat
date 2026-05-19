"""3D sparsity routines for IFU spectral cube denoising.

Usage
-----
    from pycs.sparsity.sparse3d.wavelet2d1dtransform import Wavelet2D1DTransform
    from pycs.sparsity.sparse3d.denoiser2d1d import Denoiser2D1D

    denoiser = Denoiser2D1D(threshold_type='soft')
    result = denoiser.denoise(noisy_cube, signal_cube, method='iterative')
"""

from pycs.sparsity.sparse3d.wavelet2d1d_transform import Wavelet2D1DTransform
from pycs.sparsity.sparse3d.denoiser2d1d import Denoiser2D1D, mock_noise_value

__all__ = ["Wavelet2D1DTransform", "Denoiser2D1D", "mock_noise_value"]


