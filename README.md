# cosmostat
Software package for cosmostatistics

## Requirements

- Armadillo
- CFITSIO
- CMake
- FFTW
- GSL
- HEALPix
- OpenMP
- Python

### macOS Set Up

To install the package requirements it is recommended to use [Homebrew](https://brew.sh/).

```bash
brew tap sfarrens/sf
brew install armadillo bigmac cfitsio fftw gsl healpix libomp
```

## Install

To install without thinking:
```
$ pip install git+https://github.com/CosmoStat/cosmostat.git
```

To install by cloning:
```
$ git clone git@github.com:CosmoStat/cosmostat.git
$ cd cosmostat
$ pip install -e .
```

## Example

```python
import numpy as np
import matplotlib.pyplot as plt
import pycs

# Generate a test image
img = np.random.randn(256, 256)

# Take the starlet transform with 5 wavelet scales
st = pycs.astro.wl.image.transforms.starlet2d(img, nscales=5)

# Compute the aperture mass map at scale 4 using the starlet filter
apm = pycs.astro.wl.image.filters.aperture_mass(img, theta=2**4, filter='starlet')

# Plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6, 5))
ax1.imshow(st[3], cmap='magma')
ax2.imshow(apm, cmap='magma')
for ax in (ax1, ax2, ax3):
    ax.set_axis_off()
```
