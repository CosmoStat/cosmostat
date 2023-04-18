<a href="http://www.cosmostat.org/" target_="blank"><img src="http://www.cosmostat.org/wp-content/uploads/2017/07/CosmoStat-Logo_WhiteBK.jpg" width="400"></a>

# CosmoStat

Software developed at the [CosmoStat](https://www.cosmostat.org/) lab at [CEA](https://www.cea.fr/english) Paris-Saclay.

## Basic installation

The package can be installed with `pip` as follows.

```bash
python -m pip install .
```

### Requirements

This package installs [Sparse2D](https://github.com/CosmoStat/Sparse2D) as a backend, which requires the follow dependencies.

- [Armadillo](https://arma.sourceforge.net/)
- [BigMac](https://github.com/sfarrens/bigmac) (macOS only)
- [CFITSIO](https://heasarc.gsfc.nasa.gov/fitsio/)
- [CMake](https://cmake.org/)
- [FFTW](https://www.fftw.org/)
- [GSL](https://www.gnu.org/software/gsl/)
- [HEALPix](https://healpix.jpl.nasa.gov/)
- [OpenMP](https://www.openmp.org/)
- [Python](https://www.python.org/)

### macOS Set Up

To install the package requirements on macOS, it is recommended to use [Homebrew](https://brew.sh/).

```bash
brew tap sfarrens/sf
brew install armadillo bigmac cfitsio fftw gsl healpix libomp
```

## Docker installation

### Pull the Docker image

If you have [Docker](https://www.docker.com/) installed, you can pull the latest build of the [CosmoStat image](https://github.com/cosmostat/cosmostat/pkgs/container/cosmostat) as follows:

```bash
docker pull ghcr.io/cosmostat/cosmostat:master
```

No further installation is required.

### Run a Docker container

To run a container on data in your current working directory, simply run:

```bash
docker run -v ${PWD}:/workdir --rm ghcr.io/cosmostat/cosmostat:master
```

The reference to `${PWD}` can be replaced by the path to any directory on your system.

Additionally, you can run a Sparse2D executables dirextly from the image. For example, to run a bspline wavelet transform on a FITS image called `myfile.fits` you would run:

```bash
docker run -v ${PWD}:/workdir --rm ghcr.io/cosmostat/cosmostat:master mr_transform -t 2 myfile.fits myoutput.mr
```

> Tip: If you don't want to constantly write the full Docker run command you can create an alias *e.g.*:
> ```bash
> alias cosmostat="docker run -v ${PWD}:/workdir --rm ghcr.io/cosmostat/cosmostat:master"
> ```
> then you can simply run *e.g.*:
> ```bash
> cosmostat mr_transform -h
> ```

### Launch a Jupyter notebook

You can also run a Jupyter notebook with a Docker container as the backend in order to use the cosmostat package.

```bash
docker run -p 8888:8888 -v ${PWD}:/workdir --rm ghcr.io/cosmostat/cosmostat:master notebook
```
