# CosmoStat
Software package for cosmostatistics

This branch is stripped down to provide the binaries `cb_mca1d` and `cb_mmca1d` exclusively.

## Requirements

- CFITSIO (installed automatically)
- CMake
- OpenMP

### macOS Set Up

To install the package requirements it is recommended to use [Homebrew](https://brew.sh/).

```bash
brew tap sfarrens/sf
brew install bigmac libomp
```

Note that this package installs `CFITSIO` from source. Alternative
installations of this package may cause issues.

## Install

To install the package:

```bash
git clone -b cb_mca1d --depth 1 git@github.com:CosmoStat/shapepipe.git
cd cosmostat
mkdir build
cd build
cmake ../src/
make
```

The binaries will be available in the `build` directory by default. These binaries can be installed to the path of your choosing by passing the `-DCMAKE_INSTALL_PREFIX:PATH` option to CMake. e.g.

```bash
cmake ../src/ -DCMAKE_INSTALL_PREFIX:PATH=<PATH>
make
make install
```

where `<PATH>` is the path of your choice.
