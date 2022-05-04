#========================================================#
# Build the CfitsIO dependencies for the project         #
#========================================================#

# Set CFITSIO Version
# set(CFITSIOVersion 4.0.0)
# set(cfitsioSHA256 b2a8efba0b9f86d3e1bd619f662a476ec18112b4f27cc441cc680a4e3777425e)

set(CFITSIOVersion 3.49)
set(cfitsioSHA256 5b65a20d5c53494ec8f638267fca4a629836b7ac8dd0ef0266834eab270ed4b3)

# Set C/C++ compiler and flags
set(CFITSIO_COMPILE
  CC=${CMAKE_C_COMPILER}
  CXX=${CMAKE_CXX_COMPILER}
)

# Set CFITSIO configuration flags
set(CFITSIO_CONFIG_FLAGS
  --prefix=${CMAKE_BINARY_DIR}
  --enable-reentrant
)

# Download and build CFITSIO
ExternalProject_Add(cfitsio
    PREFIX cfitsio
    URL  http://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio-${CFITSIOVersion}.tar.gz
    URL_HASH  SHA256=${cfitsioSHA256}
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND ${CFITSIO_COMPILE} ./configure ${CFITSIO_CONFIG_FLAGS}
    BUILD_COMMAND make shared
    INSTALL_COMMAND make install
)

set(CFITSIO_LIBRARY_DIRS ${CMAKE_BINARY_DIR}/lib/ )
set(CFITSIO_INCLUDE_DIRS ${CMAKE_BINARY_DIR}/include/ )
set(CFITSIO_LIBRARIES -lcfitsio)

include_directories(${CFITSIO_INCLUDE_DIRS})
link_directories(${CFITSIO_LIBRARY_DIRS})
