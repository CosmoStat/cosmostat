#================#
# Build HEALPix  #
#================#

# Set HEALPix Version
set(HEALPixVersion 3.80)
set(HEALPixFileDate 2021Jun22)
set(HEALPixHash 923d31845716014e38f34c4de59264e1)

# set(HEALPixVersion 3.60)
# set(HEALPixFileDate 2019Dec18)
# set(HEALPixHash 9b51b2fc919f4e70076d296826eebee0)

# Set C/C++ compiler and flags
set(HEALPIX_COMPILE
  CC=${CMAKE_C_COMPILER}
  CXX=${CMAKE_CXX_COMPILER}
  FITSDIR=${CFITSIO_LIBRARY_DIRS}
  FITSINC=${CFITSIO_INCLUDE_DIRS}
)

# Set FFTW configuration flags
set(HEALPIX_CONFIG_FLAGS
  -L
  --auto=profile,sharp,c,cxx
)

# Download and build FFTW
ExternalProject_Add(healpix
    URL               https://sourceforge.net/projects/healpix/files/Healpix_${HEALPixVersion}/Healpix_${HEALPixVersion}_${HEALPixFileDate}.tar.gz
    URL_HASH          MD5=${HEALPixHash}
    SOURCE_DIR        ${CMAKE_BINARY_DIR}/healpix
    BUILD_IN_SOURCE   1
    CONFIGURE_COMMAND ${HEALPIX_COMPILE} ./configure ${HEALPIX_CONFIG_FLAGS}
    BUILD_COMMAND     make -j 4
    COMMAND           make clean
    DEPENDS           cfitsio
)

set(HEALPIX_LIBRARY_DIRS ${CMAKE_BINARY_DIR}/healpix/lib/ )
set(HEALPIX_INCLUDE_DIRS ${CMAKE_BINARY_DIR}/healpix/include/healpix_cxx/ )
set(HEALPIX_LIBRARIES -lhealpix_cxx)

include_directories(${HEALPIX_INCLUDE_DIRS})
link_directories(${HEALPIX_LIBRARY_DIRS})
