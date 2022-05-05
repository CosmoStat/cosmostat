#=================#
# Build Sparse2D  #
#=================#

# Set Sparse2D Version
set(Sparse2DVersion v2.1.7)

# Set Sparse2D source dir
set(SPARSE2D_SOURCE ${CMAKE_CURRENT_BINARY_DIR}/sparse2d)

# Check for CMAKE_PREFIX_PATH
if(NOT "${CMAKE_PREFIX_PATH}" STREQUAL "")
  # Create a list with an alternate separator e.g. pipe symbol
  string(REPLACE ";" "|" SPARSE2D_PREFIX_PATH "${CMAKE_PREFIX_PATH}")
endif()

if(NOT "${SPARSE2D_PREFIX_PATH}" STREQUAL "")
  message(WARNING "Using this prefix Path for Sparse2D build: ${SPARSE2D_PREFIX_PATH}")
endif()

# Download and build Sparse2D
ExternalProject_Add(sparse2d-git
  GIT_REPOSITORY https://github.com/CosmoStat/Sparse2D.git
  GIT_TAG ${Sparse2DVersion}
  PREFIX sparse2d
  LIST_SEPARATOR | # Use the alternate list separator
  CMAKE_ARGS ${CMAKE_ARGS}
  -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
  -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
  -DBUILD_MSVST=ON
  -DUSE_FFTW=OFF
  -DBUILD_FFTW=OFF
  -DBUILD_NFFT=OFF
  -DCMAKE_PREFIX_PATH:PATH=${CMAKE_BINARY_DIR}
  -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}
  -DCMAKE_INSTALL_RPATH=${CMAKE_BINARY_DIR}/lib
  DEPENDS cfitsio
)

# List Sparse2D libs
set(sparse2d_libs mga2d sparse3d sparse2d sparse1d tools)


# Include external headers
include_directories(${CMAKE_INSTALL_PREFIX}/include/sparse2d)

set(external_libs ${sparse2d_libs})

# Add external libraries
foreach(library ${external_libs})
  add_library(${library} STATIC IMPORTED)
  set_property(TARGET ${library} PROPERTY IMPORTED_LOCATION ${CMAKE_INSTALL_PREFIX}/lib/lib${library}.a)
  add_dependencies(${library} sparse2d-git)
endforeach()
