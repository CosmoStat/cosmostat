##################
# COSMOSTAT 2020 #
##################

cmake_minimum_required(VERSION 3.0.0)
project(cosmostat)

# Include modules in CMake module path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/modules/")

# Load CMake tools
include(ExternalProject)
include(FindPkgConfig)

# Load custom CMake function
include(functions)

# Set build type
if(CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE "Release")
endif()
message(STATUS "CMake Build Type: ${CMAKE_BUILD_TYPE}")

# ----------------- #
# Find Dependencies #
# ----------------- #

# Locate OpenMP
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
  find_package(BigMac 0.0.1 REQUIRED)
else()
  find_package(OpenMP REQUIRED)
endif()

# Locate CFITSIO
find_cfitsio()


# Set Sparse2D source dir
set(SPARSE2D_SOURCE ${CMAKE_CURRENT_BINARY_DIR}/sparse2d)

# Build Sparse2D
include(BuildSparse2D)

# set(FFTW_CXX_FLAGS "-DUSE_FFTW")

# Set compilation flags
set(CMAKE_CXX_FLAGS
  "${CMAKE_CXX_FLAGS} -DNO_DISP_IO -fPIC -Wno-write-strings\
  ${OpenMP_CXX_FLAGS} ${BigMac_NOWARN} ${FFTW_CXX_FLAGS}"
)
set(CMAKE_CXX_FLAGS_RELEASE
  "${CMAKE_CXX_FLAGS_RELEASE} -g0 -fomit-frame-pointer"
)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")

# -------- #
# Astro WL #
# -------- #

set(astro_wl_path astro/weaklensing)

# Build astro wl and wls libs
set(astro_wl_libs divcurlwt wl)
foreach(library ${astro_wl_libs})
  build_lib(${library} ${astro_wl_path} cc h)
endforeach()

# Find all astro wl targets
find_targets(astro_wl_targets ${astro_wl_path}/mainwl cc)

# # Build binaries
foreach(program ${astro_wl_targets})
  build_bin(${program} "${astro_wl_libs}" ${astro_wl_path}/mainwl cc)
endforeach(program)

# Install libraries
install(TARGETS ${astro_wl_libs} DESTINATION lib)

# Install binaraies
install(TARGETS ${astro_wl_targets} DESTINATION bin)
