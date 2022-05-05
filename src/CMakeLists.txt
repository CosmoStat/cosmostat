##################
# COSMOSTAT 2020 #
##################

cmake_minimum_required(VERSION 3.12.0)
project(cosmostat)

# Include modules in CMake module path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/modules/")

# Load CMake tools
include(ExternalProject)
include(FindPkgConfig)

# Load custom CMake function
# - find_pkg
# - build_package
# - build_pybind_list
include(functions)

# Set build type
if(CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE "Release")
endif()
message(STATUS "CMake Build Type: ${CMAKE_BUILD_TYPE}")

# Set RPATH options
set(CMAKE_INSTALL_RPATH "${CMAKE_BINARY_DIR}/lib")

# ----------------- #
# Find Dependencies #
# ----------------- #

# Locate OpenMP
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
  find_package(BigMac 0.0.1 REQUIRED)
else()
  find_package(OpenMP REQUIRED)
endif()

# Set compilation flags
set(CMAKE_CXX_FLAGS
  "${CMAKE_CXX_FLAGS} -w -DNO_DISP_IO -fPIC -Wno-write-strings -std=c++17\
  ${OpenMP_CXX_FLAGS} ${BigMac_NOWARN}"
)
set(CMAKE_CXX_FLAGS_RELEASE
  "${CMAKE_CXX_FLAGS_RELEASE} ${CMAKE_CXX_FLAGS} -g0 -fomit-frame-pointer"
)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS} -O0")

# Build CFITSIO
include(BuildCFITSIO)

# Build Sparse2D
include(BuildSparse2D)

# Set basic package dependency list including Sparse2D, CFITSIO and FFTW
set(basic_libs ${sparse2d_libs} ${CFITSIO_LIBRARIES})

# -------------- #
# Build packages #
# -------------- #

# ---- #
# MISC #
# ---- #

# List package libs
set(misc_libs misc1d misc2d misc3d)
# List package mains
set(misc_mains "")
# List package dependencies
# - Sparse2D
# - CFITSIO
# - FFTW
set(misc_deps ${basic_libs})
# Build package
build_package(
  "${misc_mains}"       #List of mains
  "${misc_libs}"        #List of libs
  "${misc_deps}"        #List of dependencies
  misc                  #Path to package
  cc                    #Source extension
  h                     #Header extension
)

# -- #
# MR #
# -- #

# List package libs
set(mr_libs mr1d mr2d)
# List package mains
set(mr_mains "")
# List package dependencies
# - MISC
# - Sparse2D
# - CFITSIO
# - FFTW
set(mr_deps ${misc_libs} ${basic_libs})
# Build package
build_package(
  "${mr_mains}"         #List of mains
  "${mr_libs}"          #List of libs
  "${mr_deps}"          #List of dependencies
  mr                    #Path to package
  cc                    #Source extension
  h                     #Header extension
)

# -- #
# MC #
# -- #

# List package libs
set(mc_libs mc1d mc2d)
# List package mains
set(mc_mains "")
# List package dependencies
# - MR
# - MISC
# - Sparse2D
# - CFITSIO
# - FFTW
set(mc_deps ${mr_libs} ${misc_libs} ${basic_libs})
# Build package
build_package(
  "${mc_mains}"         #List of mains
  "${mc_libs}"          #List of libs
  "${mc_deps}"          #List of dependencies
  mc                    #Path to package
  cc                    #Source extension
  h                     #Header extension
)

# --- #
# MGA #
# --- #

# List package libs
set(mga_libs mga1d)
# List package mains
set(mga_mains mgamain1d)
# List package dependencies
# - MC
# - MR
# - MISC
# - Sparse2D
# - CFITSIO
# - FFTW
set(mga_deps ${mc_libs} ${mr_libs} ${misc_libs} ${basic_libs})
# Build package
build_package(
  "${mga_mains}"        #List of mains
  "${mga_libs}"         #List of libs
  "${mga_deps}"         #List of dependencies
  mga                   #Path to package
  cc                    #Source extension
  h                     #Header extension
)
