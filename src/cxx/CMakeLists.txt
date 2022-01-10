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

# ----------------- #
# Find Dependencies #
# ----------------- #

# Locate OpenMP
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
  find_package(BigMac 0.0.1 REQUIRED)
else()
  find_package(OpenMP REQUIRED)
endif()

# Locate Python
find_package(Python COMPONENTS Interpreter Development)
include_directories(${Python_INCLUDE_DIRS})
link_directories(${Python_LIBRARY_DIRS})

# Locate CFITSIO
find_pkg(CFITSIO cfitsio)

# Locate HEALPix
find_pkg(HEALPIX healpix_cxx>=3.80)

# Locate FFTW
find_package(FFTW3 QUIET)

# Build Sparse2D
include(BuildSparse2D)

# Set compilation flags
set(CMAKE_CXX_FLAGS
  "${CMAKE_CXX_FLAGS} -DNO_DISP_IO -fPIC -Wno-write-strings -std=c++11\
  ${OpenMP_CXX_FLAGS} ${BigMac_NOWARN}"
)
set(CMAKE_CXX_FLAGS_RELEASE
  "${CMAKE_CXX_FLAGS_RELEASE} -g0 -fomit-frame-pointer"
)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")

# Set basic package dependency list including Sparse2D, CFITSIO and FFTW
set(basic_libs ${sparse2d_libs} ${CFITSIO_LIBRARIES} ${fftw_libs})

# -------------- #
# Build packages #
# -------------- #

# ---- #
# MISC #
# ---- #

# List package libs
set(misc_libs misc1d misc2d misc3d)
# List package mains
set(misc_mains main2d main3d)
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
set(mr_mains mrmain1d mrmain2d mrmain3d)
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
set(mc_mains mcmain1d mcmain2d)
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
set(mga_libs mga1d mga3d)
# List package mains
set(mga_mains mgamain1d mgamain2d mgamain3d)
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

# ---- #
# MWIR #
# ---- #

option(BUILD_MWIR "Build MWIR library" ON)
message(STATUS "MWIR Build: ${BUILD_MWIR}")

if(BUILD_MWIR)

  # List package libs
  set(mwir_libs mwir)
  # List package mains
  set(mwir_mains mwirmain)
  # List package dependencies
  # - MISC
  # - Sparse2D
  # - CFITSIO
  # - FFTW
  set(mwir_deps ${misc_libs} ${basic_libs})
  # Build package
  build_package(
    "${mwir_mains}"       #List of mains
    "${mwir_libs}"        #List of libs
    "${mwir_deps}"        #List of dependencies
    mwir                  #Path to package
    cc                    #Source extension
    h                     #Header extension
  )

endif()

# ------------------- #
# DICTIONARY LEARNING #
# ------------------- #

# NOTE: Need to to include GSL and Armadillo libraries

option(BUILD_DICLEARN "Build Dictionary Learning library" ON)
message(STATUS "BUILD_DICLEARN Build: ${BUILD_DICLEARN}")

if(BUILD_DICLEARN)

  # Locate GSL using pkg-config
  pkg_check_modules(GSL gsl)
  find_pkg(GSL gsl)

  # Locate the Armadillo library
  find_package(Armadillo)
  include_directories(${ARMADILLO_INCLUDE_DIRS})

  # List package libs
  set(diclearn_libs diclearn)
  # List package mains
  set(diclearn_mains diclearnmain)
  # List package dependencies
  # - MR
  # - MISC
  # - GSL
  # - Armadillo
  # - Sparse2D
  # - CFITSIO
  # - FFTW
  set(diclearn_deps
    ${mr_libs}
    ${misc_libs}
    ${GSL_LIBRARIES}
    ${ARMADILLO_LIBRARIES}
    ${basic_libs}
  )
  # Build package
  build_package(
    "${diclearn_mains}"   #List of mains
    "${diclearn_libs}"    #List of libs
    "${diclearn_deps}"    #List of dependencies
    diclearn              #Path to package
    cc                    #Source extension
    h                     #Header extension
  )

endif()

# --- #
# MRS #
# --- #

# List package libs
set(mrs_libs mrs mrsp)
# List package mains
set(mrs_mains mainmrs mainmrsp)
# List package dependencies
  # - HEALPix
  # - Sparse2D
  # - CFITSIO
  # - FFTW
set(mrs_deps ${HEALPIX_LIBRARIES} ${basic_libs})
# Build package
build_package(
  "${mrs_mains}"        #List of mains
  "${mrs_libs}"         #List of libs
  "${mrs_deps}"         #List of dependencies
  mrs                   #Path to package
  cc                    #Source extension
  h                     #Header extension
)

# -------- #
# Astro WL #
# -------- #

option(BUILD_ASTRO_WL "Build ASTRO WL library" ON)
message(STATUS "ASTRO WL Build: ${BUILD_ASTRO_WL}")

if(BUILD_ASTRO_WL)

  # Set path to astro weak lensing
  set(astro_wl_path astro/weaklensing)
  # List package libs
  set(astro_wl_libs divcurlwt wl wls)
  # List package mains
  set(astro_wl_mains mainwl mainwls)
  # List package dependencies
  # - MRS
  # - HEALPix
  # - Sparse2D
  # - CFITSIO
  # - FFTW
  set(astro_wl_deps ${mrs_libs} ${HEALPIX_LIBRARIES} ${basic_libs})
  # Build package
  build_package(
    "${astro_wl_mains}"   #List of mains
    "${astro_wl_libs}"    #List of libs
    "${astro_wl_deps}"    #List of dependencies
    ${astro_wl_path}      #Path to package
    cc                    #Source extension
    h                     #Header extension
  )

endif()

# --------- #
# Astro Gal #
# --------- #

option(BUILD_ASTRO_GAL "Build ASTRO GAL library" ON)
message(STATUS "ASTRO GAL Build: ${BUILD_ASTRO_GAL}")

if(BUILD_ASTRO_GAL)

  # Set path to astro weak lensing
  set(astro_gal_path astro/galaxies)
  # List package libs
  set(astro_gal_libs gal)
  # List package mains
  set(astro_gal_mains galmain)
  # List package dependencies
  # - MISC
  # - Sparse2D
  # - CFITSIO
  # - FFTW
  set(astro_gal_deps ${misc_libs} ${basic_libs})
  # Build package
  build_package(
    "${astro_gal_mains}"  #List of mains
    "${astro_gal_libs}"   #List of libs
    "${astro_gal_deps}"   #List of dependencies
    ${astro_gal_path}     #Path to package
    cc                    #Source extension
    h                     #Header extension
  )

endif()

# --------------- #
# Python Bindings #
# --------------- #

option(BUILD_BINDING "Build Python Bindings" ON)
message(STATUS "Bindings Build: ${BUILD_BINDING}")

if(BUILD_BINDING)

  # Include Pybind11 headers
  if(NOT DEFINED PYBIND11_INCLUDE_DIR)
    message(FATAL_ERROR "Path to Pybind11 headers must be specified.")
  else()
    include_directories(${PYBIND11_INCLUDE_DIR})
  endif()

  # Define binding targets
  set(binding_targets pymrs)
  # List binding dependencies
  # - MRS
  # - HEALPix
  # - OpenMP
  # - Sparse2D
  # - CFITSIO
  # - FFTW
  set(pybind_deps
    ${mrs_libs}
    ${HEALPIX_LIBRARIES}
    ${OpenMP_CXX_LIBRARIES}
    ${basic_libs}
  )
  # Build bindings
  build_pybind_list(
    "${binding_targets}"  #List of bindings
    "${pybind_deps}"      #List of binding depdendencies
    python_binding        #Path to bindings
    cpp                   #Source extension
  )

endif()
