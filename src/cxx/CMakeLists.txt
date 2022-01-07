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

# ---- #
# MISC #
# ---- #

# Build libs
set(misc_libs misc1d misc2d misc3d)
build_libs_cch("${misc_libs}" misc)

# Find all targets
find_targets(misc_targets misc/main2d cc)
find_targets(misc_3d_targets misc/main3d cc)

# Build binaries
build_bins_cc("${misc_targets}" "${misc_libs}" misc/main2d)
build_bins_cc("${misc_3d_targets}" "${misc_libs}" misc/main3d)

# Install libraries and binaries
install(TARGETS ${misc_libs} DESTINATION lib)
install(TARGETS ${misc_targets} DESTINATION bin)

# -- #
# MR #
# -- #

# Build libs
set(mr_libs mr1d mr2d)
build_libs_cch("${mr_libs}" mr)

# Set list of all required libs
set(mr_all_libs ${mr_libs} ${misc_libs})

# Find all targets
find_targets(mr1d_targets mr/mrmain1d cc)
find_targets(mr2d_targets mr/mrmain2d cc)
find_targets(mr3d_targets mr/mrmain3d cc)

# Build binaries
build_bins_cc("${mr1d_targets}" "${mr_all_libs}" mr/mrmain1d)
build_bins_cc("${mr2d_targets}" "${mr_all_libs}" mr/mrmain2d)
build_bins_cc("${mr3d_targets}" "${mr_all_libs}" mr/mrmain3d)

# Install libraries and binaries
install(TARGETS ${mr_libs} DESTINATION lib)
install(TARGETS ${mr1d_targets} ${mr2d_targets} ${mr3d_targets} DESTINATION bin)

# -- #
# MC #
# -- #

# Build libs
set(mc_libs mc1d mc2d)
build_libs_cch("${mc_libs}" mc)

# Set list of all required libs
set(mc_all_libs ${mc_libs} ${mr_libs} ${misc_libs})

# Find all targets
find_targets(mc1d_targets mc/mcmain1d cc)
find_targets(mc2d_targets mc/mcmain2d cc)

# Build binaries
build_bins_cc("${mc1d_targets}" "${mc_all_libs}" mc/mcmain1d)
build_bins_cc("${mc2d_targets}" "${mc_all_libs}" mc/mcmain2d)

# Install libraries and binaries
install(TARGETS ${mc_libs} DESTINATION lib)
install(TARGETS ${mc1d_targets} ${mc2d_targets} DESTINATION bin)

# --- #
# MGA #
# --- #

# Build libs
set(mga_libs mga1d mga3d)
build_libs_cch("${mga_libs}" mga)

# Set list of all required libs
set(mga_all_libs ${mga_libs} ${mc_libs} ${mr_libs} ${misc_libs})

# Find all targets
find_targets(mga1d_targets mga/mgamain1d cc)
find_targets(mga2d_targets mga/mgamain2d cc)
find_targets(mga3d_targets mga/mgamain3d cc)

# Build binaries
build_bins_cc("${mga1d_targets}" "${mga_all_libs}" mga/mgamain1d)
build_bins_cc("${mga2d_targets}" "${mga_all_libs}" mga/mgamain2d)
build_bins_cc("${mga3d_targets}" "${mga_all_libs}" mga/mgamain3d)

# Install libraries and binaries
install(TARGETS ${mga_libs} DESTINATION lib)
install(TARGETS ${mga1d_targets} ${mga2d_targets} ${mga3d_targets} DESTINATION bin)

# ---- #
# MWIR #
# ---- #

option(BUILD_MWIR "Build MWIR library" ON)
message(STATUS "MWIR Build: ${BUILD_MWIR}")

if(BUILD_MWIR)

  # Build libs
  set(mwir_libs mwir)
  build_libs_cch("${mwir_libs}" mwir)

  # Set list of all required libs
  set(mwir_all_libs ${mwir_libs} ${misc_libs})

  # Find all targets
  find_targets(mwir_targets mwir/mwirmain cc)

  # Build binaries
  build_bins_cc("${mwir_targets}" "${mwir_all_libs}" mwir/mwirmain)

  # Install libraries and binaries
  install(TARGETS ${mwir_libs} DESTINATION lib)
  install(TARGETS ${mwir_targets} DESTINATION bin)

endif()

# ------------------- #
# DICTIONARY LEARNING #
# ------------------- #

# ISSUE: Need to to include GSL and Armadillo libraries

option(BUILD_DICLEARN "Build Dictionary Learning library" ON)
message(STATUS "BUILD_DICLEARN Build: ${BUILD_DICLEARN}")

if(BUILD_DICLEARN)

  # Locate GSL
  find_pkg(GSL gsl)

  # Locate Armadillo
  find_pkg(Armadillo Armadillo)

  # Build libs
  set(learn_libs diclearn)
  build_libs_cch("${diclearn_libs}" diclearn)

  # Set list of all required libs
  set(diclearn_all_libs ${diclearn_libs} ${misc_libs} ${GSL_LIBRARIES} ${ARMADILLO_LIBRARIES})

  # Find all targets
  find_targets(diclearn_targets diclearn/diclearnmain cc)

  # Build binaries
  build_bins_cc("${diclearn_targets}" "${diclearn_all_libs}" diclearn/diclearnmain)

  # Install libraries and binaries
  install(TARGETS ${diclearn_libs} DESTINATION lib)
  install(TARGETS ${diclearn_targets} DESTINATION bin)

endif()

# --- #
# MRS #
# --- #

# Build libs
set(mrs_libs mrs mrsp)
build_libs_cch("${mrs_libs}" mrs)

# Find all targets
find_targets(mrs_targets mrs/mainmrs cc)
find_targets(mrsp_targets mrs/mainmrsp cc)

# Build binaries
build_bins_cc("${mrs_targets}" "${mrs_libs}" mrs/mainmrs)
build_bins_cc("${mrsp_targets}" "${mrs_libs}" mrs/mainmrsp)

# Install libraries and binaries
install(TARGETS ${mrs_libs} DESTINATION lib)
install(TARGETS ${mrs_targets} ${mrsp_targets} DESTINATION bin)

# -------- #
# Astro WL #
# -------- #

option(BUILD_ASTRO_WL "Build ASTRO WL library" ON)
message(STATUS "ASTRO WL Build: ${BUILD_ASTRO_WL}")

if(BUILD_ASTRO_WL)

  # Set path to astro weak lensing
  set(astro_wl_path astro/weaklensing)

  # Build libs
  set(astro_wl_libs divcurlwt wl wls)
  build_libs_cch("${astro_wl_libs}" ${astro_wl_path})

  # Set list of all required libs
  set(astro_wl_all_libs ${astro_wl_libs} ${mrs_libs})

  # Find all targets
  find_targets(astro_wl_targets ${astro_wl_path}/mainwl cc)
  find_targets(astro_wls_targets ${astro_wl_path}/mainwls cc)

  # Build binaries
  build_bins_cc("${astro_wl_targets}" "${astro_wl_all_libs}" ${astro_wl_path}/mainwl)
  build_bins_cc("${astro_wls_targets}" "${astro_wl_all_libs}" ${astro_wl_path}/mainwls)

  # Install libraries and binaries
  install(TARGETS ${astro_wl_libs} DESTINATION lib)
  install(TARGETS ${astro_wl_targets} ${astro_wls_targets} DESTINATION bin)

endif()

# --------- #
# Astro Gal #
# --------- #

option(BUILD_ASTRO_GAL "Build ASTRO GAL library" ON)
message(STATUS "ASTRO GAL Build: ${BUILD_ASTRO_GAL}")

if(BUILD_ASTRO_GAL)

  # Set path to astro galaxies
  set(astro_gal_path astro/galaxies)

  # Build libs
  set(astro_gal_libs gal)
  build_libs_cch("${astro_gal_libs}" ${astro_gal_path})

  # Set list of all required libs
  set(astro_gal_all_libs ${astro_gal_libs} ${misc_libs})

  # Find all targets
  find_targets(astro_gal_targets ${astro_gal_path}/galmain cc)

  # Build binaries
  build_bins_cc("${astro_gal_targets}" "${astro_gal_all_libs}" ${astro_gal_path}/galmain)

  # Install libraries and binaries
  install(TARGETS ${astro_gal_libs} DESTINATION lib)
  install(TARGETS ${astro_gal_targets} DESTINATION bin)

endif()

# --------------- #
# Python Bindings #
# --------------- #
option(BUILD_BINDING "Build BINDING library" ON)
message(STATUS "BINDING Build: ${BUILD_BINDING}")

if(BUILD_BINDING)

  # Include Pybind11 headers
  if(NOT DEFINED PYBIND11_INCLUDE_DIR)
    message(FATAL_ERROR "Path to Pybind11 headers must be specified.")
  else()
    include_directories(${PYBIND11_INCLUDE_DIR})
  endif()

  # Define binding targets
  set(binding_targets pymrs)

  # Build bindings
  foreach(program ${binding_targets})
    build_pybind(${program} python_binding cpp)
  endforeach(program)

endif()
