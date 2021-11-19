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
#find_pkg(HEALPIX healpix_cxx)
pkg_check_modules(HEALPIX healpix_cxx>=3.80)

include_directories(${HEALPIX_INCLUDE_DIRS})
if (HEALPIX_FOUND)
  include_directories(${HEALPIX_INCLUDE_DIRS})
  include_directories(${HEALPIX_INCLUDE_DIRS})
  include_directories(/usr/local/Cellar/healpix/3.80/include/healpix_cxx)
  link_directories(${HEALPIX_LIBRARY_DIRS})

  message("HEALPIX package: == YES == ")
  message("DIRLIB = ${HEALPIX_LINK_LIBRARIES}")
  message("INCLIB = ${HEALPIX_INCLUDE_DIRS}")
  message("${HEALPIX_LIBRARIES}")
else(HEALPIX_FOUND)
  if(DO_NOT_COMPILE_HEALPIX_BIN MATCHES "0")
    message("======================================================")
    message("NO HEALPIX package. MRS binaries will not be compiled.")
    message("======================================================")
  endif(DO_NOT_COMPILE_HEALPIX_BIN MATCHES "0")
#  set(DO_NOT_COMPILE_HEALPIX_BIN 1)
endif(HEALPIX_FOUND)


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

# --- #
# MR #
# --- #

# Build misc libs
set(misc_libs misc1d misc2d misc3d)
foreach(library ${misc_libs})
  build_lib(${library} misc cc h)
endforeach()
set(mr_libs mr1d mr2d)
foreach(library ${mr_libs})
  build_lib(${library} mr cc h)
endforeach()
set(mr_alllibs ${misc_libs} ${mr_libs})

# Find all mrs targets
find_targets(misc_targets misc/main2d cc)
find_targets(misc3d_targets misc/main3d cc)
find_targets(mr1d_targets mr/mrmain1d cc)
find_targets(mr2d_targets mr/mrmain2d cc)
find_targets(mr3d_targets mr/mrmain3d cc)

set(misc_all_targets ${misc_targets} ${misc3d_targets} ${mr1d_targets} ${mr2d_targets} ${mr3d_targets})

# Build  binaries
foreach(program ${misc_targets})
  build_bin(${program} "${misc_libs}" misc/main2d cc)
endforeach(program)
foreach(program ${misc3d_targets})
  build_bin(${program} "${misc_libs}" misc/main3d cc)
endforeach(program)
foreach(program ${mr1d_targets})
  build_bin(${program} "${mr_alllibs}" mr/mrmain1d cc)
endforeach(program)
foreach(program ${mr2d_targets})
  build_bin(${program} "${mr_alllibs}" mr/mrmain2d cc)
endforeach(program)
foreach(program ${mr3d_targets})
  build_bin(${program} "${mr_alllibs}" mr/mrmain3d cc)
endforeach(program)

# Install libraries
install(TARGETS ${misc_libs} DESTINATION lib)
install(TARGETS ${mr_libs} DESTINATION lib)

# Install binaraies
install(TARGETS ${misc_all_targets} DESTINATION bin)

# --- #
# MC #
# --- #

# Build misc libs
set(mc_libs mc1d mc2d)
foreach(library ${mc_libs})
  build_lib(${library} mr cc h)
endforeach()
set(mc_alllibs ${misc_libs} ${mr_libs} ${mc_libs})

# Find all mrs targets
find_targets(mc1d_targets mr/mcmain1d cc)
find_targets(mc2d_targets mr/mcmain2d cc)

set(mc_all_targets  ${mc1d_targets} ${mc2d_targets} )

# Build  binaries
foreach(program ${mc1d_targets})
  build_bin(${program} "${mc_alllibs}" mr/mcmain1d cc)
endforeach(program)
foreach(program ${mc2d_targets})
  build_bin(${program} "${mc_alllibs}" mr/mcmain2d cc)
endforeach(program)
 
# Install libraries
install(TARGETS ${mc_libs} DESTINATION lib)

# Install binaraies
install(TARGETS ${mc_all_targets} DESTINATION bin)

# --- #
# MGA #
# --- #

# Build misc libs
set(mga_libs mga1d mga3d)
foreach(library ${mga_libs})
  build_lib(${library} mga cc h)
endforeach()

set(mga_alllibs ${misc_libs} ${mr_libs} ${mc_alllibs} ${mga_libs})

# Find all mrs targets
find_targets(mga1d_targets mga/mgamain1d cc)
find_targets(mga2d_targets mga/mgamain2d cc)
find_targets(mga3d_targets mga/mgamain3d cc)

set(gma_all_targets  ${mga1d_targets} ${mga2d_targets} ${mga3d_targets})

# Build  binaries
foreach(program ${mga1d_targets})
  build_bin(${program} "${mga_alllibs}" mga/mgamain1d cc)
endforeach(program)
foreach(program ${mga2d_targets})
  build_bin(${program} "${mga_alllibs}" mga/mgamain2d cc)
endforeach(program)
foreach(program ${mga3d_targets})
  build_bin(${program} "${mga_alllibs}" mga/mgamain3d cc)
endforeach(program)

# Install libraries
install(TARGETS ${mga_alllibs} DESTINATION lib)

# Install binaraies
install(TARGETS ${gma_all_targets} DESTINATION bin)


# ------------- #
# MWIR #
# ------------- #

option(BUILD_MWIR "Build MWIR library" ON)
message(STATUS "MWIR Build: ${BUILD_MWIR}")

if(BUILD_MWIR)

# Build misc libs
set(mwir_libs mwir)
foreach(library ${mwir_libs})
  build_lib(${library} mr cc h)
endforeach()
set(mwir_alllibs ${misc_libs} ${mwir_libs})

# Find all mrs targets
find_targets(mwir_targets mr/mwirmain cc)

set(mwir_all_targets  ${mwir_targets})

# Build  binaries
foreach(program ${mwir_all_targets})
  build_bin(${program} "${mwir_alllibs}" mr/mwirmain cc)
endforeach(program)
 
# Install libraries
install(TARGETS ${mwir_libs} DESTINATION lib)

# Install binaraies
install(TARGETS ${mwir_all_targets} DESTINATION bin)

endif()

# ------------------- #
# DICTIONARY LEARNING #
# ------------------- #

# ISSUE NEED TO INCLUDE GSL LIBRARY

option(BUILD_DICLEARN "Build Dictionary Learning library" ON)
message(STATUS "BUILD_DICLEARN Build: ${BUILD_DICLEARN}")

if(BUILD_DICLEARN)

  # Locate GSL using pkg-config
  pkg_check_modules(GSL gsl)
  if (GSL_FOUND)
      message(STATUS "GSL_FOUND: ${BUILD_DICLEARN}")
      include_directories(${GSL_INCLUDE_DIRS})
      link_directories(${GSL_LIBRARY_DIRS})
  endif()
  # Locate the Armadillo library
  find_package(Armadillo)
  if (ARMADILLO_FOUND)
     message(STATUS "ARMADILLO_FOUND: ${ARMADILLO_FOUND} ${GSL_LIBRARIES} ${ARMADILLO_LIBRARIES} ")
     include_directories(${ARMADILLO_INCLUDE_DIRS})
  endif (ARMADILLO_FOUND)

 if (GSL_FOUND)
  if (ARMADILLO_FOUND)

    # Build misc libs
    set(learn_libs diclearn)
    foreach(library ${learn_libs})
        build_lib(${library} mr cc h)
    endforeach()
    set(learn_alllibs ${misc_libs} ${learn_libs} ${GSL_LIBRARIES} ${ARMADILLO_LIBRARIES} )

    # Install libraries
    install(TARGETS ${learn_libs} DESTINATION lib)
    
    # Find all mrs targets
    find_targets(learn_targets mr/diclearnmain cc)

    set(learn_all_targets  ${learn_targets})

    # Build  binaries
    foreach(program ${learn_all_targets})
        build_bin(${program} "${learn_alllibs}" mr/diclearnmain cc)
    endforeach(program)

    # Install binaraies
    install(TARGETS ${learn_all_targets} DESTINATION bin)

    endif()
  endif(GSL_FOUND)
endif()

# --- #
# MRS #
# --- #
option(BUILD_MRS "Build MRS library" ON)
message(STATUS "MRS Build: ${BUILD_MRS}")

if(BUILD_MRS)

# Build mrs libs
set(mrs_libs mrs mrsp)
foreach(library ${mrs_libs})
  build_lib(${library} mrs cc h)
endforeach()

# Find all mrs targets
find_targets(mrs_targets mrs/mainmrs cc)
find_targets(mrsp_targets mrs/mainmrsp cc)
set(mrs_all_targets ${mrs_targets} ${mrsp_targets})

# Build MRS binaries
foreach(program ${mrs_targets})
  build_bin(${program} "${mrs_libs}" mrs/mainmrs cc)
endforeach(program)

# Build MRSP binariess
foreach(program ${mrsp_targets})
  build_bin(${program} "${mrs_libs}" mrs/mainmrsp cc)
endforeach(program)

# Install libraries
install(TARGETS ${mrs_libs} DESTINATION lib)

# Install binaraies
install(TARGETS ${mrs_all_targets} DESTINATION bin)

endif()

# -------- #
# Astro WL #
# -------- #
option(BUILD_ASTRO_WL "Build ASTRO WL library" ON)
message(STATUS "ASTRO WL Build: ${BUILD_ASTRO_WL}")

if(BUILD_ASTRO_WL)
if(BUILD_MRS)

# Set path to astro WEAKLENSING
set(astro_wl_path astro/weaklensing)

# Build astro weaklensing libs
set(astro_wl_libs divcurlwt wl wls)
foreach(library ${astro_wl_libs})
  build_lib(${library} ${astro_wl_path} cc h)
  add_dependencies(${library} ${mrs_libs})
endforeach()
set(astro_wl_mrs_libs ${astro_wl_libs} ${mrs_libs})

# Find all astro wl targets
find_targets(astro_wl_targets ${astro_wl_path}/mainwl cc)
find_targets(astro_wls_targets ${astro_wl_path}/mainwls cc)
set(astro_wl_all_targets ${astro_wl_targets} ${astro_wls_targets})

# Build astro wl binaries
foreach(program ${astro_wl_targets})
  build_bin(${program} "${astro_wl_libs}" ${astro_wl_path}/mainwl cc)
endforeach(program)

# Build astro wls binaries
foreach(program ${astro_wls_targets})
  build_bin(${program} "${astro_wl_mrs_libs}" ${astro_wl_path}/mainwls cc)
endforeach(program)

# Install libraries
install(TARGETS ${astro_wl_libs} DESTINATION lib)

# Install binaraies
install(TARGETS ${astro_wl_all_targets} DESTINATION bin)

endif()
endif()

# -------- #

option(BUILD_ASTRO_GAL "Build ASTRO GAL library" ON)
message(STATUS "ASTRO GAL Build: ${BUILD_ASTRO_GAL}")

if(BUILD_ASTRO_GAL)

# Set path to astro Galaxies
set(astro_gal_path astro/galaxies)

# Build astro galaxies libs
set(astro_gal_libs gal)
foreach(library ${astro_gal_libs})
  build_lib(${library} ${astro_gal_path} cc h)
  add_dependencies(${library} ${misc_libs} ${mr_libs})
endforeach()
set(astrogal_alllibs ${misc_libs} ${astro_gal_libs})

# Find all astro gal targets
find_targets(astro_gal_targets ${astro_gal_path}/galmain cc)
set(astro_gal_all_targets ${astro_gal_targets})

# Build astro gal binaries
foreach(program ${astro_gal_targets})
  build_bin(${program} "${astrogal_alllibs}" ${astro_gal_path}/galmain cc)
endforeach(program)

# Install libraries
install(TARGETS ${astro_gal_libs} DESTINATION lib)

# Install binaraies
install(TARGETS ${astro_gal_all_targets} DESTINATION bin)

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
