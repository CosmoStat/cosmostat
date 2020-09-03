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
find_pkg(HEALPIX healpix_cxx)

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
# MRS #
# --- #

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

# -------- #
# Astro WL #
# -------- #

# Set path to astro weaklensing
set(astro_wl_path astro/weaklensing)

# Build astro weaklensing libs
set(astro_wl_libs divcurlwt wl wls)
foreach(library ${astro_wl_libs})
  build_lib(${library} ${astro_wl_path} cc h)
  add_dependencies(${library} ${mrs_libs})
endforeach()
set(astro_wl_mrs_libs ${mrs_libs} ${astro_wl_libs})

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

# --------------- #
# Python Bindings #
# --------------- #

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
