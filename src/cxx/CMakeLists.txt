##################
# COSMOSTAT 2020 #
##################

cmake_minimum_required(VERSION 3.0.0)
include(ExternalProject)
include(FindPkgConfig)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/modules/")

project(cosmostat)
message(STATUS "Installing to: ${CMAKE_INSTALL_PREFIX}")

option(SPARSE3D "Build Sparse3D library" OFF)
option(BUILD_MSVST "Build MSVST library" OFF)
option(USE_FFTW "Use FFTW for Fourier transformation" ON)
option(COMPILE_NFFT "Get and compile NFFT" OFF)
option(MRS "Build MRS library" ON)
option(MRSP "Build MRS library" ON)
enable_language (Fortran OPTIONAL)
set(COMPILE_FORTRAN_BIN 1)
set(COMPILE_BINDING 1)

  # Set build type
  if(CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE "Release")
  endif()
  message(STATUS "CMake Build Type: ${CMAKE_BUILD_TYPE}")

  # Set the default installation path
  if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_BINARY_DIR}/../cosmostat"
        CACHE PATH "Setting new default path" FORCE)
  endif()
  message(STATUS "Installing to: ${CMAKE_INSTALL_PREFIX}")

  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
    set(OpenMP_INCLUDE_PATH "/usr/local/include")
    set(OpenMP_LIB_PATH "/usr/local/lib")
    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -lomp")
    set(OpenMP_CXX_LIB_NAMES "omp")
    set(OpenMP_CXX_LIBRARIES "${OpenMP_LIB_PATH}/libomp.dylib")
    include_directories(${OpenMP_INCLUDE_PATH})
    link_directories(${OpenMP_LIB_PATH})
  else()
    find_package(OpenMP REQUIRED)
  endif()

  # Locate CFITSIO using pkg-config or use command line arguments to configure
  # CFITSIO
  if(CFITSIO_LIBRARIES STREQUAL "" OR NOT DEFINED CFITSIO_LIBRARIES OR
     CFITSIO_LIBRARY_DIRS STREQUAL "" OR NOT DEFINED CFITSIO_LIBRARY_DIRS OR
     CFITSIO_INCLUDE_DIRS STREQUAL "" OR NOT DEFINED CFITSIO_INCLUDE_DIRS)
    pkg_check_modules(CFITSIO REQUIRED cfitsio)
  else()
    message(STATUS "Use manually configured cfitsio")
    message(STATUS "  includes: ${CFITSIO_INCLUDE_DIRS}")
    message(STATUS "  libs: ${CFITSIO_LIBRARY_DIRS}")
    message(STATUS "  flags: ${CFITSIO_LIBRARIES}")
  endif()
  include_directories(${CFITSIO_INCLUDE_DIRS})
  link_directories(${CFITSIO_LIBRARY_DIRS})

  # Build location for external modules
  set(MODULE_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/module_build/")

  # Optional build of FFTW3
  include(BuildFFTW3)

  # Optional build of NFFT
  include(BuildNFFT)

  # Set compilation flags Mac OSX or any other system
  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")

    set(CMAKE_CXX_FLAGS_RELEASE "-DNO_DISP_IO -g0 -O2 -fPIC \
-fomit-frame-pointer ${OpenMP_CXX_FLAGS} -Wno-write-strings -DNDEBUG \
${FFTW_FLAGS}")
    set(CMAKE_CXX_FLAGS_DEBUG "-DNO_DISP_IO -g -fPIC ${OpenMP_CXX_FLAGS} \
-Wno-write-strings ${FFTW_CXX_FLAGS}")

  else()

    set(CMAKE_CXX_FLAGS_RELEASE "-DNO_DISP_IO -g0 -O2 -fPIC \
-fomit-frame-pointer ${OpenMP_CXX_FLAGS} -Wno-write-strings -DNDEBUG \
${FFTW_FLAGS}"
    )
    set(CMAKE_CXX_FLAGS_DEBUG "-DNO_DISP_IO -g -fPIC ${OpenMP_CXX_FLAGS} \
-Wno-write-strings ${FFTW_CXX_FLAGS}"
    )

    message(STATUS "Using CXX Flags for ${CMAKE_SYSTEM}")

  endif()


if(MRS)
# to BE REPLACED by a module
pkg_check_modules(HEALPIX healpix_cxx>=3.60)
# find_package(healpix_cxx)
if (HEALPIX_FOUND)
  include_directories(${HEALPIX_INCLUDE_DIRS})
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
  set(DO_NOT_COMPILE_HEALPIX_BIN 1)
  set(DO_NOT_COMPILE_HEALPIX_BIN_POL 1)
endif(HEALPIX_FOUND)
endif(MRS)

  # Build tools library
  FILE(GLOB src_lib1  "${PROJECT_SOURCE_DIR}/lib/libtools/*.cc")
  include_directories("${PROJECT_SOURCE_DIR}/lib/libtools")
  add_library(tools STATIC ${src_lib1})
  target_link_libraries(tools ${CFITSIO_LIBRARIES} ${FFTW_LD_FLAGS})
  if(USE_FFTW)
    add_dependencies(tools fftw3)
  endif(USE_FFTW)

  # Build sparse1d library
  FILE(GLOB src_lib2 "${PROJECT_SOURCE_DIR}/lib/libsparse1d/*.cc")
  include_directories("${PROJECT_SOURCE_DIR}/lib/libsparse1d")
  add_library(sparse1d STATIC ${src_lib2})
  target_link_libraries(sparse1d ${CFITSIO_LIBRARIES} ${FFTW_LD_FLAGS})
  if(USE_FFTW)
    add_dependencies(sparse1d fftw3)
  endif(USE_FFTW)

  # Build sparse2d library
  FILE(GLOB src_lib2 "${PROJECT_SOURCE_DIR}/lib/libsparse2d/*.cc")
  include_directories("${PROJECT_SOURCE_DIR}/lib/libsparse2d")
  add_library(sparse2d STATIC ${src_lib2})
  target_link_libraries(sparse2d ${CFITSIO_LIBRARIES} ${FFTW_LD_FLAGS})
  if(USE_FFTW)
    add_dependencies(sparse2d fftw3)
  endif(USE_FFTW)

  if(SPARSE3D)
    # Build sparse3d library
    FILE(GLOB src_lib2 "${PROJECT_SOURCE_DIR}/lib/libsparse3d/*.cc")
    include_directories("${PROJECT_SOURCE_DIR}/lib/libsparse3d")
    add_library(sparse3d STATIC ${src_lib2})
    target_link_libraries(sparse3d ${CFITSIO_LIBRARIES} ${FFTW_LD_FLAGS})
    if(USE_FFTW)
      add_dependencies(sparse3d fftw3)
    endif(USE_FFTW)
  endif(SPARSE3D)
  message(STATUS "Sparse3D Build: ${SPARSE3D}")

  # Build mga2d library
  FILE(GLOB src_mgalib2 "${PROJECT_SOURCE_DIR}/lib/libmga2d/*.cc")
  include_directories("${PROJECT_SOURCE_DIR}/lib/libmga2d")
  add_library(mga2d STATIC ${src_mgalib2})
  target_link_libraries(mga2d ${CFITSIO_LIBRARIES} ${FFTW_LD_FLAGS})
  if(USE_FFTW)
    add_dependencies(mga2d fftw3)
  endif(USE_FFTW)

##################

  # Compile and link executables
  set(BINMR1D im1d_convert im1d_info im1d_stf
      im1d_tfreq mr1d_detect mr1d_trans im1d_der_loggamma im1d_period im1d_tend
      mr1d_filter)
  set(BINMR2D mr_transform mr_recons mr_filter mr_deconv mr_info cur_contrast
      cur_deconv cur_filter cur_stat cur_trans im_isospec)
  if(SPARSE3D)
    set(BINMR3D mr3d_trans mr3d_filter mr3d_recons mr3d_stat mr2d1d_trans)
  endif(SPARSE3D)

  foreach(program ${BINMR1D})
    add_executable(${program} ${PROJECT_SOURCE_DIR}/main1d/${program}.cc)
    target_link_libraries(${program} sparse1d tools)
  endforeach(program)

  foreach(program ${BINMR2D})
    add_executable(${program} ${PROJECT_SOURCE_DIR}/main2d/${program}.cc)
    target_link_libraries(${program} mga2d sparse2d sparse1d tools)
  endforeach(program)

  if(SPARSE3D)
    foreach(program ${BINMR3D})
    	add_executable(${program} ${PROJECT_SOURCE_DIR}/main3d/${program}.cc)
    	target_link_libraries(${program} mga2d sparse3d sparse2d sparse1d tools)
    endforeach(program)
  endif(SPARSE3D)

################## MSVST ##################

  if(BUILD_MSVST)
    FILE(GLOB msvst_src_lib "${PROJECT_SOURCE_DIR}/msvst/libmsvst/*.cc")
    include_directories("${PROJECT_SOURCE_DIR}/msvst/libmsvst")
    add_library(msvst STATIC ${msvst_src_lib} )
  endif(BUILD_MSVST)
  message(STATUS "MSVST Build: ${BUILD_MSVST}")


  if(BUILD_MSVST)
    set(MSVST msvst_uwt2d msvst_iwt2d msvst_2d1d msvst_iwt2d_coupled)
    foreach(program ${MSVST})
      add_executable(${program} ${PROJECT_SOURCE_DIR}/msvst/${program}.cc)
      target_link_libraries(${program} msvst tools)
    endforeach(program)
    FILE(GLOB inc_lib "${PROJECT_SOURCE_DIR}/src/libmsvst/*.h")
    INSTALL(FILES ${inc_lib} DESTINATION include/cosmostat)
    INSTALL(TARGETS ${MSVST} DESTINATION bin)
  endif(BUILD_MSVST)

################## MRS ##################
if(MRS)
  # Build mrs library
  FILE(GLOB src_lib  "${PROJECT_SOURCE_DIR}/mrs/libmrs/*.cc")
  include_directories("${PROJECT_SOURCE_DIR}/mrs/libmrs")
  add_library(mrs STATIC ${src_lib})
  target_link_libraries(mrs ${CFITSIO_LIBRARIES} ${FFTW_LD_FLAGS})
  if(USE_FFTW)
    add_dependencies(mrs fftw3)
  endif(USE_FFTW)

  # set(MRSBIN mrs_info)
  # mrs_planck_lgmca_inv  mrs_gmca
  set(MRSBIN mrs_info  mrs_crosspowspec mrs_alm_inpainting  mrs_almrec mrs_alm_rotate  mrs_almtrans mrs_clean_cmb mrs_deriv mrs_getcmb mrs_mca mrs_inpainting mrs_iso_inpainting mrs_iter_wiener mrs_matmask mrs_owttrans mrs_powspec mrs_pwttrans mrs_resize mrs_smooth mrs_uwttrans mrs_wiener mrs_wt_filtering)

  # Build and link the module
  foreach(program ${MRSBIN})
    add_executable(${program} ${PROJECT_SOURCE_DIR}/mrs/mainmrs/${program}.cc)
    target_link_libraries(${program} mrs mga2d sparse2d sparse1d tools ${CFITSIO_LIBRARIES} ${HEALPIX_LINK_LIBRARIES})
   endforeach(program)
endif(MRS)

################## MRSP ##################

if(MRSP)
   FILE(GLOB src_lib  "${PROJECT_SOURCE_DIR}/mrs/libmrsp/*.cc")
    include_directories("${PROJECT_SOURCE_DIR}/mrs/libmrsp")
    add_library(mrsp STATIC ${src_lib})
    target_link_libraries(mrsp ${CFITSIO_LIBRARIES} ${FFTW_LD_FLAGS})
    if(USE_FFTW)
      add_dependencies(mrs fftw3)
    endif(USE_FFTW)

   if(COMPILE_FORTRAN_BIN MATCHES "1")
    if(DEFINED FORTRAN_FOUND)
     FILE(GLOB fortran_pol_lib  "${PROJECT_SOURCE_DIR}/libmrsp/f_wigner3j/*.f")
     add_library(wigner3j  ${fortran_pol_lib})
     add_executable(mrsp_matmask  ${PROJECT_SOURCE_DIR}/mrs/pol/mrsp_matmask.cc)
     target_link_libraries(mrsp_matmask libmrsp libmrs  wigner3j "-lgfortran")
    endif(DEFINED FORTRAN_FOUND)
   endif(COMPILE_FORTRAN_BIN MATCHES "1")

set(MRSPBIN mrsp_almtrans mrsp_almrec mrsp_alm2spec mrsp_spec mrsp_teb2tqu mrsp_tqu2teb mrsp_resize mrsp_smooth mrsp_wiener mrsp_alm_inpainting mrsp_alm_inpainting_spectra mrsp_eb_to_qu mrsp_matmask mrsp_qu_to_eb)
 
  # set(MRSPBIN mrsp_almtrans mrsp_spec mrsp_teb2tqu mrsp_tqu2teb mrsp_wiener)
  # Build and link the module
  foreach(program ${MRSPBIN})
    add_executable(${program} ${PROJECT_SOURCE_DIR}/mrs/mainmrsp/${program}.cc)
    target_link_libraries(${program} mrsp mrs mga2d sparse2d sparse1d tools ${CFITSIO_LIBRARIES} ${HEALPIX_LINK_LIBRARIES})
    add_dependencies(${program} mrs)
    if(APPLE)
      set_target_properties(${program} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
    endif(APPLE)
   endforeach(program)
endif(MRSP)

################## BINDING ##################

if(COMPILE_BINDING MATCHES "1")
  # Find default python libraries and interpreter
  find_package(PythonInterp REQUIRED)
  if(NOT DEFINED PYTHON_LIBRARIES OR NOT DEFINED PYTHON_INCLUDE_DIRS)
    find_package(PythonLibs REQUIRED)
  else()
    message(STATUS "Using Python Lib: ${PYTHON_LIBRARIES}")
    message(STATUS "Using Python Inc: ${PYTHON_INCLUDE_DIRS}")
  endif()
  include_directories(${PYTHON_INCLUDE_DIRS})
  link_directories(${PYTHON_LIBRARY_DIRS})

  # Include Pybind11 headers
  message(STATUS "Using Pybind11 Inc: ${PYBIND11_INCLUDE_DIR}")
  include_directories(${PYBIND11_INCLUDE_DIR})
  include_directories(${PYTHON_INCLUDE_DIRS})

 set(CSBinding pymrs)
 foreach(program ${CSBinding})
    add_library(${program} SHARED python_binding/${program}.cpp)
    if(APPLE)
      target_link_libraries(${program} ${OpenMP_CXX_LIBRARIES} mrs mga2d sparse2d sparse1d tools ${CFITSIO_LIBRARIES})
    else(APPLE)
      target_link_libraries(${program} ${OpenMP_CXX_LIBRARIES} mrs mga2d sparse2d sparse1d tools ${CFITSIO_LIBRARIES} ${PYTHON_LIBRARIES})
    endif(APPLE)

 #   add_dependencies(${program} cfitsio)

    set_target_properties(${program} PROPERTIES SUFFIX .so)
    set_target_properties(${program} PROPERTIES PREFIX "")

    if(APPLE)
      set_target_properties(${program} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
    endif(APPLE)

  endforeach(program)
endif(COMPILE_BINDING MATCHES "1")

################## INSTALL ##################
# INSTALLATION

  ##### Install headers
  FILE(GLOB inc_lib  "${PROJECT_SOURCE_DIR}/lib/libsparse1d/*.h")
  INSTALL(FILES ${inc_lib} DESTINATION include)
  FILE(GLOB inc_lib  "${PROJECT_SOURCE_DIR}/lib/libsparse2d/*.h")
  INSTALL(FILES ${inc_lib} DESTINATION include)
  FILE(GLOB inc_lib  "${PROJECT_SOURCE_DIR}/lib/libmga2d/*.h")
  INSTALL(FILES ${inc_lib} DESTINATION include)
  FILE(GLOB inc_lib  "${PROJECT_SOURCE_DIR}/lib/libtools/*.h")
  INSTALL(FILES ${inc_lib} DESTINATION include)
  if(SPARSE3D)
    FILE(GLOB inc_lib  "${PROJECT_SOURCE_DIR}/lib/libsparse3d/*.h")
    INSTALL(FILES ${inc_lib} DESTINATION include)
  endif(SPARSE3D)

  #####  Install library
  INSTALL(TARGETS sparse1d sparse2d mga2d tools DESTINATION lib)
  if(SPARSE3D)
    INSTALL(TARGETS sparse3d DESTINATION lib)
  endif(SPARSE3D)
  if(BUILD_MSVST)
 	 INSTALL(TARGETS msvst DESTINATION lib)
  endif(BUILD_MSVST)
  if(MRS)
 	 INSTALL(TARGETS mrs DESTINATION lib)
  endif(MRS)
  if(MRSP)
 	 INSTALL(TARGETS mrsp DESTINATION lib)
  endif(MRSP)
  if(MRSP)
 	 INSTALL(TARGETS mrsp DESTINATION lib)
  endif(MRSP)
  if(COMPILE_BINDING MATCHES "1")
 	 INSTALL(TARGETS ${CSBinding} DESTINATION lib)
  endif(COMPILE_BINDING MATCHES "1")

  #####  install  executables
  INSTALL(TARGETS ${BINMR1D} ${BINMR2D} DESTINATION bin)
  if(SPARSE3D)
 	 INSTALL(TARGETS ${BINMR3D} DESTINATION bin)
  endif(SPARSE3D)
  if(BUILD_MSVST)
 	 INSTALL(TARGETS ${MSVST} DESTINATION bin)
  endif(BUILD_MSVST)

if(MRS)
  INSTALL(TARGETS ${MRSBIN} DESTINATION bin)
endif(MRS)

if(MRSP)
 INSTALL(TARGETS ${MRSPBIN} DESTINATION bin)
 if(COMPILE_FORTRAN_BIN MATCHES "1")
    if(DEFINED FORTRAN_FOUND)
	 INSTALL(TARGETS mrsp_matmask DESTINATION bin)
       endif(DEFINED FORTRAN_FOUND)
 endif(COMPILE_FORTRAN_BIN MATCHES "1")
endif(MRSP)






