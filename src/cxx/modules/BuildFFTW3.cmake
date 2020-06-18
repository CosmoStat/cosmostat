#=============#
# Build FFTW3 #
#=============#

set(FFTWVersion 3.3.8)

if(USE_FFTW)

  ExternalProject_Add(fftw3
      URL               http://www.fftw.org/fftw-${FFTWVersion}.tar.gz
      CONFIGURE_COMMAND ""
      BUILD_COMMAND     "${CMAKE_SOURCE_DIR}/modules/fftw.cmd"
      SOURCE_DIR        "${MODULE_BUILD_DIR}/fftw"
      INSTALL_COMMAND   ""
      BUILD_IN_SOURCE   1
      )

  set(FFTW_CXX_FLAGS "-DUSE_FFTW")
  set(FFTW_LD_FLAGS "-L ${MODULE_BUILD_DIR}/fftw/lib -lfftw3f_omp -lfftw3_omp \
-lfftw3f -lfftw3 -lm")

else(USE_FFTW)

  set(FFTW_CXX_FLAGS "")
  set(FFTW_LD_FLAGS "")

endif(USE_FFTW)

message(STATUS "FFTW3 Build: ${USE_FFTW}")
