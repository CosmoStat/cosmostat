#=================#
# Build Sparse2D  #
#=================#

# Set Sparse2D Version
set(Sparse2DVersion v2.1.5_beta)

message(STATUS "SPARSE2D_SOURCE: ${SPARSE2D_SOURCE}")

# Download and build Sparse2D
ExternalProject_Add(sparse2d-git
  GIT_REPOSITORY https://github.com/sfarrens/Sparse2D.git
  GIT_TAG ${Sparse2DVersion}
  PREFIX sparse2d
  CMAKE_ARGS ${CMAKE_ARGS} -DBUILD_MSVST=ON -DUSE_FFTW=ON -DBUILD_NFFT=OFF
  # INSTALL_COMMAND make install
)

set(FFTW_LD_FLAGS "-L ${SPARSE2D_SOURCE}/sparse2d-git-build/module_build/fftw/lib -lfftw3f_omp -lfftw3_omp \
-lfftw3f -lfftw3 -lm")

set(sparse2d_lib_names mga2d sparse3d sparse2d sparse1d tools)
set(sparse2d_libs "")
foreach(library ${sparse2d_lib_names})
  include_directories(${SPARSE2D_SOURCE}/sparse2d-git/src/lib${library})
  find_library(${library} NAMES ${library} PATHS ${SPARSE2D_SOURCE}/sparse2d-git-build)
  list(APPEND sparse2d_libs ${${library}})
endforeach()
