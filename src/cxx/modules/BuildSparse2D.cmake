#=================#
# Build Sparse2D  #
#=================#

# Set Sparse2D Version
set(Sparse2DVersion v2.1.5_beta2)

# Download and build Sparse2D
ExternalProject_Add(sparse2d-git
  GIT_REPOSITORY https://github.com/sfarrens/Sparse2D.git
  GIT_TAG ${Sparse2DVersion}
  PREFIX sparse2d
  CMAKE_ARGS ${CMAKE_ARGS} -DBUILD_MSVST=ON -DUSE_FFTW=ON -DBUILD_NFFT=OFF
-DCMAKE_INSTALL_PREFIX:PATH=${SPARSE2D_SOURCE}
)

# Extract FFTW libraries from Sparse2D build
ExternalProject_Add_Step(sparse2d-git fftw
  DEPENDEES install
  COMMAND cp -r ${SPARSE2D_SOURCE}/src/sparse2d-git-build/module_build/lib/. ${SPARSE2D_SOURCE}/lib
  COMMAND cp -r ${SPARSE2D_SOURCE}/src/sparse2d-git-build/module_build/include/. ${SPARSE2D_SOURCE}/include/sparse2d
  COMMENT "Extracting FFTW build"
)

# List Sparse2D libs
set(sparse2d_libs mga2d sparse3d sparse2d sparse1d tools)

# List FFTW libs
set(fftw_libs fftw3f_omp fftw3_omp fftw3f fftw3)

# Include external headers
include_directories(${SPARSE2D_SOURCE}/include/sparse2d)

# Add external libraries
set(external_libs ${sparse2d_libs} ${fftw_libs})
foreach(library ${external_libs})
  add_library(${library} STATIC IMPORTED)
  set_property(TARGET ${library} PROPERTY IMPORTED_LOCATION ${SPARSE2D_SOURCE}/lib/lib${library}.a)
  add_dependencies(${library} sparse2d-git)
endforeach()
