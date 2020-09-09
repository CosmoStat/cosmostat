#=================#
# Build Sparse2D  #
#=================#

# Set Sparse2D Version
set(Sparse2DVersion v2.1.5_beta3)

# Set Sparse2D source dir
set(SPARSE2D_SOURCE ${CMAKE_CURRENT_BINARY_DIR}/sparse2d)

option(BUILD_FFTW "Build FFTW libraries" OFF)
if(NOT FFTW3_FOUND)
  set(BUILD_FFTW ON)
endif()

# Download and build Sparse2D
ExternalProject_Add(sparse2d-git
  GIT_REPOSITORY https://github.com/sfarrens/Sparse2D.git
  GIT_TAG ${Sparse2DVersion}
  PREFIX sparse2d
  CMAKE_ARGS ${CMAKE_ARGS} -DBUILD_MSVST=ON -DUSE_FFTW=ON -DBUILD_FFTW=${BUILD_FFTW} -DBUILD_NFFT=OFF
-DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}
)

# List Sparse2D libs
set(sparse2d_libs mga2d sparse3d sparse2d sparse1d tools)

# Include external headers
include_directories(${CMAKE_INSTALL_PREFIX}/include/sparse2d)

# List FFTW libs
set(fftw_libs fftw3f_omp fftw3_omp fftw3f fftw3)

if(NOT FFTW3_FOUND)
  # Extract FFTW libraries from Sparse2D build
  ExternalProject_Add_Step(sparse2d-git fftw
    DEPENDEES install
    COMMAND cp -r ${SPARSE2D_SOURCE}/src/sparse2d-git-build/module_build/lib/. ${CMAKE_INSTALL_PREFIX}/lib
    COMMAND cp -r ${SPARSE2D_SOURCE}/src/sparse2d-git-build/module_build/include/. ${CMAKE_INSTALL_PREFIX}/include/sparse2d
    COMMENT "Extracting FFTW build"
  )
  # Add external libraries
  set(external_libs ${sparse2d_libs} ${fftw_libs})
else()
  # Add external libraries
  set(external_libs ${sparse2d_libs})
endif()

foreach(library ${external_libs})
  add_library(${library} STATIC IMPORTED)
  set_property(TARGET ${library} PROPERTY IMPORTED_LOCATION ${CMAKE_INSTALL_PREFIX}/lib/lib${library}.a)
  add_dependencies(${library} sparse2d-git)
endforeach()
