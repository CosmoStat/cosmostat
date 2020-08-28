#=================#
# Build Sparse2D  #
#=================#

# Set Sparse2D Version
set(Sparse2DVersion v2.1.5_beta)

message(STATUS ${SPARSE2D_SOURCE})

# Download and build Sparse2D
ExternalProject_Add(sparse2d
  GIT_REPOSITORY https://github.com/sfarrens/Sparse2D.git
  GIT_TAG ${Sparse2DVersion}
  PREFIX sparse2d
  CMAKE_ARGS ${CMAKE_ARGS} -DBUILD_MSVST=ON -DUSE_FFTW=OFF -DBUILD_NFFT=OFF
  # INSTALL_COMMAND make install
)
