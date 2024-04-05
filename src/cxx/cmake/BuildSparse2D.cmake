# =================================================#
# Build the Sparse2D dependencies for the project #
# =================================================#

# Set Sparse2D Version
set(sparse2dVersion 3.0.1)

# Download and build Sparse2D
ExternalProject_Add(sparse2d
    PREFIX sparse2d
    GIT_REPOSITORY https://github.com/CosmoStat/Sparse2D.git
    GIT_TAG v${sparse2dVersion}
    CONFIGURE_COMMAND cmake ../sparse2d
    -DCMAKE_INSTALL_PREFIX=${SPARSE2D_INSTALL_DIR}/..
    -DCMAKE_BUILD_TYPE=RELEASE
    -DPYBIND_INSTALL_PATH=${Python_SITELIB}
    BUILD_COMMAND make install -j8
    BUILD_IN_SOURCE 0
)
