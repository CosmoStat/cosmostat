# Locate CFITSIO using pkg-config or use command line arguments to configure
# CFITSIO
function(find_cfitsio)
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
endfunction()

# Extract target names from source files
function(find_targets targets target_path ext)
  file(GLOB src_targets "${PROJECT_SOURCE_DIR}/${target_path}/*.${ext}")
  list(TRANSFORM src_targets REPLACE "${PROJECT_SOURCE_DIR}/${target_path}/" "")
  list(TRANSFORM src_targets REPLACE ".${ext}" "")
  set(${targets} "${src_targets}" PARENT_SCOPE)
endfunction()

# Build library
function(build_lib library lib_path src_ext head_ext)
  file(GLOB src_${library} "${PROJECT_SOURCE_DIR}/${lib_path}/lib${library}/*.${src_ext}")
  file(GLOB inc_${library} "${PROJECT_SOURCE_DIR}/${lib_path}/lib${library}/*.${head_ext}")
  include_directories("${PROJECT_SOURCE_DIR}/${lib_path}/lib${library}")
  add_library(${library} STATIC ${src_${library}})
  add_dependencies(${library} sparse2d-git)
  INSTALL(FILES ${inc_${library}} DESTINATION include/sparse2d)
endfunction()

# Build binary
function(build_bin program libs target_path ext)
  add_executable(${program} "${PROJECT_SOURCE_DIR}/${target_path}/${program}.${ext}")
  target_link_libraries(${program} ${CFITSIO_LIBRARIES} ${fftw_libs} ${libs} ${sparse2d_libs})
  add_dependencies(${program} sparse2d-git)
endfunction()
