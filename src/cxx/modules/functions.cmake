# Locate pkgconfig package
function(find_pkg package module)
  if(${package}_LIBRARIES STREQUAL "" OR NOT DEFINED ${package}_LIBRARIES OR
     ${package}_LIBRARY_DIRS STREQUAL "" OR NOT DEFINED ${package}_LIBRARY_DIRS OR
     ${package}_INCLUDE_DIRS STREQUAL "" OR NOT DEFINED ${package}_INCLUDE_DIRS)
    pkg_check_modules(${package} REQUIRED ${module})
  else()
    message(STATUS "Use manually configured ${package}")
    message(STATUS "  includes: ${${package}_INCLUDE_DIRS}")
    message(STATUS "  libs: ${${package}_LIBRARY_DIRS}")
    message(STATUS "  flags: ${${package}_LIBRARIES}")
  endif()
  include_directories(${${package}_INCLUDE_DIRS})
  link_directories(${${package}_LIBRARY_DIRS})
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
  INSTALL(FILES ${inc_${library}} DESTINATION include/cosmostat)
endfunction()

# Build binary
function(build_bin program libs target_path ext)
  add_executable(${program} "${PROJECT_SOURCE_DIR}/${target_path}/${program}.${ext}")
  target_link_libraries(${program} ${CFITSIO_LIBRARIES} ${HEALPIX_LIBRARIES} ${fftw_libs} ${libs} ${sparse2d_libs})
endfunction()

# Build Python binding
function(build_pybind program bind_path ext)
  add_library(${program} SHARED ${bind_path}/${program}.${ext})
  target_link_libraries(${program} ${OpenMP_CXX_LIBRARIES} ${sparse2d_libs} ${mrs_libs} ${HEALPIX_LIBRARIES})
  if(APPLE)
    set_target_properties(${program} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
  else(APPLE)
    target_link_libraries(${program} ${Python_LIBRARIES})
  endif(APPLE)
  set_target_properties(${program} PROPERTIES SUFFIX .so)
  set_target_properties(${program} PROPERTIES PREFIX "")
endfunction()
