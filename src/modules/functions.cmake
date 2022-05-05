# Locate pkgconfig package
function(find_pkg package module)
  if(${package}_LIBRARIES STREQUAL "" OR NOT DEFINED ${package}_LIBRARIES OR
     ${package}_LIBRARY_DIRS STREQUAL "" OR NOT DEFINED ${package}_LIBRARY_DIRS OR
     ${package}_INCLUDE_DIRS STREQUAL "" OR NOT DEFINED ${package}_INCLUDE_DIRS)
    pkg_check_modules(${package} REQUIRED ${module})
  endif()
  include_directories(${${package}_INCLUDE_DIRS})
  link_directories(${${package}_LIBRARY_DIRS})
  if(${package}_FOUND)
    message(STATUS "  include dirs: ${${package}_INCLUDE_DIRS}")
    message(STATUS "  lib dirs: ${${package}_LIBRARY_DIRS}")
    message(STATUS "  libs: ${${package}_LIBRARIES}")
  endif()
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
  install(FILES ${inc_${library}} DESTINATION include/cosmostat)
  install(TARGETS ${library} DESTINATION lib)
endfunction()

# Build list of libraries
function(build_lib_list library_list lib_path src_ext head_ext)
  foreach(library ${library_list})
    build_lib(${library} ${lib_path} ${src_ext} ${head_ext})
  endforeach()
endfunction()

# Build binary
function(build_bin program libs target_path ext)
  add_executable(${program} "${PROJECT_SOURCE_DIR}/${target_path}/${program}.${ext}")
  target_link_libraries(${program} ${libs})
endfunction()

# Build list of binaries
function(build_bin_list program_list libs target_path ext)
  foreach(program ${program_list})
    build_bin(${program} "${libs}" ${target_path} ${ext})
  endforeach()
endfunction()

# Build list of main directories
function(build_main_list main_list libs target_path ext)
  foreach(main ${main_list})
    find_targets(main_targets ${target_path}/${main} ${ext})
    build_bin_list("${main_targets}" "${libs}" ${target_path}/${main} ${ext})
    install(TARGETS ${main_targets} DESTINATION bin)
  endforeach()
endfunction()

# Build package
function(build_package main_list lib_list lib_deps_list package_path src_ext head_ext)
  build_lib_list("${lib_list}" ${package_path} ${src_ext} ${head_ext})
  set(package_all_libs ${lib_list} ${lib_deps_list})
  build_main_list("${main_list}" "${package_all_libs}" ${package_path} ${src_ext})
endfunction()

# Build Python binding
function(build_pybind program libs bind_path ext)
  add_library(${program} SHARED ${bind_path}/${program}.${ext})
  target_link_libraries(${program} "${libs}")
  if(APPLE)
    set_target_properties(${program} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
  else(APPLE)
    target_link_libraries(${program} ${Python_LIBRARIES})
  endif(APPLE)
  set_target_properties(${program} PROPERTIES SUFFIX .so)
  set_target_properties(${program} PROPERTIES PREFIX "")
endfunction()

# Build list of Python binding
function(build_pybind_list bind_list libs bind_path ext)
  foreach(binding ${bind_list})
    build_pybind(${binding} "${libs}" ${bind_path} ${ext})
  endforeach()
endfunction()
