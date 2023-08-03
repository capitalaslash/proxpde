macro(proxpde_add_test name)
  # check arguments
  set(options EXCLUDE_FROM_ALL WILL_FAIL)
  set(one_value_args HEADER_ROOT)
  set(multi_value_args
      COMPILE_FLAGS
      COMPONENT_DEPENDENCIES
      DATA
      DEPENDENCIES
      HEADERS
      LABELS
      LINK_FLAGS
      LIBRARIES
      SOURCES
  )
  cmake_parse_arguments(${name} "${options}" "${one_value_args}" "${multi_value_args}" ${ARGN})

  # message(STATUS "name = ${name}")
  # message(STATUS "options = 1- ${${name}_EXCLUDE_FROM_ALL} 2- ${${name}_WILL_FAIL}")
  # message(STATUS "one_value-args = 1- ${${name}_HEADER_ROOT}")
  # message(STATUS "multi_value-args = 1- ${${name}_SOURCES} 2- ${${name}_HEADERS} 3- ${${name}_DATA} 4- ${${name}_DEPENDENCIES} 5- ${${name}_COMPONENT_DEPENDENCIES} 6- ${${name}_COMPILE_FLAGS} 7- ${${name}_LINK_FLAGS}")

  if(NOT ${name}_SOURCES)
    set(${name}_SOURCES ${name}.cpp)
  endif()
  # if(NOT ${name}_COMPILE_FLAGS) set(${name}_COMPILE_FLAGS "-Wall -Wpedantic") endif()

  add_executable(${name} ${${name}_SOURCES})
  set_target_properties(${name} PROPERTIES EXCLUDE_FROM_ALL TRUE)
  add_dependencies(build_tests ${name})
  # TODO: support COMPONENTS
  target_link_libraries(${name} PUBLIC ProXPDE::proxpde)
  # TODO: manage optional compile flags
  target_compile_options(${name}
    PRIVATE ${PROXPDE_COMPILE_FLAGS}
  )
  target_link_options(${name}
    PRIVATE ${PROXPDE_LINK_FLAGS}
  )
  if(${name}_LIBRARIES)
    target_link_libraries(${name} PUBLIC ${${name}_LIBRARIES})
  endif()
  add_test(
    NAME ProXPDE.${name}
    COMMAND ${name}
  )
  if(${name}_LABELS)
    set_property(TEST ProXPDE.${name} PROPERTY LABELS ${${name}_LABELS})
  endif()
  if(${${name}_WILL_FAIL})
    set_tests_properties(ProXPDE.${name} PROPERTIES WILL_FAIL TRUE)
  endif()

  foreach(datafile ${${name}_DATA})
    message(STATUS "copying ${datafile}")
    file(COPY ${PROJECT_SOURCE_DIR}/data/${datafile}
         DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
    )
  endforeach()

endmacro()
