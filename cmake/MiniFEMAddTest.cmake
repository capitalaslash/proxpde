macro(minifem_add_test name)
  # check arguments
  set(options EXCLUDE_FROM_ALL)
  set(one_value_args HEADER_ROOT)
  set(multi_value_args SOURCES HEADERS DATA DEPENDENCIES COMPONENT_DEPENDENCIES COMPILE_FLAGS LINK_FLAGS)
  cmake_parse_arguments(${name} "${options}" "${one_value_args}" "${multi_value_args}" ${ARGN})

  # message(STATUS "name = ${name}")
  # message(STATUS "options = 1- ${${name}_EXCLUDE_FROM_ALL}")
  # message(STATUS "one_value-args = 1- ${${name}_HEADER_ROOT}")
  # message(STATUS "multi_value-args = 1- ${${name}_SOURCES} 2- ${${name}_HEADERS} 3- ${${name}_DATA} 4- ${${name}_DEPENDENCIES} 5- ${${name}_COMPONENT_DEPENDENCIES} 6- ${${name}_COMPILE_FLAGS} 7- ${${name}_LINK_FLAGS}")

  if(NOT ${name}_SOURCES)
    set(${name}_SOURCES ${name}.cpp)
  endif()
  # if(NOT ${name}_COMPILE_FLAGS)
  #   set(${name}_COMPILE_FLAGS "-Wall -Wpedantic")
  # endif()

  # message(STATUS "name = ${name}")
  # message(STATUS "options = 1- ${${name}_EXCLUDE_FROM_ALL}")
  # message(STATUS "one_value-args = 1- ${${name}_HEADER_ROOT}")
  # message(STATUS "multi_value-args = 1- ${${name}_SOURCES} 2- ${${name}_HEADERS} 3- ${${name}_DATA} 4- ${${name}_DEPENDENCIES} 5- ${${name}_COMPONENT_DEPENDENCIES} 6- ${${name}_COMPILE_FLAGS} 7- ${${name}_LINK_FLAGS}")

  add_executable(${name} ${${name}_SOURCES})
  # TODO: support COMPONENTS
  target_link_libraries(${name} MiniFEM::minifem)
  # TODO: manage optional compile flags
  target_compile_options(${name}
    PRIVATE ${MINIFEM_COMPILE_FLAGS}
  )
  target_link_options(${name}
    PRIVATE ${MINIFEM_LINK_FLAGS}
  )
  add_test(
    NAME MiniFEM.${name}
    COMMAND ${name}
  )

  foreach(datafile ${${name}_DATA})
    message(STATUS "copying ${datafile}")
    file(COPY ${PROJECT_SOURCE_DIR}/data/${datafile}
         DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  endforeach()

endmacro()

