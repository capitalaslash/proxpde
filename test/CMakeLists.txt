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
  target_compile_options(${name} PRIVATE $<$<CXX_COMPILER_ID:GNU>:-Wall -Wpedantic>)
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

# basic tests
minifem_add_test(blockmatrix)
minifem_add_test(hdf5write)
minifem_add_test(qr)
minifem_add_test(normal)
minifem_add_test(inside_test)
minifem_add_test(reffe_test)
minifem_add_test(solver_test)
minifem_add_test(vectorfe1d)
minifem_add_test(assemblyvector)
minifem_add_test(timeder)

# 1d
minifem_add_test(poisson1d)
minifem_add_test(poisson1d_p2)
minifem_add_test(advection1d DATA advection1d.yaml)

# 2d tri
minifem_add_test(poisson2dtri)
minifem_add_test(poisson2dtri_p2)
minifem_add_test(poisson2dmixedtri)
minifem_add_test(advection2dtri)
minifem_add_test(dual)

# 2d quad
minifem_add_test(poisson2dquad)
minifem_add_test(poisson2dquad_q2)
minifem_add_test(poisson2dmixedquad)
minifem_add_test(dr2dquad)
minifem_add_test(advection2dquad DATA square_q.msh)

# 3d tet
minifem_add_test(poisson3dtet DATA poisson3dtet.yaml cube_uns.msh)
minifem_add_test(poisson3dtet_p2)

# 3d hex
minifem_add_test(poisson3dhex)
minifem_add_test(poisson3dhex_q2)

minifem_add_test(poisson2dproj)
minifem_add_test(neumann)

# 2d stokes
minifem_add_test(stokes2dquad)
minifem_add_test(stokes2dtri DATA square_uns.msh)
minifem_add_test(stokes2dcircle)
minifem_add_test(cavity)
minifem_add_test(cavitytime DATA cavitytime.yaml)
minifem_add_test(nipcs)
minifem_add_test(ipcs)

# 3d stokes
minifem_add_test(stokes3dtet)

# bc
minifem_add_test(robin_bc)
minifem_add_test(fourier_nb)
minifem_add_test(projection_p1top2)
minifem_add_test(projection_p2top1)
minifem_add_test(projection_p1tort0)

# io
minifem_add_test(io)
minifem_add_test(iotime)
