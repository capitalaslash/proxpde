macro(minifem_add_test test_name)
  add_executable(${test_name} ${test_name}.cpp)
  target_link_libraries(${test_name} MiniFEM::minifem)
  target_compile_options(${test_name} PRIVATE $<$<CXX_COMPILER_ID:GNU>:-Wall -Wpedantic>)
  add_test(
    NAME MiniFEM.${test_name}
    COMMAND ${test_name}
  )
endmacro()

include(CTest)
# enable_testing()

minifem_add_test(blockmatrix)
minifem_add_test(hdf5write)

minifem_add_test(poisson1d)
minifem_add_test(poisson1d_p2)
minifem_add_test(advection1d)

minifem_add_test(poisson2dtri)
minifem_add_test(poisson2dtri_p2)
minifem_add_test(advection2dtri)

minifem_add_test(poisson2dquad)
minifem_add_test(poisson2dquad_q2)
minifem_add_test(dr2dquad)
minifem_add_test(advection2dquad)

minifem_add_test(poisson3dtet)

minifem_add_test(poisson2dproj)

minifem_add_test(vectorfe1d)
minifem_add_test(assemblyvector)

minifem_add_test(stokes2dquad)
file(COPY ${PROJECT_SOURCE_DIR}/data/square_uns.msh
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
minifem_add_test(stokes2dtri)
minifem_add_test(stokes2dcircle)
minifem_add_test(ipcs)
minifem_add_test(cavity)
file(COPY ${PROJECT_SOURCE_DIR}/data/cavitytime.yaml
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
minifem_add_test(cavitytime)

minifem_add_test(gmres)

minifem_add_test(robin_bc)
minifem_add_test(projection_p1top2)
minifem_add_test(projection_p2top1)
minifem_add_test(projection_p1tort0)

