macro(minifem_add_test test_name)
  add_executable(${test_name} ${test_name}.cpp)
  target_link_libraries(${test_name} ${LIBS})
  add_test(
    NAME ${test_name}
    COMMAND ${test_name}
  )
endmacro()

minifem_add_test(blockmatrix)

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

minifem_add_test(vectorfe1d)
minifem_add_test(assemblyvector)

minifem_add_test(stokes2dquad)
minifem_add_test(stokes2dcircle)
minifem_add_test(cavity)
minifem_add_test(cavitytime)
