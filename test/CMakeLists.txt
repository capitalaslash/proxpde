add_executable(poisson1d poisson1d.cpp)
target_link_libraries(poisson1d ${LIBS})
add_test(
  NAME poisson1d
  COMMAND poisson1d
)

add_executable(poisson1d_p2 poisson1d_p2.cpp)
target_link_libraries(poisson1d_p2 ${LIBS})
add_test(
  NAME poisson1d_p2
  COMMAND poisson1d_p2
)

add_executable(poisson2dtri poisson2dtri.cpp)
target_link_libraries(poisson2dtri ${LIBS})
add_test(
  NAME poisson2dtri
  COMMAND poisson2dtri
)

add_executable(poisson2dtri_p2 poisson2dtri_p2.cpp)
target_link_libraries(poisson2dtri_p2 ${LIBS})
add_test(
  NAME poisson2dtri_p2
  COMMAND poisson2dtri_p2
)

add_executable(poisson2dquad poisson2dquad.cpp)
target_link_libraries(poisson2dquad ${LIBS})
add_test(
  NAME poisson2dquad
  COMMAND poisson2dquad
)

add_executable(poisson2dquad_q2 poisson2dquad_q2.cpp)
target_link_libraries(poisson2dquad_q2 ${LIBS})
add_test(
  NAME poisson2dquad_q2
  COMMAND poisson2dquad_q2
)

add_executable(dr2dquad dr2dquad.cpp)
target_link_libraries(dr2dquad ${LIBS})
add_test(
  NAME dr2dquad
  COMMAND dr2dquad
)

add_executable(stokes2dquad stokes2dquad.cpp)
target_link_libraries(stokes2dquad ${LIBS})
add_test(
  NAME stokes2dquad
  COMMAND stokes2dquad
)
