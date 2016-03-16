add_executable(poisson1d poisson1d.cpp)
target_link_libraries(poisson1d ${LIBS})
add_test(
  NAME poisson1d
  COMMAND poisson1d
)

add_executable(poisson2dtri poisson2dtri.cpp)
target_link_libraries(poisson2dtri ${LIBS})
add_test(
  NAME poisson2dtri
  COMMAND poisson2dtri
)

add_executable(poisson2dquad poisson2dquad.cpp)
target_link_libraries(poisson2dquad ${LIBS})
add_test(
  NAME poisson2dquad
  COMMAND poisson2dquad
)
