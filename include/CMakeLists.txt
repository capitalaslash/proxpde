set(MINIFEM_HEADERS_LOCAL
  def.hpp
  array.hpp
  blockmatrix.hpp
  geo.hpp
  mesh.hpp
  reffe.hpp
  dof.hpp
  fe.hpp
  fespace.hpp
  curfe.hpp
  qr.hpp
  bc.hpp
  var.hpp
  assembly.hpp
  builder.hpp
  iomanager.hpp
  xdmf_traits.hpp
  timer.hpp
)

set(MINIFEM_HEADERS ${MINIFEM_HEADERS_LOCAL} PARENT_SCOPE)

install(FILES ${MINIFEM_HEADERS_LOCAL} DESTINATION include)

