#include "fespace.hpp"

#include "qr.hpp"

// explicit instantiations
#ifdef MINIFEM_EXPLICIT_INSTANTIATION
template struct FESpace<Mesh<Line>, RefLineP0, GaussQR<Line, 1>>;
template struct FESpace<Mesh<Line>, RefLineP1, GaussQR<Line, 2>>;
template struct FESpace<Mesh<Line>, RefLineP2, GaussQR<Line, 3>>;
template struct FESpace<Mesh<Triangle>, RefTriangleP0, GaussQR<Triangle, 1>>;
template struct FESpace<Mesh<Triangle>, RefTriangleP1, GaussQR<Triangle, 3>>;
template struct FESpace<Mesh<Triangle>, RefTriangleP2, GaussQR<Triangle, 6>>;
template struct FESpace<Mesh<Triangle>, RefTriangleRT0, GaussQR<Triangle, 3>>;
template struct FESpace<Mesh<Quad>, RefQuadP0, GaussQR<Quad, 1>>;
template struct FESpace<Mesh<Quad>, RefQuadQ1, GaussQR<Quad, 4>>;
template struct FESpace<Mesh<Quad>, RefQuadQ2, GaussQR<Quad, 9>>;
template struct FESpace<Mesh<Quad>, RefQuadRT0, GaussQR<Quad, 4>>;
template struct FESpace<Mesh<Tetrahedron>, RefTetrahedronP0, GaussQR<Tetrahedron, 1>>;
template struct FESpace<Mesh<Tetrahedron>, RefTetrahedronP1, GaussQR<Tetrahedron, 4>>;
template struct FESpace<Mesh<Tetrahedron>, RefTetrahedronP2, GaussQR<Tetrahedron, 11>>;
template struct FESpace<Mesh<Hexahedron>, RefHexahedronP0, GaussQR<Hexahedron, 1>>;
template struct FESpace<Mesh<Hexahedron>, RefHexahedronQ1, GaussQR<Hexahedron, 8>>;
template struct FESpace<Mesh<Hexahedron>, RefHexahedronQ2, GaussQR<Hexahedron, 27>>;
#endif
