#include "qr.hpp"

namespace proxpde
{

// GaussQR =============================================================================

// 0D ----------------------------------------------------------------------------------
template <>
GaussQR<NullElem, 0U>::Weights_T const GaussQR<NullElem, 0U>::weight = {};
template <>
std::array<GaussQR<NullElem, 0U>::Vec_T, 0U> const GaussQR<NullElem, 0U>::node = {};
template <>
short_T const GaussQR<NullElem, 0>::bestPt = 0;

// 1D
// - Line ------------------------------------------------------------------------------
template <>
GaussQR<Line, 1U>::Weights_T const GaussQR<Line, 1U>::weight =
    GaussQR<Line, 1>::Weights_T::Constant(2.L);
template <>
std::array<GaussQR<Line, 1U>::Vec_T, 1> const GaussQR<Line, 1U>::node = {
    {GaussQR<Line, 1>::Vec_T::Constant(0.0L)}};
template <>
short_T const GaussQR<Line, 1U>::bestPt = 0;

static GaussQR<NullElem, 0U>::Real_T constexpr sqrt13rd = 0.5773502691896258L;
template <>
GaussQR<Line, 2U>::Weights_T const GaussQR<Line, 2U>::weight = {1.L, 1.L};
template <>
std::array<GaussQR<Line, 2U>::Vec_T, 2U> const GaussQR<Line, 2U>::node = {{
    GaussQR<Line, 2U>::Vec_T::Constant(-sqrt13rd),
    GaussQR<Line, 2U>::Vec_T::Constant(sqrt13rd),
}};
template <>
short_T const GaussQR<Line, 2U>::bestPt = 0;

static GaussQR<NullElem, 0>::Real_T constexpr sqrt35th = 0.774596669241483L;
template <>
GaussQR<Line, 3U>::Weights_T const GaussQR<Line, 3U>::weight = {
    5.L / 9, 8.L / 9, 5.L / 9};
template <>
std::array<GaussQR<Line, 3U>::Vec_T, 3U> const GaussQR<Line, 3U>::node = {
    {GaussQR<Line, 3U>::Vec_T::Constant(-sqrt35th),
     GaussQR<Line, 3U>::Vec_T::Constant(0.0L),
     GaussQR<Line, 3U>::Vec_T::Constant(sqrt35th)}};
template <>
short_T const GaussQR<Line, 3U>::bestPt = 1;

// sqrt((3.L - 2.L * sqrt(1.2L)) / 7)
static GaussQR<NullElem, 0>::Real_T constexpr sqrt37thm = 0.3399810435848563L;
// sqrt((3.L + 2.L * sqrt(1.2L)) / 7)
static GaussQR<NullElem, 0>::Real_T constexpr sqrt37thp = 0.8611363115940526L;
template <>
GaussQR<Line, 4U>::Weights_T const GaussQR<Line, 4U>::weight = {
    (18.L - sqrt(30.L)) / 36,
    (18.L + sqrt(30.L)) / 36,
    (18.L + sqrt(30.L)) / 36,
    (18.L - sqrt(30.L)) / 36};
template <>
std::array<GaussQR<Line, 4U>::Vec_T, 4U> const GaussQR<Line, 4U>::node = {{
    GaussQR<Line, 4U>::Vec_T::Constant(-sqrt37thp),
    GaussQR<Line, 4U>::Vec_T::Constant(-sqrt37thm),
    GaussQR<Line, 4U>::Vec_T::Constant(sqrt37thm),
    GaussQR<Line, 4U>::Vec_T::Constant(sqrt37thp),
}};
template <>
short_T const GaussQR<Line, 4U>::bestPt = 1;

// 2D
// - Triangle --------------------------------------------------------------------------
template <>
GaussQR<Triangle, 1U>::Weights_T const GaussQR<Triangle, 1U>::weight =
    GaussQR<Triangle, 1U>::Weights_T::Constant(0.5L);
template <>
std::array<GaussQR<Triangle, 1U>::Vec_T, 1U> const GaussQR<Triangle, 1U>::node = {
    {GaussQR<Triangle, 1U>::Vec_T(1.L / 3, 1.L / 3)}};
template <>
short_T const GaussQR<Triangle, 1U>::bestPt = 0;

template <>
GaussQR<Triangle, 3U>::Weights_T const GaussQR<Triangle, 3U>::weight =
    GaussQR<Triangle, 3U>::Weights_T::Constant(1.L / 6);
template <>
std::array<GaussQR<Triangle, 3U>::Vec_T, 3U> const GaussQR<Triangle, 3U>::node = {
    {GaussQR<Triangle, 3U>::Vec_T(1. / 6, 1. / 6),
     GaussQR<Triangle, 3U>::Vec_T(2. / 3, 1. / 6),
     GaussQR<Triangle, 3U>::Vec_T(1. / 6, 2. / 3)}};
// this node positioning seems to give smaller errors in some cases but it is
// not ideal since they coincide with the location of the dofs of the P2 element
// template<> std::array<GaussQR<Triangle,3>::Vec_T,3> const GaussQR<Triangle,3>::node =
// {{
//     GaussQR<Triangle, 3U>::Vec_T(0.5L, 0.0L),
//     GaussQR<Triangle, 3U>::Vec_T(0.5L, 0.5L),
//     GaussQR<Triangle, 3U>::Vec_T(0.0L, 0.5L)
// }};
template <>
short_T const GaussQR<Triangle, 3U>::bestPt = 0;

template <>
GaussQR<Triangle, 4U>::Weights_T const GaussQR<Triangle, 4U>::weight =
    (GaussQR<Triangle, 4U>::Weights_T() << 25.L, 25.L, 25.L, -27.L).finished() / 96.L;
template <>
std::array<GaussQR<Triangle, 4U>::Vec_T, 4U> const GaussQR<Triangle, 4U>::node = {
    {GaussQR<Triangle, 4U>::Vec_T(0.2L, 0.2L),
     GaussQR<Triangle, 4U>::Vec_T(0.6L, 0.2L),
     GaussQR<Triangle, 4U>::Vec_T(0.2L, 0.6L),
     GaussQR<Triangle, 4U>::Vec_T(1.L / 3, 1.L / 3)}};
template <>
short_T const GaussQR<Triangle, 4U>::bestPt = 3;

// clang-format off
template <>
GaussQR<Triangle, 6U>::Weights_T const GaussQR<Triangle, 6U>::weight =
    .5 * (GaussQR<Triangle, 6U>::Weights_T() <<
        0.109951743655322L, 0.109951743655322L, 0.109951743655322L,
        0.223381589678011L, 0.223381589678011L, 0.223381589678011L).finished();
// clang-format on
template <>
std::array<GaussQR<Triangle, 6U>::Vec_T, 6U> const GaussQR<Triangle, 6U>::node = {{
    GaussQR<Triangle, 6U>::Vec_T(0.091576213509771L, 0.091576213509771L),
    GaussQR<Triangle, 6U>::Vec_T(0.816847572980459L, 0.091576213509771L),
    GaussQR<Triangle, 6U>::Vec_T(0.091576213509771L, 0.816847572980459L),
    GaussQR<Triangle, 6U>::Vec_T(0.445948490915965L, 0.445948490915965L),
    GaussQR<Triangle, 6U>::Vec_T(0.108103018168070L, 0.445948490915965L),
    GaussQR<Triangle, 6U>::Vec_T(0.445948490915965L, 0.108103018168070L),
}};
template <>
short_T const GaussQR<Triangle, 6U>::bestPt = 3;

// clang-format off
template <>
GaussQR<Triangle, 7U>::Weights_T const GaussQR<Triangle, 7U>::weight =
    .5 * (GaussQR<Triangle, 7U>::Weights_T() <<
        0.22500000000000L,
        0.13239415278851L, 0.13239415278851L, 0.13239415278851L,
        0.12593918054483L, 0.12593918054483L, 0.12593918054483L).finished();
// clang-format on
template <>
std::array<GaussQR<Triangle, 7U>::Vec_T, 7U> const GaussQR<Triangle, 7U>::node = {{
    GaussQR<Triangle, 7U>::Vec_T(1.L / 3, 1.L / 3),
    GaussQR<Triangle, 7U>::Vec_T(0.47014206410511L, 0.47014206410511L),
    GaussQR<Triangle, 7U>::Vec_T(0.47014206410511L, 0.05971587178977L),
    GaussQR<Triangle, 7U>::Vec_T(0.05971587178977L, 0.47014206410511L),
    GaussQR<Triangle, 7U>::Vec_T(0.10128650732346L, 0.10128650732346L),
    GaussQR<Triangle, 7U>::Vec_T(0.10128650732346L, 0.79742698535309L),
    GaussQR<Triangle, 7U>::Vec_T(0.79742698535309L, 0.10128650732346L),
}};
template <>
short_T const GaussQR<Triangle, 7U>::bestPt = 0;

// - Quad ------------------------------------------------------------------------------
template <>
GaussQR<Quad, 1U>::Weights_T const GaussQR<Quad, 1U>::weight =
    GaussQR<Quad, 1U>::Weights_T::Constant(4.L);
template <>
std::array<GaussQR<Quad, 1U>::Vec_T, 1U> const GaussQR<Quad, 1U>::node = {
    {GaussQR<Quad, 1U>::Vec_T(0.L, 0.L)}};
template <>
short_T const GaussQR<Quad, 1U>::bestPt = 0;

template <>
GaussQR<Quad, 4U>::Weights_T const GaussQR<Quad, 4U>::weight =
    GaussQR<Quad, 4U>::Weights_T::Constant(1.L);
// clang-format off
template <>
std::array<GaussQR<Quad, 4U>::Vec_T, 4U> const GaussQR<Quad, 4U>::node = {{
    GaussQR<Quad, 4U>::Vec_T(-sqrt13rd, -sqrt13rd),
    GaussQR<Quad, 4U>::Vec_T( sqrt13rd, -sqrt13rd),
    GaussQR<Quad, 4U>::Vec_T(-sqrt13rd,  sqrt13rd),
    GaussQR<Quad, 4U>::Vec_T( sqrt13rd,  sqrt13rd),
}};
// clang-format on
template <>
short_T const GaussQR<Quad, 4U>::bestPt = 0;

// clang-format off
template <>
GaussQR<Quad, 9U>::Weights_T const
    GaussQR<Quad, 9U>::weight = (GaussQR<Quad, 9U>::Weights_T() <<
        25.L, 40.L, 25.L,
        40.L, 64.L, 40.L,
        25.L, 40.L, 25.L).finished() / 81.L;
template <>
std::array<GaussQR<Quad, 9U>::Vec_T, 9U> const GaussQR<Quad, 9U>::node = {
    {GaussQR<Quad, 9U>::Vec_T(-sqrt35th, -sqrt35th),
     GaussQR<Quad, 9U>::Vec_T(       0., -sqrt35th),
     GaussQR<Quad, 9U>::Vec_T( sqrt35th, -sqrt35th),
     GaussQR<Quad, 9U>::Vec_T(-sqrt35th,        0.),
     GaussQR<Quad, 9U>::Vec_T(       0.,        0.),
     GaussQR<Quad, 9U>::Vec_T( sqrt35th,        0.),
     GaussQR<Quad, 9U>::Vec_T(-sqrt35th, sqrt35th),
     GaussQR<Quad, 9U>::Vec_T(       0., sqrt35th),
     GaussQR<Quad, 9U>::Vec_T( sqrt35th, sqrt35th)}};
// clang-format on
template <>
short_T const GaussQR<Quad, 9U>::bestPt = 4;

// 3D
// - Tetrahedron -----------------------------------------------------------------------
template <>
GaussQR<Tetrahedron, 1U>::Weights_T const GaussQR<Tetrahedron, 1U>::weight =
    GaussQR<Tetrahedron, 1U>::Weights_T::Constant(1.L / 6);
template <>
std::array<GaussQR<Tetrahedron, 1U>::Vec_T, 1U> const GaussQR<Tetrahedron, 1U>::node = {
    {GaussQR<Tetrahedron, 1U>::Vec_T(.25L, .25L, .25L)}};
template <>
short_T const GaussQR<Tetrahedron, 1U>::bestPt = 0;

template <>
GaussQR<Tetrahedron, 4U>::Weights_T const GaussQR<Tetrahedron, 4U>::weight =
    GaussQR<Tetrahedron, 4U>::Weights_T::Constant(.25L / 6);
// clang-format off
template <>
std::array<GaussQR<Tetrahedron, 4U>::Vec_T, 4U> const GaussQR<Tetrahedron, 4U>::node = {
    {GaussQR<Tetrahedron, 4U>::Vec_T(0.1381966011250105L, 0.1381966011250105L, 0.1381966011250105L),
     GaussQR<Tetrahedron, 4U>::Vec_T(0.5854101966249685L, 0.1381966011250105L, 0.1381966011250105L),
     GaussQR<Tetrahedron, 4U>::Vec_T(0.1381966011250105L, 0.5854101966249685L, 0.1381966011250105L),
     GaussQR<Tetrahedron, 4U>::Vec_T(0.1381966011250105L, 0.1381966011250105L, 0.5854101966249685L)}};
// clang-format on
template <>
short_T const GaussQR<Tetrahedron, 4U>::bestPt = 0;

// clang-format off
template <>
GaussQR<Tetrahedron, 11U>::Weights_T const GaussQR<Tetrahedron, 11U>::weight =
    (GaussQR<Tetrahedron, 11U>::Weights_T() <<
        -0.0131555555555556L,
         0.0076222222222222L, 0.0076222222222222L, 0.0076222222222222L, 0.0076222222222222L,
         0.024888888888889L, 0.024888888888889L, 0.024888888888889L, 0.024888888888889L, 0.024888888888889L, 0.024888888888889L).finished();
template <>
std::array<GaussQR<Tetrahedron, 11U>::Vec_T, 11U> const GaussQR<Tetrahedron, 11U>::node = {{
    GaussQR<Tetrahedron, 11U>::Vec_T(              0.25L,               0.25L,               0.25L),
    GaussQR<Tetrahedron, 11U>::Vec_T(0.0714285714285714L, 0.0714285714285714L, 0.0714285714285714L),
    GaussQR<Tetrahedron, 11U>::Vec_T(0.7857142857142857L, 0.0714285714285714L, 0.0714285714285714L),
    GaussQR<Tetrahedron, 11U>::Vec_T(0.0714285714285714L, 0.7857142857142857L, 0.0714285714285714L),
    GaussQR<Tetrahedron, 11U>::Vec_T(0.0714285714285714L, 0.0714285714285714L, 0.7857142857142857L),
    GaussQR<Tetrahedron, 11U>::Vec_T(0.399403576166799L,  0.399403576166799L,  0.100596423833201L ),
    GaussQR<Tetrahedron, 11U>::Vec_T(0.399403576166799L,  0.100596423833201L,  0.399403576166799L ),
    GaussQR<Tetrahedron, 11U>::Vec_T(0.100596423833201L,  0.399403576166799L,  0.399403576166799L ),
    GaussQR<Tetrahedron, 11U>::Vec_T(0.399403576166799L,  0.100596423833201L,  0.100596423833201L ),
    GaussQR<Tetrahedron, 11U>::Vec_T(0.100596423833201L,  0.399403576166799L,  0.100596423833201L ),
    GaussQR<Tetrahedron, 11U>::Vec_T(0.100596423833201L,  0.100596423833201L,  0.399403576166799L),
}};
// clang-format on
template <>
short_T const GaussQR<Tetrahedron, 11U>::bestPt = 0;

// 3D
// - Hexahedron ------------------------------------------------------------------------
template <>
GaussQR<Hexahedron, 1U>::Weights_T const GaussQR<Hexahedron, 1U>::weight =
    GaussQR<Hexahedron, 1U>::Weights_T::Constant(8.L);
template <>
std::array<GaussQR<Hexahedron, 1U>::Vec_T, 1U> const GaussQR<Hexahedron, 1U>::node = {{
    GaussQR<Hexahedron, 1U>::Vec_T(0.L, 0.L, 0.L),
}};
template <>
short_T const GaussQR<Hexahedron, 1U>::bestPt = 0;

template <>
GaussQR<Hexahedron, 8U>::Weights_T const GaussQR<Hexahedron, 8U>::weight =
    GaussQR<Hexahedron, 8U>::Weights_T::Constant(1.L);
// clang-format off
template <>
std::array<GaussQR<Hexahedron, 8U>::Vec_T, 8U> const GaussQR<Hexahedron, 8U>::node = {{
    GaussQR<Hexahedron, 8U>::Vec_T(-sqrt13rd, -sqrt13rd, -sqrt13rd),
    GaussQR<Hexahedron, 8U>::Vec_T( sqrt13rd, -sqrt13rd, -sqrt13rd),
    GaussQR<Hexahedron, 8U>::Vec_T(-sqrt13rd,  sqrt13rd, -sqrt13rd),
    GaussQR<Hexahedron, 8U>::Vec_T( sqrt13rd,  sqrt13rd, -sqrt13rd),
    GaussQR<Hexahedron, 8U>::Vec_T(-sqrt13rd, -sqrt13rd,  sqrt13rd),
    GaussQR<Hexahedron, 8U>::Vec_T( sqrt13rd, -sqrt13rd,  sqrt13rd),
    GaussQR<Hexahedron, 8U>::Vec_T(-sqrt13rd,  sqrt13rd,  sqrt13rd),
    GaussQR<Hexahedron, 8U>::Vec_T( sqrt13rd,  sqrt13rd,  sqrt13rd),
}};
// clang-format on
template <>
short_T const GaussQR<Hexahedron, 8U>::bestPt = 0;

// clang-format off
template <>
GaussQR<Hexahedron, 27U>::Weights_T const
    GaussQR<Hexahedron, 27U>::weight = (GaussQR<Hexahedron, 27U>::Weights_T() <<
        125.L, 200.L, 125.L, 200.L, 320.L, 200.L, 125.L, 200.L, 125.L,
        200.L, 320.L, 200.L, 320.L, 512.L, 320.L, 200.L, 320.L, 200.L,
        125.L, 200.L, 125.L, 200.L, 320.L, 200.L, 125.L, 200.L, 125.L).finished() / 729.L;
template <>
std::array<GaussQR<Hexahedron, 27U>::Vec_T, 27U> const GaussQR<Hexahedron, 27U>::node = {{
    GaussQR<Hexahedron, 27U>::Vec_T(-sqrt35th, -sqrt35th, -sqrt35th),
    GaussQR<Hexahedron, 27U>::Vec_T(      0.L, -sqrt35th, -sqrt35th),
    GaussQR<Hexahedron, 27U>::Vec_T( sqrt35th, -sqrt35th, -sqrt35th),
    GaussQR<Hexahedron, 27U>::Vec_T(-sqrt35th,       0.L, -sqrt35th),
    GaussQR<Hexahedron, 27U>::Vec_T(      0.L,       0.L, -sqrt35th),
    GaussQR<Hexahedron, 27U>::Vec_T( sqrt35th,       0.L, -sqrt35th),
    GaussQR<Hexahedron, 27U>::Vec_T(-sqrt35th,  sqrt35th, -sqrt35th),
    GaussQR<Hexahedron, 27U>::Vec_T(      0.L,  sqrt35th, -sqrt35th),
    GaussQR<Hexahedron, 27U>::Vec_T( sqrt35th,  sqrt35th, -sqrt35th),
    //
    GaussQR<Hexahedron, 27U>::Vec_T(-sqrt35th, -sqrt35th, 0.L),
    GaussQR<Hexahedron, 27U>::Vec_T(      0.L, -sqrt35th, 0.L),
    GaussQR<Hexahedron, 27U>::Vec_T( sqrt35th, -sqrt35th, 0.L),
    GaussQR<Hexahedron, 27U>::Vec_T(-sqrt35th,       0.L, 0.L),
    GaussQR<Hexahedron, 27U>::Vec_T(      0.L,       0.L, 0.L),
    GaussQR<Hexahedron, 27U>::Vec_T( sqrt35th,       0.L, 0.L),
    GaussQR<Hexahedron, 27U>::Vec_T(-sqrt35th,  sqrt35th, 0.L),
    GaussQR<Hexahedron, 27U>::Vec_T(      0.L,  sqrt35th, 0.L),
    GaussQR<Hexahedron, 27U>::Vec_T( sqrt35th,  sqrt35th, 0.L),
    //
    GaussQR<Hexahedron, 27U>::Vec_T(-sqrt35th, -sqrt35th, sqrt35th),
    GaussQR<Hexahedron, 27U>::Vec_T(      0.L, -sqrt35th, sqrt35th),
    GaussQR<Hexahedron, 27U>::Vec_T( sqrt35th, -sqrt35th, sqrt35th),
    GaussQR<Hexahedron, 27U>::Vec_T(-sqrt35th,       0.L, sqrt35th),
    GaussQR<Hexahedron, 27U>::Vec_T(      0.L,       0.L, sqrt35th),
    GaussQR<Hexahedron, 27U>::Vec_T( sqrt35th,       0.L, sqrt35th),
    GaussQR<Hexahedron, 27U>::Vec_T(-sqrt35th,  sqrt35th, sqrt35th),
    GaussQR<Hexahedron, 27U>::Vec_T(      0.L,  sqrt35th, sqrt35th),
    GaussQR<Hexahedron, 27U>::Vec_T( sqrt35th,  sqrt35th, sqrt35th),
}};
// clang-format on
template <>
short_T const GaussQR<Hexahedron, 27U>::bestPt = 13;

// SideGaussQR =========================================================================

// 1D
template <>
SideGaussQR<Line, 1u>::Weights_T const SideGaussQR<Line, 1u>::weight =
    SideGaussQR<Line, 1u>::Weights_T::Constant(1.L);
template <>
std::array<SideGaussQR<Line, 1u>::Vec_T, 1u * 2u> const SideGaussQR<Line, 1u>::node = {{
    FVec<1u>{1.0},
    FVec<1u>{1.0},
}};

// 2D
// - Triangle --------------------------------------------------------------------------
template <>
SideGaussQR<Triangle, 2U>::Weights_T const SideGaussQR<Triangle, 2U>::weight =
    (SideGaussQR<Triangle, 2U>::Weights_T{} << 1.L, 1.L, 1.L, 1.L, 1.L, 1.L).finished();
template <>
std::array<SideGaussQR<Triangle, 2U>::Vec_T, 2U * 3U> const
    SideGaussQR<Triangle, 2U>::node = {{
        SideGaussQR<Triangle, 2U>::Vec_T{0.5 * (1. - sqrt13rd), 0.0},
        SideGaussQR<Triangle, 2U>::Vec_T{0.5 * (1. + sqrt13rd), 0.0},

        SideGaussQR<Triangle, 2U>::Vec_T{0.5 * (1. + sqrt13rd), 0.5 * (1. - sqrt13rd)},
        SideGaussQR<Triangle, 2U>::Vec_T{0.5 * (1. - sqrt13rd), 0.5 * (1. + sqrt13rd)},

        SideGaussQR<Triangle, 2U>::Vec_T{0.0, 0.5 * (1. - sqrt13rd)},
        SideGaussQR<Triangle, 2U>::Vec_T{0.0, 0.5 * (1. + sqrt13rd)},
    }};

// clang-format off
template <>
SideGaussQR<Triangle, 3U>::Weights_T const SideGaussQR<Triangle, 3U>::weight =
    (SideGaussQR<Triangle, 3U>::Weights_T{} << 5.L, 8.L, 5.L, 5.L, 8.L, 5.L, 5.L, 8.L, 5.L).finished() / 9.L;
template <>
std::array<SideGaussQR<Triangle, 3U>::Vec_T, 3U * 3U> const
    SideGaussQR<Triangle, 3U>::node = {{
        SideGaussQR<Triangle, 3U>::Vec_T{0.5 * (1. - sqrt35th), 0.0},
        SideGaussQR<Triangle, 3U>::Vec_T{0.5, 0.0},
        SideGaussQR<Triangle, 3U>::Vec_T{0.5 * (1. + sqrt35th), 0.0},

        SideGaussQR<Triangle, 3U>::Vec_T{0.5 * (1. + sqrt35th), 0.5 * (1. - sqrt35th)},
        SideGaussQR<Triangle, 3U>::Vec_T{0.5, 0.5},
        SideGaussQR<Triangle, 3U>::Vec_T{0.5 * (1. - sqrt35th), 0.5 * (1. + sqrt35th)},

        SideGaussQR<Triangle, 3U>::Vec_T{0.0, 0.5 * (1. - sqrt35th)},
        SideGaussQR<Triangle, 3U>::Vec_T{0.0, 0.5},
        SideGaussQR<Triangle, 3U>::Vec_T{0.0, 0.5 * (1. + sqrt35th)},
    }};
// clang-format on

// - Quad ------------------------------------------------------------------------------
template <>
SideGaussQR<Quad, 2U>::Weights_T const SideGaussQR<Quad, 2U>::weight =
    (SideGaussQR<Quad, 2U>::Weights_T{} << 1.L, 1.L, 1.L, 1.L, 1.L, 1.L, 1.L, 1.L)
        .finished();
template <>
std::array<SideGaussQR<Quad, 2U>::Vec_T, 2U * 4U> const SideGaussQR<Quad, 2U>::node = {{
    // bottom
    SideGaussQR<Quad, 2U>::Vec_T{-sqrt13rd, -1.0},
    SideGaussQR<Quad, 2U>::Vec_T{sqrt13rd, -1.0},
    // right
    SideGaussQR<Quad, 2U>::Vec_T{1.0, -sqrt13rd},
    SideGaussQR<Quad, 2U>::Vec_T{1.0, sqrt13rd},
    // top
    SideGaussQR<Quad, 2U>::Vec_T{-sqrt13rd, 1.0},
    SideGaussQR<Quad, 2U>::Vec_T{sqrt13rd, 1.0},
    // left
    SideGaussQR<Quad, 2U>::Vec_T{-1.0, -sqrt13rd},
    SideGaussQR<Quad, 2U>::Vec_T{-1.0, sqrt13rd},
}};

// clang-format off
template <>
SideGaussQR<Quad, 3U>::Weights_T const
    SideGaussQR<Quad, 3U>::weight = (SideGaussQR<Quad, 3U>::Weights_T{} <<
        5.L, 8.L, 5.L,
        5.L, 8.L, 5.L,
        5.L, 8.L, 5.L,
        5.L, 8.L, 5.L).finished() / 9.L;
template <>
std::array<SideGaussQR<Quad, 3U>::Vec_T, 3U * 4U> const SideGaussQR<Quad, 3U>::node = {{
    // bottom
    SideGaussQR<Quad, 3U>::Vec_T{-sqrt35th, -1.0},
    SideGaussQR<Quad, 3U>::Vec_T{     0.0L, -1.0},
    SideGaussQR<Quad, 3U>::Vec_T{ sqrt35th, -1.0},
    // right
    SideGaussQR<Quad, 3U>::Vec_T{1.0, -sqrt35th},
    SideGaussQR<Quad, 3U>::Vec_T{1.0,      0.0L},
    SideGaussQR<Quad, 3U>::Vec_T{1.0,  sqrt35th},
    // top
    SideGaussQR<Quad, 3U>::Vec_T{-sqrt35th, 1.0},
    SideGaussQR<Quad, 3U>::Vec_T{     0.0L, 1.0},
    SideGaussQR<Quad, 3U>::Vec_T{ sqrt35th, 1.0},
    // left
    SideGaussQR<Quad, 3U>::Vec_T{-1.0, -sqrt35th},
    SideGaussQR<Quad, 3U>::Vec_T{-1.0,      0.0L},
    SideGaussQR<Quad, 3U>::Vec_T{-1.0,  sqrt35th},
}};
// clang-format on

// 3D
// - Tetrahedron -----------------------------------------------------------------------
// clang-format off
template <>
SideGaussQR<Tetrahedron, 3U>::Weights_T const SideGaussQR<Tetrahedron, 3U>::weight =
    (SideGaussQR<Tetrahedron, 3U>::Weights_T{} <<
     1.L, 1.L, 1.L,
     1.L, 1.L, 1.L,
     1.L, 1.L, 1.L,
     sqrt(3) / 4., sqrt(3) / 4., sqrt(3) / 4.)
        .finished() /
    6.L;
template <>
std::array<SideGaussQR<Tetrahedron, 3U>::Vec_T, 3U * 4U> const
    SideGaussQR<Tetrahedron, 3U>::node = {{
        SideGaussQR<Tetrahedron, 3U>::Vec_T{1. / 6, 1. / 6, 0.0},
        SideGaussQR<Tetrahedron, 3U>::Vec_T{2. / 3, 1. / 6, 0.0},
        SideGaussQR<Tetrahedron, 3U>::Vec_T{1. / 6, 2. / 3, 0.0},

        SideGaussQR<Tetrahedron, 3U>::Vec_T{1. / 6, 0.0, 1. / 6},
        SideGaussQR<Tetrahedron, 3U>::Vec_T{2. / 3, 0.0, 1. / 6},
        SideGaussQR<Tetrahedron, 3U>::Vec_T{1. / 6, 0.0, 2. / 3},

        SideGaussQR<Tetrahedron, 3U>::Vec_T{0.0, 1. / 6, 1. / 6},
        SideGaussQR<Tetrahedron, 3U>::Vec_T{0.0, 2. / 3, 1. / 6},
        SideGaussQR<Tetrahedron, 3U>::Vec_T{0.0, 1. / 6, 2. / 3},

        SideGaussQR<Tetrahedron, 3U>::Vec_T{1. / 6, 1. / 6, 2. / 3},
        SideGaussQR<Tetrahedron, 3U>::Vec_T{2. / 3, 1. / 6, 1. / 6},
        SideGaussQR<Tetrahedron, 3U>::Vec_T{1. / 6, 2. / 3, 1. / 6},
}};
// clang-format on

// 3D
// - Hexahedron ------------------------------------------------------------------------
// clang-format off
template <>
SideGaussQR<Hexahedron, 4U>::Weights_T const SideGaussQR<Hexahedron, 4U>::weight =
    (SideGaussQR<Hexahedron, 4U>::Weights_T{} <<
     0.5L, 0.5L, 0.5L, 0.5L,
     0.5L, 0.5L, 0.5L, 0.5L,
     0.5L, 0.5L, 0.5L, 0.5L,
     0.5L, 0.5L, 0.5L, 0.5L,
     0.5L, 0.5L, 0.5L, 0.5L,
     0.5L, 0.5L, 0.5L, 0.5L)
        .finished();
template <>
std::array<SideGaussQR<Hexahedron, 4U>::Vec_T, 4U * 6U> const
    SideGaussQR<Hexahedron, 4U>::node = {{
        // back
        SideGaussQR<Hexahedron, 4U>::Vec_T{-sqrt13rd, -sqrt13rd, -1.0},
        SideGaussQR<Hexahedron, 4U>::Vec_T{+sqrt13rd, -sqrt13rd, -1.0},
        SideGaussQR<Hexahedron, 4U>::Vec_T{+sqrt13rd, +sqrt13rd, -1.0},
        SideGaussQR<Hexahedron, 4U>::Vec_T{-sqrt13rd, +sqrt13rd, -1.0},
        // left
        SideGaussQR<Hexahedron, 4U>::Vec_T{-1.0, -sqrt13rd, -sqrt13rd},
        SideGaussQR<Hexahedron, 4U>::Vec_T{-1.0, +sqrt13rd, -sqrt13rd},
        SideGaussQR<Hexahedron, 4U>::Vec_T{-1.0, +sqrt13rd, +sqrt13rd},
        SideGaussQR<Hexahedron, 4U>::Vec_T{-1.0, -sqrt13rd, +sqrt13rd},
        // bottom
        SideGaussQR<Hexahedron, 4U>::Vec_T{-sqrt13rd, -1.0, -sqrt13rd},
        SideGaussQR<Hexahedron, 4U>::Vec_T{+sqrt13rd, -1.0, -sqrt13rd},
        SideGaussQR<Hexahedron, 4U>::Vec_T{+sqrt13rd, -1.0, +sqrt13rd},
        SideGaussQR<Hexahedron, 4U>::Vec_T{-sqrt13rd, -1.0, +sqrt13rd},
        // right
        SideGaussQR<Hexahedron, 4U>::Vec_T{+1.0, -sqrt13rd, -sqrt13rd},
        SideGaussQR<Hexahedron, 4U>::Vec_T{+1.0, +sqrt13rd, -sqrt13rd},
        SideGaussQR<Hexahedron, 4U>::Vec_T{+1.0, +sqrt13rd, +sqrt13rd},
        SideGaussQR<Hexahedron, 4U>::Vec_T{+1.0, -sqrt13rd, +sqrt13rd},
        // top
        SideGaussQR<Hexahedron, 4U>::Vec_T{-sqrt13rd, +1.0, -sqrt13rd},
        SideGaussQR<Hexahedron, 4U>::Vec_T{+sqrt13rd, +1.0, -sqrt13rd},
        SideGaussQR<Hexahedron, 4U>::Vec_T{+sqrt13rd, +1.0, +sqrt13rd},
        SideGaussQR<Hexahedron, 4U>::Vec_T{-sqrt13rd, +1.0, +sqrt13rd},
        // front
        SideGaussQR<Hexahedron, 4U>::Vec_T{-sqrt13rd, -sqrt13rd, +1.0},
        SideGaussQR<Hexahedron, 4U>::Vec_T{+sqrt13rd, -sqrt13rd, +1.0},
        SideGaussQR<Hexahedron, 4U>::Vec_T{+sqrt13rd, +sqrt13rd, +1.0},
        SideGaussQR<Hexahedron, 4U>::Vec_T{-sqrt13rd, +sqrt13rd, +1.0},
}};
// clang-format on

// TrapQR ==============================================================================

// 1D ----------------------------------------------------------------------------------
template <>
FVec<2U> const TrapQR<Line>::weight = FVec<2U>::Constant(1.L);
template <>
std::array<TrapQR<Line>::Vec_T, 2U> const TrapQR<Line>::node = {
    {TrapQR<Line>::Vec_T::Constant(-1.L), TrapQR<Line>::Vec_T::Constant(1.L)}};

// 2D ----------------------------------------------------------------------------------
template <>
FVec<3U> const TrapQR<Triangle>::weight = FVec<3U>::Constant(1.L / 6);
template <>
std::array<TrapQR<Triangle>::Vec_T, 3U> const TrapQR<Triangle>::node = {
    {TrapQR<Triangle>::Vec_T(0., 0.),
     TrapQR<Triangle>::Vec_T(1., 0.),
     TrapQR<Triangle>::Vec_T(0., 1.)}};

template <>
FVec<4U> const TrapQR<Quad>::weight = FVec<4U>::Constant(1.L);
template <>
std::array<TrapQR<Quad>::Vec_T, 4U> const TrapQR<Quad>::node = {
    {TrapQR<Quad>::Vec_T(-1., -1.),
     TrapQR<Quad>::Vec_T(-1., 1.),
     TrapQR<Quad>::Vec_T(1., -1.),
     TrapQR<Quad>::Vec_T(1., 1.)}};

// 3D ----------------------------------------------------------------------------------
template <>
FVec<4U> const TrapQR<Tetrahedron>::weight = FVec<4U>::Constant(1.L);
template <>
std::array<TrapQR<Tetrahedron>::Vec_T, 4U> const TrapQR<Tetrahedron>::node = {
    {TrapQR<Tetrahedron>::Vec_T(0., 0., 0.),
     TrapQR<Tetrahedron>::Vec_T(1., 0., 0.),
     TrapQR<Tetrahedron>::Vec_T(0., 1., 0.),
     TrapQR<Tetrahedron>::Vec_T(0., 0., 1.)}};

template <>
FVec<8U> const TrapQR<Hexahedron>::weight = FVec<8U>::Constant(1.L);
// clang-format off
template <>
std::array<TrapQR<Hexahedron>::Vec_T, 8U> const TrapQR<Hexahedron>::node = {
    {TrapQR<Hexahedron>::Vec_T(-1., -1., -1.),
     TrapQR<Hexahedron>::Vec_T(-1.,  1., -1.),
     TrapQR<Hexahedron>::Vec_T( 1., -1., -1.),
     TrapQR<Hexahedron>::Vec_T( 1.,  1., -1.),
     TrapQR<Hexahedron>::Vec_T(-1., -1.,  1.),
     TrapQR<Hexahedron>::Vec_T(-1.,  1.,  1.),
     TrapQR<Hexahedron>::Vec_T( 1., -1.,  1.),
     TrapQR<Hexahedron>::Vec_T( 1.,  1.,  1.)}};
// clang-format on

// SimpsonQR ===========================================================================

// 1D ----------------------------------------------------------------------------------
template <>
FVec<3U> const SimpsonQR<Line>::weight = (FVec<3U>() << 1., 4., 1.).finished() / 3.;
// clang-format off
template <>
std::array<SimpsonQR<Line>::Vec_T, 3U> const SimpsonQR<Line>::node = {{
    SimpsonQR<Line>::Vec_T{-1.},
    SimpsonQR<Line>::Vec_T{ 0.},
    SimpsonQR<Line>::Vec_T{ 1.},
}};
// clang-format on

// 2D ----------------------------------------------------------------------------------
template <>
FVec<9> const SimpsonQR<Quad>::weight =
    (FVec<9>() << 1., 4., 1., 4., 16., 4., 1., 4., 1.).finished() / 9.;
// clang-format off
template <>
std::array<SimpsonQR<Quad>::Vec_T, 9U> const SimpsonQR<Quad>::node = {{
    SimpsonQR<Quad>::Vec_T{-1., -1.},
    SimpsonQR<Quad>::Vec_T{ 0., -1.},
    SimpsonQR<Quad>::Vec_T{ 1., -1.},
    SimpsonQR<Quad>::Vec_T{-1.,  0.},
    SimpsonQR<Quad>::Vec_T{ 0.,  0.},
    SimpsonQR<Quad>::Vec_T{ 1.,  0.},
    SimpsonQR<Quad>::Vec_T{-1.,  1.},
    SimpsonQR<Quad>::Vec_T{ 0.,  1.},
    SimpsonQR<Quad>::Vec_T{ 1.,  1.},
}};
// clang-format on

// 3D ----------------------------------------------------------------------------------
// clang-format off
template <>
FVec<27> const SimpsonQR<Hexahedron>::weight = (FVec<27>() <<
    1.,  4., 1.,  4., 16.,  4., 1.,  4., 1.,
    4., 16., 4., 16., 64., 16., 4., 16., 4.,
    1.,  4., 1.,  4., 16.,  4., 1.,  4., 1.).finished() / 27.;
template <>
std::array<SimpsonQR<Hexahedron>::Vec_T, 27U> const SimpsonQR<Hexahedron>::node = {{
    SimpsonQR<Hexahedron>::Vec_T{-1., -1., -1.},
    SimpsonQR<Hexahedron>::Vec_T{ 0., -1., -1.},
    SimpsonQR<Hexahedron>::Vec_T{ 1., -1., -1.},
    SimpsonQR<Hexahedron>::Vec_T{-1.,  0., -1.},
    SimpsonQR<Hexahedron>::Vec_T{ 0.,  0., -1.},
    SimpsonQR<Hexahedron>::Vec_T{ 1.,  0., -1.},
    SimpsonQR<Hexahedron>::Vec_T{-1.,  1., -1.},
    SimpsonQR<Hexahedron>::Vec_T{ 0.,  1., -1.},
    SimpsonQR<Hexahedron>::Vec_T{ 1.,  1., -1.},
    SimpsonQR<Hexahedron>::Vec_T{-1., -1.,  0.},
    SimpsonQR<Hexahedron>::Vec_T{ 0., -1.,  0.},
    SimpsonQR<Hexahedron>::Vec_T{ 1., -1.,  0.},
    SimpsonQR<Hexahedron>::Vec_T{-1.,  0.,  0.},
    SimpsonQR<Hexahedron>::Vec_T{ 0.,  0.,  0.},
    SimpsonQR<Hexahedron>::Vec_T{ 1.,  0.,  0.},
    SimpsonQR<Hexahedron>::Vec_T{-1.,  1.,  0.},
    SimpsonQR<Hexahedron>::Vec_T{ 0.,  1.,  0.},
    SimpsonQR<Hexahedron>::Vec_T{ 1.,  1.,  0.},
    SimpsonQR<Hexahedron>::Vec_T{-1., -1.,  1.},
    SimpsonQR<Hexahedron>::Vec_T{ 0., -1.,  1.},
    SimpsonQR<Hexahedron>::Vec_T{ 1., -1.,  1.},
    SimpsonQR<Hexahedron>::Vec_T{-1.,  0.,  1.},
    SimpsonQR<Hexahedron>::Vec_T{ 0.,  0.,  1.},
    SimpsonQR<Hexahedron>::Vec_T{ 1.,  0.,  1.},
    SimpsonQR<Hexahedron>::Vec_T{-1.,  1.,  1.},
    SimpsonQR<Hexahedron>::Vec_T{ 0.,  1.,  1.},
    SimpsonQR<Hexahedron>::Vec_T{ 1.,  1.,  1.},
}};
// clang-format on

} // namespace proxpde
