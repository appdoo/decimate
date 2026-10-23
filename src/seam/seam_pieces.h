// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

/// @file seam_pieces.h
/// @brief Internal: exact piecewise-polynomial form of the seam mismatch (not installed).
///
/// A seam pair is split at every edge parameter gamma where either half-edge crosses a
/// texel-grid line.  On each piece both bilinear reconstructions (and their inward
/// normal derivatives) are polynomials of degree <= 2 in gamma, so one `SeamPiece`
/// carries the 8-texel coefficient vectors of the difference (or sum) polynomial.
/// `seamEnergy` and the erase solver both consume these pieces -- one unit, one set of
/// sampling rules.

#include <array>
#include <decimate/math/vec2.h>
#include <decimate/math/vec3.h>
#include <decimate/math/vec4.h>
#include <decimate/seam/seam_topology.h>
#include <decimate/types/textured_mesh.h>
#include <span>
#include <vector>

namespace decimate::seam::detail {

/// @brief Texture extent in texels.
struct TexelGrid {
    int width = 0;
    int height = 0;
};


/// @brief A UV segment in (unclamped) texel coordinates: `uv * N - 0.5`.
struct Segment {
    double x0 = 0.0;
    double y0 = 0.0;
    double x1 = 0.0;
    double y1 = 0.0;
};


/// @brief Lower-left texel of a 2x2 bilinear cell.
struct Cell {
    int px = 0;
    int py = 0;
};


/// @brief Which polynomial a piece describes.
enum class PieceKind {
    Value,    ///< Bilinear(A) - Bilinear(B): the seam mismatch (paper eq. 3-8).
    Gradient, ///< d/dn_in Bilinear(A) + d/dn_in Bilinear(B): the C1 term.
};


/// @brief The tangent frames of a mesh: what a tangent-space normal map is decoded with.
///
/// `tangents[3*f+k]` is T (xyz) and the handedness w of corner k of face f; `normals[v]` is the
/// vertex normal of surface vertex v.  The engine decodes `F s` with `F = [T | w N x T | N]`.
struct FrameField {
    std::span<const math::Vec4f> tangents;
    std::span<const math::Vec3f> normals;
};


/// @brief T and N of one half-edge along gamma: linear and unnormalised, as the engine interpolates them.
/// The default is the identity frame (T = x, N = z, w = 1: F = I), so an unframed piece is the scalar case.
struct EdgeFrame {
    math::Vec3d t0 { 1.0, 0.0, 0.0 }, t1 { 1.0, 0.0, 0.0 }; ///< T at gamma 0 and 1.
    math::Vec3d n0 { 0.0, 0.0, 1.0 }, n1 { 0.0, 0.0, 1.0 }; ///< N at gamma 0 and 1.
    double      w = 1.0;
};


/// @brief One gamma-interval of a seam pair: phi(g) = c2 g^2 + c1 g + c0 over 8 texels.
///
/// `texel[0..3]` are the cell of half-edge A (p00, p10, p01, p11), `texel[4..7]` that of
/// half-edge B.  Texel indices are row-major (`y * width + x`).
struct SeamPiece {
    std::array<int, 8>       texel {};
    std::array<double, 8>    c2 {};
    std::array<double, 8>    c1 {};
    std::array<double, 8>    c0 {};
    double                   a = 0.0;        ///< Interval start in [0,1].
    double                   b = 0.0;        ///< Interval end in [0,1].
    double                   length3D = 0.0; ///< 3D length of the edge this piece belongs to.
    std::array<EdgeFrame, 2> frame {};       ///< Frames of side A and B; the identity unless `seamPieces` got a `FrameField`.
};


/// @brief UV point -> texel coordinate (`uv * N - 0.5`, no clamping).
[[nodiscard]] Segment toTexelSegment(math::Vec2f uvStart, math::Vec2f uvEnd, TexelGrid grid);


/// @brief The UV segment of half-edge @p he, start -> end, in texel coordinates.
[[nodiscard]] Segment halfEdgeSegment(const TexturedMesh &tm, HalfEdgeRef he, TexelGrid grid);


/// @brief Bilinear cells visited by @p seg, one per sub-interval (duplicates possible).
[[nodiscard]] std::vector<Cell> cellsAlong(const Segment &seg, TexelGrid grid);


/// @brief Pieces of every seam pair of @p topo.
///
/// Half-edge B is traversed end-to-start so that gamma addresses the same 3D point on
/// both sides.  Intervals shorter than 1e-12 are dropped.  With @p frames every piece also
/// carries the frames of its two half-edges (the tangent-space mismatch, see `frameResidual`).
[[nodiscard]] std::vector<SeamPiece> seamPieces(const TexturedMesh &tm, const SeamTopology &topo,
                                                TexelGrid grid, PieceKind kind,
                                                const FrameField *frames = nullptr);


/// @brief Sum of the 3D lengths of all seam pairs.
[[nodiscard]] double totalSeamLength(const SeamTopology &topo);


/// @brief True when @p topo has seams of positive total 3D length.
///
/// Seams whose edges all have length 0 (coincident vertices) carry no energy and would divide the
/// normalization by 0, so callers treat them as "no seam".
[[nodiscard]] bool hasSeamLength(const SeamTopology &topo);


/// @brief Gauss-Legendre 5-point rule on [a,b] (exact to degree 9): the framed mismatch has degree 8.
struct GaussRule5 {
    std::array<double, 5> gamma {};
    std::array<double, 5> weight {};
};


/// @brief Nodes and weights of the 5-point rule on [@p a, @p b].
[[nodiscard]] GaussRule5 gauss5On(double a, double b);


/// @brief The 8 row-major 3x3 blocks Phi_i = phi_i(gamma) F_side of a Value piece.
///
/// The residual of the framed mismatch is `sum_i Phi_i x_i`, x_i the texel's 3 values (B carries the
/// minus sign inside phi).
[[nodiscard]] std::array<std::array<double, 9>, 8> frameBlocks(const SeamPiece &piece, double gamma);


/// @brief Per world component: seam-length weighted integral of the (framed) mismatch squared, over @p total
///        (paper eq. 9).  With the default identity frames this is the per-channel D_total.
///
/// @param pieces       Output of `seamPieces(..., PieceKind::Value)`.
/// @param totalLength  `totalSeamLength(topo)`; must be > 0.
/// @param texels       `width * height * 3` interleaved values.
[[nodiscard]] std::array<double, 3> integrateEnergy(std::span<const SeamPiece> pieces, double totalLength,
                                                    std::span<const double> texels);


/// @brief phi(gamma) of @p piece, all eight coefficients.
[[nodiscard]] std::array<double, 8> evaluatePiece(const SeamPiece &piece, double gamma);

} // namespace decimate::seam::detail
