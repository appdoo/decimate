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


/// @brief One gamma-interval of a seam pair: phi(g) = c2 g^2 + c1 g + c0 over 8 texels.
///
/// `texel[0..3]` are the cell of half-edge A (p00, p10, p01, p11), `texel[4..7]` that of
/// half-edge B.  Texel indices are row-major (`y * width + x`).
struct SeamPiece {
    std::array<int, 8>    texel {};
    std::array<double, 8> c2 {};
    std::array<double, 8> c1 {};
    std::array<double, 8> c0 {};
    double                a = 0.0;        ///< Interval start in [0,1].
    double                b = 0.0;        ///< Interval end in [0,1].
    double                length3D = 0.0; ///< 3D length of the edge this piece belongs to.
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
/// both sides.  Intervals shorter than 1e-12 are dropped.
[[nodiscard]] std::vector<SeamPiece> seamPieces(const TexturedMesh &tm, const SeamTopology &topo,
                                                TexelGrid grid, PieceKind kind);


/// @brief Sum of the 3D lengths of all seam pairs.
[[nodiscard]] double totalSeamLength(const SeamTopology &topo);


/// @brief True when @p topo has seams of positive total 3D length.
///
/// Seams whose edges all have length 0 (coincident vertices) carry no energy and would divide the
/// normalization by 0, so callers treat them as "no seam".
[[nodiscard]] bool hasSeamLength(const SeamTopology &topo);


/// @brief Per-channel D_total of @p texels over precomputed Value pieces (paper eq. 9).
///
/// @param pieces       Output of `seamPieces(..., PieceKind::Value)`.
/// @param totalLength  `totalSeamLength(topo)`; must be > 0.
/// @param texels       `width * height * 3` interleaved values.
[[nodiscard]] std::array<double, 3> integrateEnergy(std::span<const SeamPiece> pieces,
                                                    double                     totalLength,
                                                    std::span<const double>    texels);


/// @brief Gauss-Legendre 3-point rule on [a,b]: nodes and weights (exact to degree 5).
struct GaussRule {
    std::array<double, 3> gamma {};
    std::array<double, 3> weight {};
};


/// @brief Nodes and weights of the 3-point rule on [@p a, @p b].
[[nodiscard]] GaussRule gaussOn(double a, double b);


/// @brief phi(gamma) of @p piece, all eight coefficients.
[[nodiscard]] std::array<double, 8> evaluatePiece(const SeamPiece &piece, double gamma);

} // namespace decimate::seam::detail
