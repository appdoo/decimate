// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

/// @file chart_cut.h
/// @brief Private: cut graphs inside a chart (Spec 20 phase 2, Gu, Gortler, Hoppe 2002 "Geometry Images", section 3).
///
/// A cut is a set of mesh edges inside one chart.  Both half-edges of every cut edge go into
/// `Chart::seamHalfEdges`, so `flatten()` does not merge the corners across them and the chart opens
/// along the cut.  Edges are named by their lower half-edge id.

#include <decimate/chart/chart.h>
#include <decimate/chart/parametrize.h>
#include <decimate/mesh/half_edge_mesh.h>
#include <decimate/types/mesh_data.h>
#include <vector>

namespace decimate::chart {

/// True when @p chartId has no boundary at all: every edge of its faces is shared by two of its faces
/// and no seam opens it, so `flatten()` would build a closed local mesh that cannot be flattened or packed.
[[nodiscard]] bool isClosedChart(const HalfEdgeMesh &mesh, const ChartSet &charts, int chartId);


/// Cut graph that opens a closed chart into a disk with at least three boundary vertices.
///
/// Genus 0 (Euler characteristic 2): the shortest edge path (3D length, Dijkstra, ties to the lower
/// vertex id) between the two mutually farthest vertices (double sweep), extended by one edge when it
/// is a single edge.  Genus > 0: a dual spanning tree grown from the chart's seed face (lowest face id
/// first) is removed; the edges it does not cross form the cut graph, pruned of dangling edges
/// (degree-1 end vertices, repeatedly).
///
/// @pre `isClosedChart(mesh, charts, chartId)`.
/// @returns Cut edges, ascending; empty when the chart has fewer than 4 vertices.
[[nodiscard]] std::vector<int> closedChartCut(const HalfEdgeMesh &mesh, const MeshData &data, const ChartSet &charts,
                                              int chartId);


/// Up to @p maxPaths vertex-disjoint shortest interior-edge paths to the boundary (3D edge length, Dijkstra,
/// ties to the lower vertex id), one per face of @p faces taken in order (worst first): from the corner of
/// the face farthest from the boundary.  A face whose corners all lie on the boundary, or whose path
/// shares a vertex with an earlier one, is skipped.  The boundary is every vertex on a chart outline edge
/// or an existing seam; interior edges only: both faces in the chart and not yet a seam.
///
/// @returns The edges of the paths, each corner first; empty when no face yields a path.
[[nodiscard]] std::vector<int> augmentationPaths(const HalfEdgeMesh &mesh, const MeshData &data, const ChartSet &charts,
                                                 int chartId, const std::vector<int> &faces, int maxPaths);


/// Add @p edges (both halves of each) to @p chart's seams, keeping the list ascending and unique.
void addCutEdges(const HalfEdgeMesh &mesh, Chart &chart, const std::vector<int> &edges);

} // namespace decimate::chart
