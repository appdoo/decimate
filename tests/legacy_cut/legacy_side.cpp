// Legacy side of the chart_legacy_golden A/B test (Spec 19-01).
//
// Drives the UNMODIFIED 2003 CChartification (legacy2003/chartification/*.cpp, compiled verbatim by
// the `legacy_chart_cut` object library) and recovers the per-face chart assignment.
//
// Chartify() deletes its chartMap on exit (legacy2003/chartification/cchartification.cpp:112), so
// faceToChart is rebuilt from the public CChart::Faces()/Vertices(): each chart face's global vertex
// triple (Vertex::index) is looked up, sorted, in a triple -> input-face map.  A duplicate triangle
// would make that map ambiguous, so it is a hard error.
//
// clang-tidy: legacy type/stub idioms (Vec3 P/F/N naming, trivial stub dtors) are intentional here.
// This TU is built as C++14 with -fms-extensions -fdelayed-template-parsing (see tests/CMakeLists.txt).
#include "math/algorithm.h"
#include "chartification/chartification.h"
#include "chartification/cconformalmapper.h"
#include "chartification/catlaspacker.h"

#include <algorithm>
#include <array>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <vector>

/// Symbols Chartify() never reaches (only TextureVisualization uses them).  Stubbed so the cone
/// links without the conformal mapper / atlas packer; they abort() if ever called.
void setup_timer() {
}


void start_timer() {
}


double stop_timer() {
    return 0.0;
}


void show_progress() {
}


/// Link stubs for the mapper and packer (never called by Chartify, see above).
namespace chartification {
CConformalMapper::CConformalMapper(Vec3 *) {
}
CConformalMapper::~CConformalMapper() = default;


void CConformalMapper::Parametrize(CChartSet *) {
    std::abort();
}


CAtlasPacker::CAtlasPacker(Vec3 *, int) {
}


CAtlasPacker::~CAtlasPacker() = default;


int CAtlasPacker::Pack(CChartSet *) {
    std::abort();
}

} // namespace chartification

using namespace chartification;

/// Runs the 2003 chartification with cap @p k.  Returns the chart count, or -1 on a duplicate
/// triangle.  @p faceToChart (nf ints) gets a chart ordinal per face (-1 = face never reached).
/// @p localVerts (k+1 ints, may be null) gets CChart::VertexCount() per chart in iteration order.
// NOLINTNEXTLINE(readability-function-size,readability-identifier-naming,bugprone-implicit-widening-of-multiplication-result): legacy-typed harness
int legacyChartify(int np, const float *pts, int nf, const int *idx, int k,
                   int *faceToChart, int *localVerts) {
    // NOLINTBEGIN(readability-identifier-naming,bugprone-implicit-widening-of-multiplication-result)

    // Repack the flat arrays into the legacy Vec3 / Face types.
    std::vector<Vec3> P(static_cast<size_t>(np));
    for (int i = 0; i < np; ++i) {
        P[i].x = pts[3 * i];
        P[i].y = pts[3 * i + 1];
        P[i].z = pts[3 * i + 2];
    }

    // The sorted vertex triple identifies a face for the faceToChart rebuild below.
    std::vector<Face>                 F(static_cast<size_t>(nf));
    std::map<std::array<int, 3>, int> key;
    for (int i = 0; i < nf; ++i) {
        F[i].a = idx[3 * i];
        F[i].b = idx[3 * i + 1];
        F[i].c = idx[3 * i + 2];
        std::array<int, 3> t = { F[i].a, F[i].b, F[i].c };
        std::sort(t.begin(), t.end());
        if (!key.emplace(t, i).second) {
            return -1; // duplicate triangle: the triple lookup would be ambiguous
        }
    }

    CChartification cf(np, P.data(), nf, F.data());
    // CChartSet reads the raw face normals only in Create(seed) (legacy2003/chartification/cchartset.cpp:57-64).
    std::vector<Vec3> N(static_cast<size_t>(nf));
    for (int i = 0; i < nf; ++i) {
        N[i] = normal_raw(P.data(), &F[i]);
    }
    CChartSet *cs = new CChartSet(P.data(), N.data(), nullptr);
    cf.Chartify(k, cs);

    // Chartify() has freed its chartMap, so the assignment is recovered from the chart faces.
    std::fill(faceToChart, faceToChart + nf, -1);
    int id = 0;
    for (CChartSet::iterator it = cs->Begin(); it != cs->End(); ++it, ++id) {
        CChart *c = *it;
        Vertex *V = c->Vertices();
        Face   *CF = c->Faces();
        for (int j = 0; j < c->FaceCount(); ++j) {
            std::array<int, 3> t = { V[CF[j].a].index, V[CF[j].b].index, V[CF[j].c].index };
            std::sort(t.begin(), t.end());
            const auto f = key.find(t);
            if (f != key.end()) {
                faceToChart[f->second] = id;
            }
        }
        if (localVerts != nullptr) {
            localVerts[id] = c->VertexCount();
        }
    }

    // NOLINTEND(readability-identifier-naming,bugprone-implicit-widening-of-multiplication-result)
    const int count = cs->Size();
    delete cs;
    return count;
}
