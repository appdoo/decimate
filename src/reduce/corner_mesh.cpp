// This file is part of the decimate library; see decimate.h for version/license details.
/// @file corner_mesh.cpp
/// @brief 8-byte corner-table build + guarded edge-collapse operator (see
///        corner_mesh.h). The collapse/ring-walk are the author's 2003
///        `CHalfedgeCollapse` (chedgecollapse.cpp:122-178), ported from raw
///        `Halfedge*` pointers to integer handles, plus a link/flip guard.

#include <decimate/reduce/corner_mesh.h>

#include <algorithm>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace decimate::corner {

using math::Vec3f;

namespace {

// Repoint every half-edge ending at the collapsed vertex to p1, both ring
// directions (boundary-aware). Port of chedgecollapse.cpp:122-151.
void replacePoint(std::vector<He8>& he, HeIndex start, HeIndex p1) {
    HeIndex loop = start;
    do {
        he[loop].end = p1;
        loop         = he[nextA(loop)].sym;  // loop->next->sym
    } while (loop != kNoHE && loop != start);
    if (loop == start) {
        return;
    }
    loop = start;
    do {
        he[loop].end    = p1;
        const HeIndex s = he[loop].sym;
        loop            = (s != kNoHE) ? nextA(nextA(s)) : kNoHE;  // loop->sym->next->next
    } while (loop != kNoHE && loop != start);
}

// Re-sym the two surviving edges of the face owning h, then tombstone the face.
// Port of the "Fläche ausklinken" relink, chedgecollapse.cpp:162-176.
void unhookFace(std::vector<He8>& he, HeIndex h, std::vector<char>& deadFace) {
    const HeIndex sn  = he[nextA(h)].sym;
    const HeIndex snn = he[nextA(nextA(h))].sym;
    if (sn != kNoHE) {
        he[sn].sym = snn;
    }
    if (snn != kNoHE) {
        he[snn].sym = sn;
    }
    deadFace[h / 3u] = 1;
}

bool isNeighborOf(const std::vector<He8>& he, HeIndex vIncoming, HeIndex w) {
    bool found = false;
    forEachIncident(he, vIncoming, [&](HeIndex loop) {
        if (originOf(he, loop) == w) {
            found = true;
        }
    });
    return found;
}

bool faceHasVerts(const std::vector<He8>& he, HeIndex vIncoming, HeIndex y, HeIndex z) {
    bool found = false;
    forEachIncident(he, vIncoming, [&](HeIndex loop) {
        const HeIndex f0 = (loop / 3u) * 3u;
        const HeIndex a  = he[f0].end;
        const HeIndex b  = he[f0 + 1].end;
        const HeIndex c  = he[f0 + 2].end;
        if ((a == y || b == y || c == y) && (a == z || b == z || c == z)) {
            found = true;
        }
    });
    return found;
}

// Link condition (Dey-Edelsbrunner-Guha): vertex part (no fin) + edge part (no
// tetra) + shared-apex guard.
bool linkOk(const std::vector<He8>& he, HeIndex h, HeIndex p1, HeIndex p2) {
    const HeIndex vl = he[nextA(h)].end;
    const HeIndex hs = he[h].sym;
    const HeIndex vr = (hs != kNoHE) ? he[nextA(hs)].end : kNoHE;
    if (hs != kNoHE && vl == vr) {
        return false;
    }
    const HeIndex p1in = nextA(nextA(h));
    bool          bad  = false;
    forEachIncident(he, h, [&](HeIndex loop) {
        const HeIndex w = originOf(he, loop);
        if (w != p1 && w != vl && w != vr && isNeighborOf(he, p1in, w)) {
            bad = true;
        }
    });
    if (bad) {
        return false;
    }
    if (hs != kNoHE && faceHasVerts(he, p1in, vl, vr)) {
        return false;  // tetra-like
    }
    return true;
}

// Geometric quality gate. A collapse is rejected if, for any surviving face
// incident to either endpoint, moving the endpoint(s) to `target` would:
//   (a) reverse the face (winding flip / degeneracy), OR
//   (b) rotate the face normal past kMaxNormalDevCos (a sharp CREASE — the
//       dark slit / red gash the bare reversal test missed: two neighbours each
//       rotate <90deg but end up >120deg apart), OR
//   (c) make the face a sliver thinner than kMinAspect (needle triangle).
// (b)+(c) are the quality strengthening over the original reversal-only test.
constexpr float kMaxNormalDevCos = 0.30F;  // reject rotation > ~72.5 deg per face
constexpr float kMinAspect       = 0.05F;  // reject aspect worse than ~1:20

bool degradesShape(const std::vector<He8>& he, const std::vector<Vec3f>& pos, HeIndex h, HeIndex p1,
                   HeIndex p2, const Vec3f& target) {
    bool       bad   = false;
    const auto check = [&](HeIndex loop) {
        const HeIndex f0 = (loop / 3u) * 3u;
        const HeIndex va = he[f0].end;
        const HeIndex vb = he[f0 + 1].end;
        const HeIndex vc = he[f0 + 2].end;
        if ((va == p1 || vb == p1 || vc == p1) && (va == p2 || vb == p2 || vc == p2)) {
            return;  // a face on the collapsed edge disappears — don't test it
        }
        const auto posOf = [&](HeIndex v, bool moved) {
            return (moved && (v == p1 || v == p2)) ? target : pos[v];
        };
        const Vec3f before =
            (posOf(vb, false) - posOf(va, false)).cross(posOf(vc, false) - posOf(va, false));
        const Vec3f a     = posOf(va, true);
        const Vec3f eb    = posOf(vb, true) - a;
        const Vec3f ec    = posOf(vc, true) - a;
        const Vec3f after = eb.cross(ec);
        const float lenB  = before.length();
        const float lenA  = after.length();
        if (lenA < 1e-12F || before.dot(after) < 0.0F) {
            bad = true;  // (a) reversal or degeneracy
            return;
        }
        if (lenB > 1e-12F && before.dot(after) < kMaxNormalDevCos * lenB * lenA) {
            bad = true;  // (b) sharp crease — normal rotated too far
            return;
        }
        const float longest = std::max({eb.length(), ec.length(), (ec - eb).length()});
        if (longest > 1e-12F && (lenA / (longest * longest)) < kMinAspect) {
            bad = true;  // (c) sliver: minAltitude/longest = |after|/longest^2
        }
    };
    forEachIncident(he, h, check);
    forEachIncident(he, nextA(nextA(h)), check);
    return bad;
}

}  // namespace

std::vector<He8> build(const MeshData& mesh) {
    const int        numFaces = static_cast<int>(mesh.indices.size()) / 3;
    const int        numHE    = numFaces * 3;
    std::vector<He8> he(static_cast<std::size_t>(numHE));
    for (int f = 0; f < numFaces; ++f) {
        const std::size_t base = static_cast<std::size_t>(f) * 3;
        // Reversed corner-end winding (c, b, a) — the inverse of faceVertices().
        he[base].end     = static_cast<HeIndex>(mesh.indices[base + 2]);  // ac -> c
        he[base + 1].end = static_cast<HeIndex>(mesh.indices[base + 1]);  // cb -> b
        he[base + 2].end = static_cast<HeIndex>(mesh.indices[base]);      // ba -> a
    }
    std::unordered_map<std::uint64_t, HeIndex> edgeMap;
    edgeMap.reserve(static_cast<std::size_t>(numHE));
    for (int i = 0; i < numHE; ++i) {
        const HeIndex       h     = static_cast<HeIndex>(i);
        const HeIndex       b     = tipOf(he, h);
        const HeIndex       a     = originOf(he, h);
        const std::uint64_t keyBA = (static_cast<std::uint64_t>(b) << 32) | a;
        const auto          it    = edgeMap.find(keyBA);
        if (it != edgeMap.end()) {
            he[h].sym          = it->second;
            he[it->second].sym = h;
            edgeMap.erase(it);
        } else {
            edgeMap.emplace((static_cast<std::uint64_t>(a) << 32) | b, h);
        }
    }
    return he;
}

void collapse(std::vector<He8>& he, HeIndex h, HeIndex p1, std::vector<char>& deadFace) {
    replacePoint(he, h, p1);
    unhookFace(he, h, deadFace);
    const HeIndex hs = he[h].sym;
    if (hs != kNoHE) {
        unhookFace(he, hs, deadFace);
    }
}

bool collapseLegal(const std::vector<He8>& he, const std::vector<Vec3f>& pos, HeIndex h, HeIndex p1,
                   HeIndex p2, const Vec3f& target) {
    return linkOk(he, h, p1, p2) && !degradesShape(he, pos, h, p1, p2, target);
}

MeshData reconstruct(const std::vector<He8>& he, const std::vector<Vec3f>& pos,
                     const std::vector<char>& deadFace, const std::vector<char>& deadVert,
                     int numFaces) {
    MeshData         out;
    std::vector<int> newId(pos.size(), -1);
    int              counter = 0;
    for (std::size_t v = 0; v < pos.size(); ++v) {
        if (deadVert[v] == 0) {
            newId[v] = counter++;
            out.positions.push_back(pos[v]);
        }
    }
    out.indices.reserve(static_cast<std::size_t>(numFaces) * 3);
    for (int f = 0; f < numFaces; ++f) {
        if (deadFace[static_cast<std::size_t>(f)] != 0) {
            continue;
        }
        const HeIndex                base = static_cast<HeIndex>(f) * 3u;
        const std::array<HeIndex, 3> v    = faceVertices(he, base);  // CCW, normal-consistent
        const int                    a    = newId[v[0]];
        const int                    b    = newId[v[1]];
        const int                    c    = newId[v[2]];
        if (a < 0 || b < 0 || c < 0 || a == b || b == c || a == c) {
            continue;
        }
        out.indices.push_back(a);
        out.indices.push_back(b);
        out.indices.push_back(c);
    }
    return out;
}

}  // namespace decimate::corner
