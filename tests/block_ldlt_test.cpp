// This file is part of the decimate library; see decimate.h for version/license details.
/// @file block_ldlt_test.cpp
/// @brief Spec 31: the 3x3-block LDL^T of the framed erase.
///
/// a [S] spot 256 framed K: block solution vs Eigen SimplicialLDLT <= 1e-8 relative, residual no worse
/// b [S] a random 3x3-block SPD matrix (seed 31, 50 blocks) against a dense Cholesky reference
/// c [S] an indefinite diagonal block makes factorize() refuse (erase then retries once with a ridge)
/// d [I] seam metric: tangent erase of spot 256 / cow 256 within 0.005 deg of the Spec 33 whole-domain values (spot 1024 hidden)

#include "block_ldlt.h"
#include "erase_domain.h"
#include "erase_system.h"
#include "model_fixture.h"
#include "seam_metric.h"
#include "seam_pieces.h"

#include <Eigen/Sparse>
#include <catch2/catch_test_macros.hpp>
#include <decimate/bake/normal_bake.h>
#include <decimate/bake/normals.h>
#include <decimate/pipeline.h>
#include <decimate/reduce/reducer.h>
#include <decimate/seam/erase.h>
#include <decimate/seam/seam_topology.h>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <random>

using namespace decimate;
using namespace decimate::seam::detail;

namespace {

/// The viewer's pipeline on `high`: reduce to 24 % of the vertices, 16 charts.
struct Fixture {
    MeshData                 high;
    Atlas                    atlas;
    std::vector<math::Vec3f> vn;
    seam::SeamTopology       topo;
};


Fixture makeFixture(const std::string &path) {
    Fixture f;
    f.high = load(path);
    const MeshData low = reduce(f.high, static_cast<int>(static_cast<float>(f.high.positions.size()) * 0.24f));
    ChartOptions   co;
    co.maxCharts = 16;
    f.atlas = chartify(low, co);
    f.vn = bake::computeVertexNormals(f.atlas.mesh);
    f.topo = seam::findSeams(seam::weld(f.atlas.mesh, f.atlas.uv));
    return f;
}


const Fixture &spot() {
    static const Fixture f = makeFixture(std::string(MESH_MODEL_DIR) + "/spot.obj");
    return f;
}


const Fixture &cow() {
    static const Fixture f = makeFixture(modelOrSkip("cow.obj"));
    return f;
}


/// The framed erase system of `f` at `res`, exactly as `seam::erase` assembles it.
EraseSystem<3> framedSystem(const Fixture &f, int res) {
    BakeOptions o;
    o.space = BakeSpace::Tangent;
    o.gutter = Gutter::Dilate;
    o.resolution = res;
    const bake::NormalMap map = bakeAtlas(f.atlas, f.high, o).map;

    const auto          tm = seam::weld(f.atlas.mesh, f.atlas.uv);
    const FrameField    frames { f.atlas.tangents, f.vn };
    const TexelGrid     grid { map.width, map.height };
    std::vector<double> texels;
    for (const bake::Rgb &p : map.pixels) {
        for (const int b : { p.r, p.g, p.b }) {
            texels.push_back(static_cast<double>(b) / 127.0 - 1.0);
        }
    }
    const auto   topo = seam::findSeams(tm);
    const auto   value = seamPieces(tm, topo, grid, PieceKind::Value, &frames);
    const auto   grad = seamPieces(tm, topo, grid, PieceKind::Gradient);
    const auto   domain = buildEraseDomain(tm, topo, grid);
    EraseWeights weights;
    weights.change = 1e4;
    return buildEraseSystem<3>({ domain, grid, texels, value, grad, totalSeamLength(topo), weights });
}


/// Both triangles of the block matrix as an Eigen sparse matrix.
Eigen::SparseMatrix<double> toEigen(const BlockMatrix<3> &m) {
    std::vector<Eigen::Triplet<double>> t;
    for (int j = 0; j < m.n; ++j) {
        for (int e = 0; e < 9; ++e) {
            t.emplace_back(3 * j + e / 3, 3 * j + e % 3, m.diag[9 * static_cast<std::size_t>(j) + static_cast<std::size_t>(e)]);
        }
        for (int s = m.start[static_cast<std::size_t>(j)]; s < m.start[static_cast<std::size_t>(j) + 1]; ++s) {
            for (int e = 0; e < 9; ++e) {
                const double v = m.off[9 * static_cast<std::size_t>(s) + static_cast<std::size_t>(e)];
                t.emplace_back(3 * m.row[static_cast<std::size_t>(s)] + e / 3, 3 * j + e % 3, v);
                t.emplace_back(3 * j + e % 3, 3 * m.row[static_cast<std::size_t>(s)] + e / 3, v);
            }
        }
    }
    Eigen::SparseMatrix<double> k(3 * m.n, 3 * m.n);
    k.setFromTriplets(t.begin(), t.end());
    return k;
}


seammetric::Stats tangentErased(const Fixture &f, int res) {
    BakeOptions o;
    o.space = BakeSpace::Tangent;
    o.gutter = Gutter::Erased;
    o.resolution = res;
    return seammetric::seamStats(f.atlas, f.vn, bakeAtlas(f.atlas, f.high, o).map, true, f.topo);
}

} // namespace


TEST_CASE("Spec 31: block LDLT matches Eigen on the spot 256 framed system", "[seam][spec31][S]") {
    const EraseSystem<3> sys = framedSystem(spot(), 256);
    REQUIRE(sys.n > 1000);

    BlockLdlt<3> ldlt(sys.blocks);
    REQUIRE(ldlt.factorize(sys.blocks));
    std::vector<double> x(sys.rhs.size());
    ldlt.solve(sys.rhs, x);

    const Eigen::SparseMatrix<double>                  k = toEigen(sys.blocks);
    const Eigen::Map<const Eigen::VectorXd>            b(sys.rhs.data(), static_cast<Eigen::Index>(sys.rhs.size()));
    const Eigen::Map<const Eigen::VectorXd>            xb(x.data(), static_cast<Eigen::Index>(x.size()));
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> eigen(k);
    REQUIRE(eigen.info() == Eigen::Success);
    const Eigen::VectorXd xe = eigen.solve(b);

    CHECK((xb - xe).norm() / xe.norm() <= 1e-8);
    CHECK((k * xb - b).norm() / b.norm() <= (k * xe - b).norm() / b.norm() * 1.01 + 1e-15);
}


TEST_CASE("Spec 31: a random SPD block matrix solves like the dense reference", "[seam][spec31][S]") {
    constexpr int                          n = 50;
    std::mt19937                           rng(31);
    std::uniform_real_distribution<double> u(-1.0, 1.0);

    // A random sparse symmetric pattern, off-diagonal blocks in [-1, 1], diagonals strictly dominant => SPD.
    BlockPattern pattern(n);
    for (int a = 0; a < n; ++a) {
        for (int b = 0; b < a; ++b) {
            if (std::uniform_int_distribution<int>(0, 5)(rng) == 0) {
                pattern.couple(a, b);
            }
        }
    }
    BlockMatrix<3>      m = pattern.build<3>();
    std::vector<double> dense(static_cast<std::size_t>(9 * n * n), 0.0);
    const auto          at = [&](int r, int c) -> double          &{
        return dense[static_cast<std::size_t>(r * 3 * n + c)];
    };
    for (int j = 0; j < n; ++j) {
        for (int s = m.start[static_cast<std::size_t>(j)]; s < m.start[static_cast<std::size_t>(j) + 1]; ++s) {
            for (int e = 0; e < 9; ++e) {
                const double v = u(rng);
                m.off[9 * static_cast<std::size_t>(s) + static_cast<std::size_t>(e)] = v;
                at(3 * m.row[static_cast<std::size_t>(s)] + e / 3, 3 * j + e % 3) = v;
                at(3 * j + e % 3, 3 * m.row[static_cast<std::size_t>(s)] + e / 3) = v;
            }
        }
    }
    for (int j = 0; j < n; ++j) {
        for (int r = 0; r < 3; ++r) {
            at(3 * j + r, 3 * j + r) = 40.0 + u(rng);
            m.diag[9 * static_cast<std::size_t>(j) + static_cast<std::size_t>(4 * r)] = at(3 * j + r, 3 * j + r);
        }
    }
    std::vector<double> b(static_cast<std::size_t>(3 * n));
    for (double &v : b) {
        v = u(rng);
    }

    // Dense reference: Gaussian elimination on the SPD matrix (no pivoting needed).
    const int           dim = 3 * n;
    std::vector<double> a = dense;
    std::vector<double> ref = b;
    for (int k = 0; k < dim; ++k) {
        for (int i = k + 1; i < dim; ++i) {
            const double f = a[static_cast<std::size_t>(i * dim + k)] / a[static_cast<std::size_t>(k * dim + k)];
            for (int j = k; j < dim; ++j) {
                a[static_cast<std::size_t>(i * dim + j)] -= f * a[static_cast<std::size_t>(k * dim + j)];
            }
            ref[static_cast<std::size_t>(i)] -= f * ref[static_cast<std::size_t>(k)];
        }
    }
    for (int i = dim - 1; i >= 0; --i) {
        for (int j = i + 1; j < dim; ++j) {
            ref[static_cast<std::size_t>(i)] -= a[static_cast<std::size_t>(i * dim + j)] * ref[static_cast<std::size_t>(j)];
        }
        ref[static_cast<std::size_t>(i)] /= a[static_cast<std::size_t>(i * dim + i)];
    }

    BlockLdlt<3> ldlt(m);
    REQUIRE(ldlt.factorize(m));
    std::vector<double> x(b.size());
    ldlt.solve(b, x);
    double err = 0.0;
    double norm = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        err += (x[i] - ref[i]) * (x[i] - ref[i]);
        norm += ref[i] * ref[i];
    }
    CHECK(std::sqrt(err / norm) <= 1e-12);
}


TEST_CASE("Spec 31: an indefinite diagonal block makes the block factorization refuse", "[seam][spec31][S]") {
    BlockPattern pattern(3);
    pattern.couple(0, 1);
    pattern.couple(1, 2);
    BlockMatrix<3> m = pattern.build<3>();
    for (int j = 0; j < 3; ++j) {
        for (int r = 0; r < 3; ++r) {
            m.diag[static_cast<std::size_t>(9 * j + 4 * r)] = 4.0;
        }
    }
    m.diag[9 * 1 + 4] = -4.0; // block 1: diag(4, -4, 4)

    BlockLdlt<3> ldlt(m);
    CHECK_FALSE(ldlt.factorize(m));
}


TEST_CASE("Spec 31: framed seam metric is unchanged by the block solver (spot 256, cow 256)", "[seam][spec31][I]") {
    // Spec 33 (whole domain, no band) values: spot 256 0.133 / 0.50, cow 256 0.204 / 0.71.
    const auto s = tangentErased(spot(), 256);
    CHECK(std::fabs(s.mean - 0.132) <= 0.005);
    CHECK(std::fabs(s.max - 0.50) <= 0.005);
    const auto c = tangentErased(cow(), 256);
    CHECK(std::fabs(c.mean - 0.204) <= 0.005);
    CHECK(std::fabs(c.max - 0.71) <= 0.005);
}


TEST_CASE("Spec 31: framed seam metric spot 1024", "[.][seam][spec31][I]") {
    const auto s = tangentErased(spot(), 1024);
    CHECK(std::fabs(s.mean - 0.133) <= 0.005);
    CHECK(std::fabs(s.max - 0.52) <= 0.005);
}
