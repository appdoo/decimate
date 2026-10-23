// This file is part of the decimate library; see decimate.h for version/license details.
/// @file tangents.cpp
/// @brief Thin wrapper over MikkTSpace (mikktspace.c, zlib licence): the only caller of genTangSpaceDefault.

#include <decimate/bake/tangents.h>

#include <mikktspace.h>

namespace decimate {

namespace {

/// What the MikkTSpace callbacks read and write: the inputs by reference, one output slot per corner.
struct TangentJob {
    const MeshData                 &mesh;
    const std::vector<math::Vec2f> &uv;
    const std::vector<math::Vec3f> &normals;
    std::vector<math::Vec4f>       &out;
};

} // namespace


static const TangentJob &job(const SMikkTSpaceContext *ctx) {
    return *static_cast<const TangentJob *>(ctx->m_pUserData);
}


static int getNumFaces(const SMikkTSpaceContext *ctx) {
    return static_cast<int>(job(ctx).mesh.indices.size() / 3);
}


static int getNumVerticesOfFace(const SMikkTSpaceContext *, int) {
    return 3;
}


static void getPosition(const SMikkTSpaceContext *ctx, float pos[], int face, int vert) {
    const math::Vec3f &p = job(ctx).mesh.positions[job(ctx).mesh.indices[3 * face + vert]];
    pos[0] = p.x;
    pos[1] = p.y;
    pos[2] = p.z;
}


static void getNormal(const SMikkTSpaceContext *ctx, float nrm[], int face, int vert) {
    const math::Vec3f &n = job(ctx).normals[job(ctx).mesh.indices[3 * face + vert]];
    nrm[0] = n.x;
    nrm[1] = n.y;
    nrm[2] = n.z;
}


static void getTexCoord(const SMikkTSpaceContext *ctx, float uv[], int face, int vert) {
    const math::Vec2f &t = job(ctx).uv[3 * face + vert];
    uv[0] = t.x;
    uv[1] = t.y;
}


static void setTSpaceBasic(const SMikkTSpaceContext *ctx, const float tan[], float sign, int face, int vert) {
    job(ctx).out[3 * face + vert] = { tan[0], tan[1], tan[2], sign };
}


std::vector<math::Vec4f> tangents(const MeshData &mesh, const std::vector<math::Vec2f> &cornerUV,
                                  const std::vector<math::Vec3f> &vertexNormals) {
    // MikkTSpace leaves the corners of faces it cannot handle untouched: (0,0,0,1) marks "no basis".
    std::vector<math::Vec4f> out(mesh.indices.size(), math::Vec4f { 0.0f, 0.0f, 0.0f, 1.0f });

    SMikkTSpaceInterface iface {};
    iface.m_getNumFaces = getNumFaces;
    iface.m_getNumVerticesOfFace = getNumVerticesOfFace;
    iface.m_getPosition = getPosition;
    iface.m_getNormal = getNormal;
    iface.m_getTexCoord = getTexCoord;
    iface.m_setTSpaceBasic = setTSpaceBasic;

    TangentJob         j { mesh, cornerUV, vertexNormals, out };
    SMikkTSpaceContext ctx { &iface, &j };
    genTangSpaceDefault(&ctx);
    return out;
}

} // namespace decimate
