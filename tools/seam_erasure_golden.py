#!/usr/bin/env python3
"""Generate the seam-erasure goldens (Spec 17-01) from the MIT reference implementation.

Reference (read-only, never vendored, never copied into src/):
    https://github.com/zfergus/seam-erasure   (MIT, Python)
    pinned commit 9b811ad6a188b28f9a5d656fecb264ec11e45243

Recreate the environment (outside the repo; nothing is installed system-wide):

    git clone https://github.com/zfergus/seam-erasure /somewhere/seam-erasure
    git -C /somewhere/seam-erasure checkout 9b811ad6a188b28f9a5d656fecb264ec11e45243
    python3 -m venv /somewhere/venv
    /somewhere/venv/bin/pip install numpy scipy recordclass pillow tqdm
        # generated with numpy 2.5.4, scipy 1.18.1, recordclass 0.24.1, pillow 13.0.0,
        # tqdm 4.70.1 on CPython 3.14; cvxopt is NOT installed (the reference then falls
        # back to scipy.sparse.linalg.spsolve, which is what produced these goldens)

Run from the repository root:

    /somewhere/venv/bin/python tools/seam_erasure_golden.py \
        --reference /somewhere/seam-erasure [--out tests/fixtures] [--png-dir DIR]

Rows written (tests/fixtures/seam_golden_<name>.txt):
    rot   two charts, 4 triangles, 2-edge seam chain; chart B = A's neighbour rotated 90 deg
    mir   same, chart B mirrored
    glob  rot solved with the reference's do_global=True

REFERENCE DEFECT, CORRECTED HERE (documented, asserted, and reported):
    seam_erasure/seam_intervals.py::compute_edgePair_intervals unions the texel-grid
    crossing parameters of BOTH half-edges, but computes each in that half-edge's OWN
    direction. The two half-edges of a seam pair run in opposite directions, while
    bilerp_energy.E_ab / seam_gradient.E_ab (which swap the second edge's endpoints:
    `((uv0, uv1), (uv1p, uv0p)) = edgePair`) parametrize both by the FIRST edge's
    direction. The second edge's crossings are therefore applied at t instead of 1 - t,
    so on some pieces the bilinear polynomial of the wrong cell is integrated and the
    quadratic form is not the seam energy of the paper (eq. 3-8). This script wraps the
    reference's own compute_edge_intervals with the mirrored parameter (1 - t) for the
    second edge; every other reference function is called unmodified. It asserts, at
    generation time, that the corrected D_before equals an independent dense numerical
    quadrature of the true seam energy, and records the unpatched D_before for the record.
    The golden therefore pins the CORRECT seam energy and the reference solver's
    minimizer of the corrected energy, not the reference's defective quadratic form.
"""

import argparse
import io
import logging
import math
import os
import subprocess
import sys

import numpy as np

PINNED_COMMIT = "9b811ad6a188b28f9a5d656fecb264ec11e45243"
N = 64
SEED = 20241018

# 3D vertices: A0 | seam chain P0 - P1 - P2 | B0 (a flat strip, +z up)
POS = {"A0": (0, 1, 0), "P0": (1, 0, 0), "P1": (1, 1, 0), "P2": (1, 2, 0), "B0": (2, 1, 0)}
POS_ORDER = ["A0", "P0", "P1", "P2", "B0"]
FACES_A = [("A0", "P0", "P1"), ("A0", "P1", "P2")]  # chart A (CCW in 3D)
FACES_B = [("B0", "P1", "P0"), ("B0", "P2", "P1")]  # chart B


def f32(v):
    return float(np.float32(v))


def chart_uvs(variant):
    """Float32 UVs per (chart, vertex). Both charts stay >= 0.07 away from the border."""
    uv_a = {
        "A0": (0.0837, 0.08 + 0.25 * 1),  # off the 45-degree diagonal through texel centres
        "P0": (0.08 + 0.25 * 1, 0.08 + 0.25 * 0),
        "P1": (0.08 + 0.25 * 1 + 0.012, 0.08 + 0.25 * 1),  # bent seam chain
        "P2": (0.08 + 0.25 * 1, 0.08 + 0.25 * 2),
    }
    # B in the plane continuing A's mapping, then rotated / mirrored and moved elsewhere.
    nat = {"P0": (0.33, 0.08), "P1": (0.342, 0.33), "P2": (0.33, 0.58), "B0": (0.58, 0.33)}
    uv_b = {}
    for k, (u, v) in nat.items():
        du, dv = u - 0.33, v - 0.33
        if variant == "rot":
            uv_b[k] = (0.68 - dv * 0.95, 0.70 + du * 0.95)
        else:
            uv_b[k] = (0.68 - du * 0.95, 0.45 + dv * 0.9)
    quant = lambda d: {k: (f32(u), f32(v)) for k, (u, v) in d.items()}
    return quant(uv_a), quant(uv_b)


def build_mesh(variant):
    """Return (positions, uvs, faceV, faceU) with one UV vertex per (chart, vertex)."""
    uv_a, uv_b = chart_uvs(variant)
    pos_id = {k: i for i, k in enumerate(POS_ORDER)}
    uvs, uv_id = [], {}
    face_v, face_u = [], []
    for chart, faces, table in (("A", FACES_A, uv_a), ("B", FACES_B, uv_b)):
        for tri in faces:
            for name in tri:
                key = (chart, name)
                if key not in uv_id:
                    uv_id[key] = len(uvs)
                    uvs.append(table[name])
                face_v.append(pos_id[name])
                face_u.append(uv_id[key])
    positions = [tuple(float(c) for c in POS[k]) for k in POS_ORDER]
    return positions, uvs, face_v, face_u


def make_texture(rng):
    """Smooth + noise, quantized to 8 bit (the fixture input is bytes)."""
    xs = np.arange(N) / (N - 1.0)
    ch = np.arange(3)
    smooth = 0.5 + 0.25 * np.sin(7.0 * xs[None, :, None] + ch[None, None, :]) * np.cos(
        5.0 * xs[:, None, None] + 2.0 * ch[None, None, :])
    tex = 0.5 * smooth + 0.5 * rng.rand(N, N, 3)
    return np.round(np.clip(tex, 0, 1) * 255).astype(np.uint8)


def to_obj_text(positions, uvs, face_v, face_u):
    lines = ["v %.17g %.17g %.17g" % p for p in positions]
    lines += ["vt %.17g %.17g" % uv for uv in uvs]
    for f in range(len(face_v) // 3):
        lines.append("f " + " ".join("%d/%d" % (face_v[3 * f + k] + 1, face_u[3 * f + k] + 1)
                                     for k in range(3)))
    return "\n".join(lines)


def check_texel_centres_generic(uvs, face_u, tol=1e-3):
    """No texel centre may sit on a UV triangle edge: reference and our inside test would differ."""
    for f in range(len(face_u) // 3):
        tri = [(uvs[face_u[3 * f + k]][0] * N - 0.5, uvs[face_u[3 * f + k]][1] * N - 0.5) for k in range(3)]
        for k in range(3):
            (x0, y0), (x1, y1) = tri[k], tri[(k + 1) % 3]
            for py in range(N):
                for px in range(N):
                    d = abs((x1 - x0) * (py - y0) - (y1 - y0) * (px - x0)) / math.hypot(x1 - x0, y1 - y0)
                    t = ((px - x0) * (x1 - x0) + (py - y0) * (y1 - y0)) / ((x1 - x0) ** 2 + (y1 - y0) ** 2)
                    assert not (d < tol and -tol <= t <= 1 + tol), "texel centre on a UV edge"


def bilinear(tex, u, v):
    """GPU convention: uv*N-0.5, clamp-to-edge, row index = v."""
    x = min(max(u * N - 0.5, 0.0), N - 1.0)
    y = min(max(v * N - 0.5, 0.0), N - 1.0)
    x0, y0 = min(int(math.floor(x)), N - 2), min(int(math.floor(y)), N - 2)
    fx, fy = x - x0, y - y0
    return ((1 - fx) * (1 - fy) * tex[y0, x0] + fx * (1 - fy) * tex[y0, x0 + 1]
            + (1 - fx) * fy * tex[y0 + 1, x0] + fx * fy * tex[y0 + 1, x0 + 1])


def true_energy(ref, mesh, tex, samples=6000):
    """Independent dense midpoint quadrature of eq. 9 (no reference interval code)."""
    from seam_erasure.find_seam import find_seam, seam_to_UV
    seam, bnd, fold = find_seam(mesh)
    uv_seam, _, _ = seam_to_UV(mesh, seam, bnd, fold)
    total, length = np.zeros(3), 0.0
    for (e1, _e2), (ua, ub) in zip(seam, uv_seam):
        fi, (i, j) = e1
        p0 = np.array(mesh.v[mesh.f[fi][i].v])
        p1 = np.array(mesh.v[mesh.f[fi][j].v])
        ln = float(np.linalg.norm(p1 - p0))
        (a0, a1), (b1, b0) = ua, ub  # b runs opposite: b0 meets a0
        acc = np.zeros(3)
        for k in range(samples):
            g = (k + 0.5) / samples
            pa = bilinear(tex, a0[0] + g * (a1[0] - a0[0]), a0[1] + g * (a1[1] - a0[1]))
            pb = bilinear(tex, b0[0] + g * (b1[0] - b0[0]), b0[1] + g * (b1[1] - b0[1]))
            acc += (pa - pb) ** 2
        total += ln * acc / samples
        length += ln
    return total / length


def install_interval_correction(ref):
    """Mirror the second half-edge's crossing parameters (t -> 1 - t); see module docstring."""
    from seam_erasure import bilerp_energy, seam_gradient, seam_intervals

    def corrected(edge_pair, width, height):
        out = set(seam_intervals.compute_edge_intervals(edge_pair[0], width, height))
        out |= {round(1.0 - t, 10) for t in seam_intervals.compute_edge_intervals(edge_pair[1], width, height)}
        return sorted(out)

    original = bilerp_energy.compute_edgePair_intervals
    bilerp_energy.compute_edgePair_intervals = corrected
    seam_gradient.compute_edgePair_intervals = corrected
    return original


def restore(original):
    from seam_erasure import bilerp_energy, seam_gradient
    bilerp_energy.compute_edgePair_intervals = original
    seam_gradient.compute_edgePair_intervals = original


def d_per_channel(Q, vec):
    return [float(vec[:, c] @ (Q @ vec[:, c])) for c in range(vec.shape[1])]


def write_fixture(path, name, do_global, data):
    positions, uvs, face_v, face_u = data["mesh"]
    with open(path, "w") as f:
        f.write("# Spec 17-01 seam-erasure golden '%s'. Generated by tools/seam_erasure_golden.py\n" % name)
        f.write("# reference: zfergus/seam-erasure @ %s (interval defect corrected, see script)\n" % PINNED_COMMIT)
        f.write("name %s\nsize %d\nglobal %d\n" % (name, N, int(do_global)))
        f.write("positions %d\n" % len(positions))
        for p in positions:
            f.write("%.9g %.9g %.9g\n" % p)
        f.write("uvs %d\n" % len(uvs))
        for u, v in uvs:
            f.write("%.9g %.9g\n" % (u, v))
        f.write("faces %d\n" % (len(face_v) // 3))
        f.write("faceV " + " ".join(map(str, face_v)) + "\n")
        f.write("faceU " + " ".join(map(str, face_u)) + "\n")
        f.write("counts seams %d boundary %d foldovers %d\n" % data["counts"])
        f.write("dBefore " + " ".join("%.17g" % d for d in data["d_before"]) + "\n")
        f.write("dAfter " + " ".join("%.17g" % d for d in data["d_after"]) + "\n")
        f.write("dBeforeUnpatched " + " ".join("%.17g" % d for d in data["d_before_unpatched"]) + "\n")
        f.write("input  # bytes, row-major, row = v, channel-interleaved\n")
        for row in data["input"]:
            f.write(" ".join(str(int(b)) for b in row.reshape(-1)) + "\n")
        f.write("reference  # float, same layout, unquantized reference solution in [0,1]\n")
        for row in data["output"]:
            f.write(" ".join("%.10g" % v for v in row.reshape(-1)) + "\n")


def run_row(ref, name, variant, do_global, tex_bytes, png_dir):
    from seam_erasure import obj_reader, seam_erasure
    from seam_erasure.find_seam import find_seam

    mesh_data = build_mesh(variant)
    check_texel_centres_generic(mesh_data[1], mesh_data[3])
    mesh = obj_reader.parse_obj(io.StringIO(to_obj_text(*mesh_data)), "golden")
    seam, boundary, foldovers = find_seam(mesh)
    texture = tex_bytes / 255.0
    vec = texture.reshape(N * N, 3)

    unpatched_e = seam_erasure.compute_energies(mesh, texture)
    d_unpatched = d_per_channel(unpatched_e.BLE.Q, vec)

    original = install_interval_correction(ref)
    try:
        energies = seam_erasure.compute_energies(mesh, texture)
        out = seam_erasure.erase_seam(mesh, texture, do_global=do_global)
    finally:
        restore(original)
    d_before = d_per_channel(energies.BLE.Q, vec)
    d_after = d_per_channel(energies.BLE.Q, out)

    truth = true_energy(ref, mesh, texture)
    for got, want in zip(d_before, truth):
        assert abs(got - want) <= 1e-5 * want, "corrected D_before disagrees with numerical truth"
    logging.warning("%s: seams %d boundary %d foldovers %d", name, len(seam), len(boundary), len(foldovers))
    logging.warning("%s: D_before %s", name, d_before)
    logging.warning("%s: D_after  %s", name, d_after)
    logging.warning("%s: D_before with the unpatched reference %s", name, d_unpatched)

    out_img = out.reshape(N, N, 3)
    if png_dir:
        from PIL import Image
        os.makedirs(png_dir, exist_ok=True)
        flip = lambda a: Image.fromarray(a).transpose(Image.FLIP_TOP_BOTTOM)  # reference flips on load
        flip(tex_bytes).save(os.path.join(png_dir, "%s_input.png" % name))
        flip(np.round(np.clip(out_img, 0, 1) * 255).astype(np.uint8)).save(
            os.path.join(png_dir, "%s_reference.png" % name))
    return {"mesh": mesh_data, "counts": (len(seam), len(boundary), len(foldovers)),
            "d_before": d_before, "d_after": d_after, "d_before_unpatched": d_unpatched,
            "input": tex_bytes, "output": out_img}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reference", required=True, help="checkout of zfergus/seam-erasure")
    ap.add_argument("--out", default="tests/fixtures")
    ap.add_argument("--png-dir", default=None, help="also write row-flipped input/reference PNGs here")
    args = ap.parse_args()

    head = subprocess.run(["git", "-C", args.reference, "rev-parse", "HEAD"], capture_output=True,
                          text=True, check=True).stdout.strip()
    if head != PINNED_COMMIT:
        sys.exit("reference is at %s, expected the pinned %s" % (head, PINNED_COMMIT))
    sys.path.insert(0, args.reference)
    logging.basicConfig(level=logging.WARNING, format="%(message)s")

    tex = make_texture(np.random.RandomState(SEED))
    os.makedirs(args.out, exist_ok=True)
    for name, variant, glob in (("rot", "rot", False), ("mir", "mir", False), ("glob", "rot", True)):
        data = run_row(args.reference, name, variant, glob, tex, args.png_dir)
        write_fixture(os.path.join(args.out, "seam_golden_%s.txt" % name), name, glob, data)


if __name__ == "__main__":
    main()
