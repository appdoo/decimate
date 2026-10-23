// W3 + W4 (spec 34). usage: node w3_w4.mjs <decimate.js> <spot.obj>
// W3: toGLB output passes the Khronos glTF validator with 0 errors; every TANGENT w is +-1.
// W4: toOBJ output re-parsed by three.js OBJLoader: same face count, uv present and == 1 - v of the core.
// Exit 0 pass, 1 fail, 77 skipped (web/test/node_modules missing: run `npm ci` in web/test).
import { readFileSync, existsSync, writeFileSync, mkdtempSync } from "node:fs";
import { execFileSync } from "node:child_process";
import { tmpdir } from "node:os";
import { join, dirname } from "node:path";
import { fileURLToPath, pathToFileURL } from "node:url";
import { inflateSync } from "node:zlib";

const here = dirname(fileURLToPath(import.meta.url));
if (!existsSync(join(here, "node_modules", "gltf-validator")) || !existsSync(join(here, "node_modules", "three"))) {
    console.log("SKIP W3/W4: run `npm ci` in web/test");
    process.exit(77);
}
const [, , glue, model] = process.argv;
let failed = 0;
const check = (ok, msg) => { console.log((ok ? "PASS " : "FAIL ") + msg); if (!ok) failed++; };

const validator = await import("gltf-validator");
const { OBJLoader } = await import("three/examples/jsm/loaders/OBJLoader.js");
const exp = await import(pathToFileURL(join(here, "..", "export.js")));

// The vendored copies the page loads must be the pinned release the test runs against.
for (const f of ["examples/jsm/loaders/OBJLoader.js", "examples/jsm/loaders/GLTFLoader.js", "examples/jsm/controls/OrbitControls.js"]) {
    check(readFileSync(join(here, "node_modules/three", f)).equals(readFileSync(join(here, "..", "vendor", f))), `vendored ${f} == three@pinned`);
}

const { default: createDecimate } = await import(pathToFileURL(glue));
const mod = await createDecimate();
const r = mod.process("spot.obj", new Uint8Array(readFileSync(model)), { targetRatio: 0.25, maxCharts: 16, resolution: 512, tangentSpace: true });
const faces = r.indices.length / 3;
console.log(`input: ${r.stats.outVerts} verts, ${faces} faces, ${r.stats.charts} charts`);

// PNG: decode our own encoder's output and compare byte for byte.
const png = await exp.encodePNG(r.map.width, r.map.height, r.map.rgb);
{
    const dv = new DataView(png.buffer, png.byteOffset);
    check(png[1] === 0x50 && dv.getUint32(16) === r.map.width && dv.getUint32(20) === r.map.height, "PNG header size");
    const idatLen = dv.getUint32(33);
    const raw = inflateSync(png.subarray(41, 41 + idatLen));
    const row = r.map.width * 3;
    let same = raw.length === (row + 1) * r.map.height;
    for (let y = 0; same && y < r.map.height; y++) {
        same = raw[y * (row + 1)] === 0 && Buffer.compare(raw.subarray(y * (row + 1) + 1, (y + 1) * (row + 1)), r.map.rgb.subarray(y * row, (y + 1) * row)) === 0;
    }
    check(same, "PNG decodes back to map.rgb exactly");
}

// ---- W3
const glb = await exp.toGLB(r);
const report = await validator.validateBytes(glb, { uri: "decimated.glb" });
const m = report.issues.messages.filter((x) => x.severity === 0).map((x) => `${x.code} @ ${x.pointer}`);
console.log(`validator: ${report.issues.numErrors} errors, ${report.issues.numWarnings} warnings, ${report.issues.numInfos} infos`);
for (const x of report.issues.messages.filter((x) => x.severity <= 1)) console.log(`  [${x.severity}] ${x.code} ${x.pointer} ${x.message}`);
check(report.issues.numErrors === 0, `W3 glTF validator: 0 errors (${m.join("; ")})`);
check(report.info.totalTriangleCount === faces, `W3 validator triangle count ${report.info.totalTriangleCount} == ${faces}`);
{
    // Read TANGENT straight back out of the GLB.
    const dv = new DataView(glb.buffer, glb.byteOffset);
    const jl = dv.getUint32(12, true);
    const json = JSON.parse(new TextDecoder().decode(glb.subarray(20, 20 + jl)));
    const bin = glb.subarray(20 + jl + 8);
    const acc = json.accessors[json.meshes[0].primitives[0].attributes.TANGENT];
    const bv = json.bufferViews[acc.bufferView];
    const t = new Float32Array(bin.buffer.slice(bin.byteOffset + bv.byteOffset, bin.byteOffset + bv.byteOffset + bv.byteLength));
    let bad = 0, nonUnit = 0;
    for (let i = 0; i < acc.count; i++) {
        if (Math.abs(t[4 * i + 3]) !== 1) bad++;
        if (Math.abs(Math.hypot(t[4 * i], t[4 * i + 1], t[4 * i + 2]) - 1) > 1e-3) nonUnit++;
    }
    check(acc.count > 0 && bad === 0, `W3 TANGENT w = +-1 for all ${acc.count} vertices (bad ${bad})`);
    check(nonUnit === 0, `W3 TANGENT xyz unit length (non-unit ${nonUnit})`);
    // Convention (export.js header): exported B = w*cross(N,T) must point UP the image = against +v, per face;
    // and the embedded PNG's green must be 254 - core green.
    const view = (a, T) => { const bv = json.bufferViews[a.bufferView]; return new T(bin.buffer.slice(bin.byteOffset + bv.byteOffset, bin.byteOffset + bv.byteOffset + bv.byteLength)); };
    const at = json.meshes[0].primitives[0].attributes;
    const P = view(json.accessors[at.POSITION], Float32Array), N = view(json.accessors[at.NORMAL], Float32Array);
    const U = view(json.accessors[at.TEXCOORD_0], Float32Array), I = view(json.accessors[json.meshes[0].primitives[0].indices], Uint32Array);
    let up = 0, down = 0;
    for (let f = 0; f < I.length; f += 3) {
        const [a, b, c] = [I[f], I[f + 1], I[f + 2]];
        const e1 = [0, 1, 2].map((k) => P[3 * b + k] - P[3 * a + k]), e2 = [0, 1, 2].map((k) => P[3 * c + k] - P[3 * a + k]);
        const du1 = U[2 * b] - U[2 * a], dv1 = U[2 * b + 1] - U[2 * a + 1], du2 = U[2 * c] - U[2 * a], dv2 = U[2 * c + 1] - U[2 * a + 1];
        const det = du1 * dv2 - du2 * dv1;
        if (Math.abs(det) < 1e-12) continue;
        const dPdv = e1.map((x, k) => (-du2 * x + du1 * e2[k]) / det);
        const n = [N[3 * a], N[3 * a + 1], N[3 * a + 2]], tt = [t[4 * a], t[4 * a + 1], t[4 * a + 2]], w = t[4 * a + 3];
        const cr = [n[1] * tt[2] - n[2] * tt[1], n[2] * tt[0] - n[0] * tt[2], n[0] * tt[1] - n[1] * tt[0]];
        ((cr[0] * dPdv[0] + cr[1] * dPdv[1] + cr[2] * dPdv[2]) * w < 0 ? up++ : down++);
    }
    check(up > 0 && down === 0, `W3 convention: B = w*cross(N,T) points up the image on all faces (up ${up}, down ${down})`);
    {
        const img = json.images[0], ibv = json.bufferViews[img.bufferView];
        const p = bin.subarray(ibv.byteOffset, ibv.byteOffset + ibv.byteLength);
        const idat = new DataView(p.buffer, p.byteOffset).getUint32(33);
        const raw = inflateSync(p.subarray(41, 41 + idat));
        const row = r.map.width * 3;
        let okG = true;
        for (let y = 0; y < r.map.height && okG; y++) for (let x = 0; x < r.map.width; x++) {
            const o = y * (row + 1) + 1 + 3 * x, g = r.map.rgb[y * row + 3 * x + 1];
            if (raw[o + 1] !== Math.max(0, 254 - g) || raw[o] !== r.map.rgb[y * row + 3 * x]) { okG = false; break; }
        }
        check(okG, "W3 embedded PNG: R unchanged, G = 254 - core G");
    }
    check(json.materials[0].normalTexture.index === 0 && json.images[0].mimeType === "image/png", "W3 normalTexture -> embedded PNG");
}

// ---- W4
{
    const { obj, mtl } = exp.toOBJ(r, { tangentSpace: true });
    const group = new OBJLoader().parse(obj);
    let tris = 0, uvCount = 0, vnCount = 0, firstMesh;
    group.traverse((o) => {
        if (!o.isMesh) return;
        firstMesh ??= o;
        tris += o.geometry.attributes.position.count / 3;
        uvCount += o.geometry.attributes.uv?.count ?? 0;
        vnCount += o.geometry.attributes.normal?.count ?? 0;
    });
    check(tris === faces, `W4 OBJLoader face count ${tris} == ${faces}`);
    check(uvCount > 0 && uvCount === 3 * faces && vnCount === 3 * faces, `W4 OBJLoader uv count ${uvCount} > 0, normals ${vnCount}`);
    // corner 0 of face 0: OBJ vt must be (u, 1 - v) of the core uv (the single flip, export.js).
    const uv = firstMesh.geometry.attributes.uv;
    check(Math.abs(uv.getX(0) - r.uv[0]) < 1e-5 && Math.abs(uv.getY(0) - (1 - r.uv[1])) < 1e-5, "W4 vt = (u, 1-v) of the core uv");
    check(/map_Bump normal\.png/.test(mtl) && /^norm normal\.png/m.test(mtl), "W4 MTL references the PNG");
    const zip = await exp.toObjZip(r, { tangentSpace: true });
    const dir = mkdtempSync(join(tmpdir(), "w4-"));
    writeFileSync(join(dir, "o.zip"), zip);
    try {
        const out = execFileSync("unzip", ["-t", join(dir, "o.zip")]).toString();
        check(/No errors detected/.test(out) && /decimated\.obj/.test(out) && /normal\.png/.test(out), "W4 zip passes `unzip -t`");
    } catch (e) { console.log("SKIP zip test (unzip missing)"); }
}
process.exit(failed ? 1 : 0);
