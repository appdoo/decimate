// Export of a `process()` result: GLB 2.0, OBJ + MTL + PNG (as one store-only ZIP). Pure data in, bytes out;
// no DOM, so the node tests (web/test/w3_w4.mjs) run the exact code the page runs.
//
// ---------------------------------------------------------------------------------------------------------
// UV / ROW CONVENTION -- the one place where it is decided.
//   Core:  normal-map texel (x, y) is the pixel at uv = ((x+0.5)/N, (y+0.5)/N); image row index == v.
//          Row 0 (the first row of map.rgb, the top of the PNG) is v = 0.
//   glTF:  texture coordinate (0,0) is the TOP-left of the image, v grows downward -> identical to the core.
//          No flip. (three.js: GLTFLoader and the viewer panes use texture.flipY = false for the same reason.)
//   OBJ:   vt (0,0) is the BOTTOM-left, v grows upward (OBJLoader/Blender flip the image on load) -> vt = 1 - v.
// Everything below goes through gltfV()/objV(); nothing else in web/ may touch v.
//
//   Normal-map green channel (tangent space). The core's bitangent is B = w*cross(N,T) and it points along +v
//   (measured on spot: 3072/3072 faces) = DOWN the image; green encodes the component along B. glTF and OBJ/
//   OpenGL both define green = +Y = UP the image. Both exporters therefore convert once, in glConvention():
//   green -> 254-green and tangent w -> -w (n = xT + yB + zN is unchanged, so a viewer that uses our TANGENT
//   renders identically; one that recomputes tangents from the uv now also gets the right sign).
//   The in-page viewer shows the raw core data (it uses the core's own tangents), so it needs no conversion.
// ---------------------------------------------------------------------------------------------------------
export const gltfV = (v) => v;
export const objV = (v) => 1 - v;

// ---- welding ---------------------------------------------------------------------------------------------
// `process()` returns an indexed position set but per-corner uv/normal/tangent. Weld corners that agree on
// (position index, uv, normal, tangent) into one vertex. Exact float equality: corners of one chart vertex are
// bit-identical, and a seam must stay split.
/** Tangent-space result in the glTF/OpenGL green-up convention (see header). Object-space maps pass through. */
export function glConvention(r, tangentSpace = true) {
    if (!tangentSpace) {
        return r;
    }
    const rgb = Uint8Array.from(r.map.rgb);
    for (let i = 1; i < rgb.length; i += 3) {
        rgb[i] = Math.max(0, 254 - rgb[i]);
    }
    const tangents = Float32Array.from(r.tangents);
    for (let i = 3; i < tangents.length; i += 4) {
        tangents[i] = -tangents[i];
    }
    return { ...r, tangents, map: { ...r.map, rgb } };
}

const DEGENERATE_TANGENT = [1, 0, 0, 1];

/** @returns {{positions: Float32Array, normals: Float32Array, uv: Float32Array, tangents: Float32Array, indices: Uint32Array}} */
export function weld(r) {
    const corners = r.indices.length;
    const remap = new Map();
    const src = []; // corner index of the first occurrence of every welded vertex
    const indices = new Uint32Array(corners);
    for (let c = 0; c < corners; c++) {
        const key = [r.indices[c],
            r.uv[2 * c], r.uv[2 * c + 1],
            r.normals[3 * c], r.normals[3 * c + 1], r.normals[3 * c + 2],
            r.tangents[4 * c], r.tangents[4 * c + 1], r.tangents[4 * c + 2], r.tangents[4 * c + 3]].join(",");
        let id = remap.get(key);
        if (id === undefined) {
            id = src.length;
            remap.set(key, id);
            src.push(c);
        }
        indices[c] = id;
    }
    const n = src.length;
    const out = {
        positions: new Float32Array(3 * n), normals: new Float32Array(3 * n),
        uv: new Float32Array(2 * n), tangents: new Float32Array(4 * n), indices,
    };
    for (let i = 0; i < n; i++) {
        const c = src[i];
        const p = 3 * r.indices[c];
        out.positions.set(r.positions.subarray(p, p + 3), 3 * i);
        out.normals.set(r.normals.subarray(3 * c, 3 * c + 3), 3 * i);
        out.uv[2 * i] = r.uv[2 * c];
        out.uv[2 * i + 1] = r.uv[2 * c + 1];
        const t = r.tangents.subarray(4 * c, 4 * c + 4);
        // Corners of degenerate faces carry (0,0,0,1): no basis. glTF wants a unit xyz and w = +-1.
        const ok = Math.hypot(t[0], t[1], t[2]) > 0.5;
        out.tangents.set(ok ? t : DEGENERATE_TANGENT, 4 * i);
    }
    return out;
}

// ---- PNG ------------------------------------------------------------------------------------------------
const CRC_TABLE = (() => {
    const t = new Uint32Array(256);
    for (let n = 0; n < 256; n++) {
        let c = n;
        for (let k = 0; k < 8; k++) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1;
        t[n] = c >>> 0;
    }
    return t;
})();

export function crc32(bytes, crc = 0) {
    let c = ~crc >>> 0;
    for (let i = 0; i < bytes.length; i++) c = CRC_TABLE[(c ^ bytes[i]) & 0xff] ^ (c >>> 8);
    return ~c >>> 0;
}

async function zlibDeflate(bytes) {
    const stream = new Blob([bytes]).stream().pipeThrough(new CompressionStream("deflate"));
    return new Uint8Array(await new Response(stream).arrayBuffer());
}

/** 8-bit RGB PNG. Shared by browser and node (CompressionStream exists in both), so the test covers it. */
export async function encodePNG(width, height, rgb) {
    const row = width * 3;
    const raw = new Uint8Array((row + 1) * height); // filter byte 0 per scanline
    for (let y = 0; y < height; y++) raw.set(rgb.subarray(y * row, (y + 1) * row), y * (row + 1) + 1);
    const idat = await zlibDeflate(raw);
    const chunk = (type, data) => {
        const b = new Uint8Array(12 + data.length);
        const dv = new DataView(b.buffer);
        dv.setUint32(0, data.length);
        for (let i = 0; i < 4; i++) b[4 + i] = type.charCodeAt(i);
        b.set(data, 8);
        dv.setUint32(8 + data.length, crc32(b.subarray(4, 8 + data.length)));
        return b;
    };
    const ihdr = new Uint8Array(13);
    const dv = new DataView(ihdr.buffer);
    dv.setUint32(0, width);
    dv.setUint32(4, height);
    ihdr.set([8, 2, 0, 0, 0], 8); // 8 bit, colour type 2 (RGB), deflate, no filter scheme, no interlace
    return concat([Uint8Array.of(137, 80, 78, 71, 13, 10, 26, 10),
        chunk("IHDR", ihdr), chunk("IDAT", idat), chunk("IEND", new Uint8Array(0))]);
}

function concat(parts) {
    const out = new Uint8Array(parts.reduce((s, p) => s + p.length, 0));
    let o = 0;
    for (const p of parts) {
        out.set(p, o);
        o += p.length;
    }
    return out;
}

// ---- GLB ------------------------------------------------------------------------------------------------
const pad4 = (n) => (n + 3) & ~3;

/**
 * GLB 2.0 of a tangent-space result: POSITION, NORMAL, TEXCOORD_0, TANGENT, one material with the PNG as
 * normalTexture. Object-space maps have no glTF slot, so the caller must not offer this for them.
 * @param raw `process()` result (tangent space, core convention); converted by glConvention() here.
 * @returns {Promise<Uint8Array>}
 */
export async function toGLB(raw) {
    const r = glConvention(raw);
    const png = await encodePNG(r.map.width, r.map.height, r.map.rgb);
    const w = weld(r);
    const uv = new Float32Array(w.uv.length);
    for (let i = 0; i < uv.length; i += 2) {
        uv[i] = w.uv[i];
        uv[i + 1] = gltfV(w.uv[i + 1]);
    }
    const blobs = [
        { name: "POSITION", data: w.positions, target: 34962, comp: 5126, type: "VEC3" },
        { name: "NORMAL", data: w.normals, target: 34962, comp: 5126, type: "VEC3" },
        { name: "TEXCOORD_0", data: uv, target: 34962, comp: 5126, type: "VEC2" },
        { name: "TANGENT", data: w.tangents, target: 34962, comp: 5126, type: "VEC4" },
        { name: "indices", data: w.indices, target: 34963, comp: 5125, type: "SCALAR" },
    ];
    const bufferViews = [];
    const accessors = [];
    const parts = [];
    let offset = 0;
    for (const b of blobs) {
        const bytes = new Uint8Array(b.data.buffer, b.data.byteOffset, b.data.byteLength);
        bufferViews.push({ buffer: 0, byteOffset: offset, byteLength: bytes.length, target: b.target });
        const nComp = { SCALAR: 1, VEC2: 2, VEC3: 3, VEC4: 4 }[b.type];
        const acc = { bufferView: bufferViews.length - 1, componentType: b.comp, count: b.data.length / nComp, type: b.type };
        if (b.name === "POSITION") {
            acc.min = [Infinity, Infinity, Infinity];
            acc.max = [-Infinity, -Infinity, -Infinity];
            for (let i = 0; i < b.data.length; i += 3) {
                for (let k = 0; k < 3; k++) {
                    acc.min[k] = Math.min(acc.min[k], b.data[i + k]);
                    acc.max[k] = Math.max(acc.max[k], b.data[i + k]);
                }
            }
        }
        accessors.push(acc);
        parts.push(bytes);
        offset += pad4(bytes.length);
        parts.push(new Uint8Array(pad4(bytes.length) - bytes.length));
    }
    bufferViews.push({ buffer: 0, byteOffset: offset, byteLength: png.length }); // the image: no target
    parts.push(png, new Uint8Array(pad4(png.length) - png.length));
    offset += pad4(png.length);
    const bin = concat(parts);

    const json = {
        asset: { version: "2.0", generator: "decimate web tool" },
        scene: 0,
        scenes: [{ nodes: [0] }],
        nodes: [{ mesh: 0 }],
        meshes: [{ primitives: [{ attributes: { POSITION: 0, NORMAL: 1, TEXCOORD_0: 2, TANGENT: 3 }, indices: 4, material: 0 }] }],
        materials: [{ pbrMetallicRoughness: { baseColorFactor: [0.8, 0.8, 0.8, 1], metallicFactor: 0, roughnessFactor: 0.6 },
            normalTexture: { index: 0 } }],
        textures: [{ sampler: 0, source: 0 }],
        samplers: [{ magFilter: 9729, minFilter: 9987, wrapS: 33071, wrapT: 33071 }],
        images: [{ bufferView: 5, mimeType: "image/png" }],
        buffers: [{ byteLength: bin.length }],
        bufferViews,
        accessors,
    };
    const jsonBytes = new TextEncoder().encode(JSON.stringify(json));
    const jsonPadded = new Uint8Array(pad4(jsonBytes.length)).fill(0x20);
    jsonPadded.set(jsonBytes);
    const total = 12 + 8 + jsonPadded.length + 8 + bin.length;
    const glb = new Uint8Array(total);
    const dv = new DataView(glb.buffer);
    dv.setUint32(0, 0x46546c67, true); // "glTF"
    dv.setUint32(4, 2, true);
    dv.setUint32(8, total, true);
    dv.setUint32(12, jsonPadded.length, true);
    dv.setUint32(16, 0x4e4f534a, true); // "JSON"
    glb.set(jsonPadded, 20);
    const o = 20 + jsonPadded.length;
    dv.setUint32(o, bin.length, true);
    dv.setUint32(o + 4, 0x004e4942, true); // "BIN\0"
    glb.set(bin, o + 8);
    return glb;
}

// ---- OBJ + MTL ------------------------------------------------------------------------------------------
/** @returns {{obj: string, mtl: string}}; the MTL points at `pngName`. */
export function toOBJ(r, { tangentSpace, name = "decimated", pngName = "normal.png" }) {
    const w = weld(r);
    const nv = w.positions.length / 3;
    const L = [`# decimate web tool -- ${nv} vertices, ${w.indices.length / 3} faces`, `mtllib ${name}.mtl`, `o ${name}`, "usemtl normalmap"];
    for (let i = 0; i < nv; i++) L.push(`v ${w.positions[3 * i]} ${w.positions[3 * i + 1]} ${w.positions[3 * i + 2]}`);
    for (let i = 0; i < nv; i++) L.push(`vt ${w.uv[2 * i]} ${objV(w.uv[2 * i + 1])}`);
    for (let i = 0; i < nv; i++) L.push(`vn ${w.normals[3 * i]} ${w.normals[3 * i + 1]} ${w.normals[3 * i + 2]}`);
    for (let f = 0; f < w.indices.length; f += 3) {
        const a = w.indices[f] + 1, b = w.indices[f + 1] + 1, c = w.indices[f + 2] + 1;
        L.push(`f ${a}/${a}/${a} ${b}/${b}/${b} ${c}/${c}/${c}`);
    }
    const mtl = [`# ${tangentSpace ? "tangent-space" : "object-space"} normal map`,
        "newmtl normalmap", "Ka 0 0 0", "Kd 0.8 0.8 0.8", "Ks 0.2 0.2 0.2", "Ns 30",
        `map_Bump ${pngName}`, `norm ${pngName}`].join("\n") + "\n";
    return { obj: L.join("\n") + "\n", mtl };
}

// ---- ZIP (store only) ----------------------------------------------------------------------------------
/** @param files {name: string, data: Uint8Array}[] */
export function toZip(files) {
    const enc = new TextEncoder();
    const local = [];
    const central = [];
    let offset = 0;
    for (const f of files) {
        const name = enc.encode(f.name);
        const crc = crc32(f.data);
        const lh = new Uint8Array(30 + name.length);
        const l = new DataView(lh.buffer);
        l.setUint32(0, 0x04034b50, true);
        l.setUint16(4, 20, true);
        l.setUint32(14, crc, true);
        l.setUint32(18, f.data.length, true);
        l.setUint32(22, f.data.length, true);
        l.setUint16(26, name.length, true);
        lh.set(name, 30);
        const ch = new Uint8Array(46 + name.length);
        const c = new DataView(ch.buffer);
        c.setUint32(0, 0x02014b50, true);
        c.setUint16(4, 20, true);
        c.setUint16(6, 20, true);
        c.setUint32(16, crc, true);
        c.setUint32(20, f.data.length, true);
        c.setUint32(24, f.data.length, true);
        c.setUint16(28, name.length, true);
        c.setUint32(42, offset, true);
        ch.set(name, 46);
        local.push(lh, f.data);
        central.push(ch);
        offset += lh.length + f.data.length;
    }
    const cd = concat(central);
    const end = new Uint8Array(22);
    const e = new DataView(end.buffer);
    e.setUint32(0, 0x06054b50, true);
    e.setUint16(8, files.length, true);
    e.setUint16(10, files.length, true);
    e.setUint32(12, cd.length, true);
    e.setUint32(16, offset, true);
    return concat([...local, cd, end]);
}

/** OBJ + MTL + PNG as one ZIP. */
export async function toObjZip(raw, { tangentSpace, name = "decimated" }) {
    const r = glConvention(raw, tangentSpace);
    const png = await encodePNG(r.map.width, r.map.height, r.map.rgb);
    const { obj, mtl } = toOBJ(r, { tangentSpace, name, pngName: "normal.png" });
    const enc = new TextEncoder();
    return toZip([{ name: `${name}.obj`, data: enc.encode(obj) }, { name: `${name}.mtl`, data: enc.encode(mtl) },
        { name: "normal.png", data: png }]);
}
