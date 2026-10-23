// UI orchestration: input, worker, two synced three.js panes, atlas view, export. No mesh algorithms here:
// the pipeline runs in the Wasm core (worker.js); this file only displays its output and draws UV edges.
import * as THREE from "three";
import { OrbitControls } from "./vendor/examples/jsm/controls/OrbitControls.js";
import { GLTFLoader } from "./vendor/examples/jsm/loaders/GLTFLoader.js";
import { toGLB, toObjZip } from "./export.js";

const $ = (id) => document.getElementById(id);
const el = {
    drop: $("drop"), file: $("file"), spot: $("spot"), run: $("run"), ratio: $("ratio"), ratioOut: $("ratioOut"),
    charts: $("charts"), chartsOut: $("chartsOut"), resolution: $("resolution"), space: $("space"),
    dlGlb: $("dlGlb"), dlObj: $("dlObj"), tNormal: $("tNormal"), tWire: $("tWire"), tAtlas: $("tAtlas"),
    stats: $("stats"), toast: $("toast"), progress: $("progress"), atlas: $("atlas"),
    emptyA: $("emptyA"), emptyB: $("emptyB"), tagA: $("tagA"), tagB: $("tagB"),
};

// ---- three.js panes ---------------------------------------------------------------------------------------
const panes = [makePane($("paneA")), makePane($("paneB"))];
let syncing = false;
let drawQueued = false;

function makePane(host) {
    const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
    renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
    host.prepend(renderer.domElement);
    const scene = new THREE.Scene();
    const camera = new THREE.PerspectiveCamera(40, 1, 0.05, 100);
    camera.position.set(0.8, 0.6, 4.6);
    scene.add(camera);
    scene.add(new THREE.HemisphereLight(0xffffff, 0x555566, 0.9));
    const key = new THREE.DirectionalLight(0xffffff, 2.2); // rides on the camera: normal maps stay lit while orbiting
    key.position.set(1, 1.5, 2);
    camera.add(key);
    const holder = new THREE.Group();
    scene.add(holder);
    const controls = new OrbitControls(camera, renderer.domElement);
    controls.addEventListener("change", () => onCameraChange(pane));
    const pane = { host, renderer, scene, camera, holder, controls, object: null };
    new ResizeObserver(() => resize(pane)).observe(host);
    return pane;
}

function resize(p) {
    const w = Math.max(1, p.host.clientWidth), h = Math.max(1, p.host.clientHeight);
    p.renderer.setSize(w, h, false);
    p.camera.aspect = w / h;
    p.camera.updateProjectionMatrix();
    requestDraw();
}

function onCameraChange(src) {
    if (!syncing) {
        syncing = true;
        const dst = panes[0] === src ? panes[1] : panes[0];
        dst.camera.position.copy(src.camera.position);
        dst.camera.quaternion.copy(src.camera.quaternion);
        dst.controls.target.copy(src.controls.target);
        dst.controls.update();
        syncing = false;
    }
    requestDraw();
}

function requestDraw() {
    if (drawQueued) {
        return;
    }
    drawQueued = true;
    requestAnimationFrame(() => {
        drawQueued = false;
        for (const p of panes) {
            p.renderer.render(p.scene, p.camera);
        }
    });
}

/** Show `object` in a pane, centred and scaled to a 2-unit box (each pane fits itself; the camera is shared). */
function setObject(p, object) {
    p.holder.clear();
    p.object = object;
    if (!object) {
        requestDraw();
        return;
    }
    const box = new THREE.Box3().setFromObject(object);
    const size = box.getSize(new THREE.Vector3());
    const s = 2 / Math.max(size.x, size.y, size.z, 1e-9);
    p.holder.scale.setScalar(s);
    p.holder.position.copy(box.getCenter(new THREE.Vector3()).multiplyScalar(-s));
    p.holder.add(object);
    applyViewToggles();
}

function materialsOf(object) {
    const out = [];
    object?.traverse((o) => o.isMesh && out.push(o.material));
    return out.flat();
}

function applyViewToggles() {
    for (const p of panes) {
        for (const m of materialsOf(p.object)) {
            m.wireframe = el.tWire.checked;
            const nm = m.userData.normalMap ?? m.normalMap; // keep the map around while it is toggled off
            m.userData.normalMap = nm;
            const want = el.tNormal.checked ? nm : null;
            if (m.normalMap !== want) {
                m.normalMap = want;
                m.needsUpdate = true;
            }
        }
    }
    updateAtlasVisibility();
    requestDraw();
}

// ---- state ------------------------------------------------------------------------------------------------
let input = null; // {name, bytes}
let current = null; // {result, options, name} of the last successful run
let jobId = 0;
let busyTimer = 0;
const setState = (s) => { document.body.dataset.state = s; };

function toast(msg) {
    el.toast.textContent = msg;
    el.toast.hidden = false;
    clearTimeout(toast.t);
    toast.t = setTimeout(() => { el.toast.hidden = true; }, 9000);
}

function currentOptions() {
    return {
        targetRatio: Number(el.ratio.value) / 100,
        maxCharts: Number(el.charts.value),
        resolution: Number(el.resolution.value),
        tangentSpace: el.space.value === "tangent",
    };
}

// ---- worker -----------------------------------------------------------------------------------------------
let worker = null;
const pending = new Map();

function startWorker() {
    worker = new Worker(new URL("./worker.js", import.meta.url), { type: "module" });
    worker.onmessage = ({ data }) => {
        const job = pending.get(data.id);
        if (!job) {
            return;
        }
        if (data.type === "status") {
            job.stage = data.stage;
            return;
        }
        pending.delete(data.id);
        data.type === "result" ? job.resolve(data.result) : job.reject(new Error(data.message));
    };
    worker.onerror = (e) => {
        // The module failed to load (decimate.js missing) -- fail every waiting job, start fresh next time.
        for (const job of pending.values()) {
            job.reject(new Error("The Wasm engine could not be loaded (decimate.js / decimate.wasm next to index.html?)."));
        }
        pending.clear();
        worker = null;
        e.preventDefault();
    };
}

function runInWorker(type, name, bytes, options) {
    worker ?? startWorker();
    const id = ++jobId;
    return new Promise((resolve, reject) => {
        pending.set(id, { resolve, reject, stage: "" });
        worker.postMessage({ type, id, name, bytes, options });
    });
}

// ---- processing -------------------------------------------------------------------------------------------
async function process() {
    if (!input || document.body.dataset.state === "busy") {
        return;
    }
    const options = currentOptions();
    const t0 = performance.now();
    setState("busy");
    el.progress.hidden = false;
    el.run.disabled = true;
    clearInterval(busyTimer);
    busyTimer = setInterval(() => {
        el.stats.textContent = `Processing ${input.name} (load, reduce, chartify, bake) ... ${((performance.now() - t0) / 1000).toFixed(1)} s`;
    }, 100);
    try {
        const result = await runInWorker("process", input.name, input.bytes.slice(), options);
        current = { result, options, name: input.name };
        showResult(current);
        setState("done");
    } catch (e) {
        toast(e.message);
        setState("error");
        el.stats.textContent = `Failed: ${input.name}`;
    } finally {
        clearInterval(busyTimer);
        el.progress.hidden = true;
        el.run.disabled = false;
        el.run.classList.remove("dirty");
    }
}

function showResult({ result: r, options, name }) {
    setObject(panes[1], buildResultMesh(r, options));
    el.tagB.textContent = "Result";
    el.emptyB.textContent = "";
    drawAtlas(r);
    const s = r.stats, t = s.seconds;
    const f = (n) => n.toLocaleString("en-US");
    el.stats.innerHTML = `<b>${name}</b> &middot; ${f(s.inVerts)} &rarr; <b>${f(s.outVerts)}</b> verts (${(100 * s.outVerts / s.inVerts).toFixed(0)} %)`
        + ` &middot; ${f(r.indices.length / 3)} faces &middot; ${s.charts} charts`
        + (s.erased ? "" : ` &middot; <b>seam erase failed (out of memory?), dilate map</b>`)
        + ` &middot; load ${t.load.toFixed(2)} reduce ${t.reduce.toFixed(2)} chartify ${t.chartify.toFixed(2)} bake ${t.bake.toFixed(2)} s`;
    const glbOk = options.tangentSpace;
    el.dlGlb.disabled = !glbOk;
    el.dlGlb.title = glbOk ? "" : "GLB has no object-space normal map slot: choose Space = Tangent and process again.";
    el.dlObj.disabled = false;
}

/** Non-indexed BufferGeometry from the per-corner arrays. Core convention throughout: row == v, flipY = false. */
function buildResultMesh(r, options) {
    const corners = r.indices.length;
    const pos = new Float32Array(3 * corners);
    for (let c = 0; c < corners; c++) {
        const p = 3 * r.indices[c];
        pos[3 * c] = r.positions[p];
        pos[3 * c + 1] = r.positions[p + 1];
        pos[3 * c + 2] = r.positions[p + 2];
    }
    const g = new THREE.BufferGeometry();
    g.setAttribute("position", new THREE.BufferAttribute(pos, 3));
    g.setAttribute("normal", new THREE.BufferAttribute(r.normals, 3));
    g.setAttribute("uv", new THREE.BufferAttribute(r.uv, 2));
    if (options.tangentSpace) {
        g.setAttribute("tangent", new THREE.BufferAttribute(r.tangents, 4));
    }
    const w = r.map.width, h = r.map.height;
    const rgba = new Uint8Array(w * h * 4);
    for (let i = 0, n = w * h; i < n; i++) {
        rgba.set(r.map.rgb.subarray(3 * i, 3 * i + 3), 4 * i);
        rgba[4 * i + 3] = 255;
    }
    const tex = new THREE.DataTexture(rgba, w, h, THREE.RGBAFormat);
    tex.flipY = false; // data row 0 is v = 0, the same as glTF; see export.js header
    tex.generateMipmaps = true;
    tex.minFilter = THREE.LinearMipmapLinearFilter;
    tex.magFilter = THREE.LinearFilter;
    tex.anisotropy = 4;
    tex.needsUpdate = true;
    const mat = new THREE.MeshStandardMaterial({ color: 0xcfd2d8, roughness: 0.55, metalness: 0, normalMap: tex });
    mat.normalMapType = options.tangentSpace ? THREE.TangentSpaceNormalMap : THREE.ObjectSpaceNormalMap;
    return new THREE.Mesh(g, mat);
}

// ---- original pane ----------------------------------------------------------------------------------------
// The input as the Wasm core loads it (OBJ, Inventor, ...): one loader for both panes, smooth vertex normals.
async function showOriginal({ name, bytes }) {
    setObject(panes[0], null);
    el.tagA.textContent = "Original";
    el.emptyA.textContent = "Loading ...";
    try {
        const m = await runInWorker("load", name, bytes.slice());
        const g = new THREE.BufferGeometry();
        g.setAttribute("position", new THREE.BufferAttribute(m.positions, 3));
        g.setIndex(new THREE.BufferAttribute(m.indices, 1));
        g.computeVertexNormals();
        const gray = new THREE.MeshStandardMaterial({ color: 0xcfd2d8, roughness: 0.6, metalness: 0 });
        setObject(panes[0], new THREE.Mesh(g, gray));
        el.emptyA.textContent = "";
        el.tagA.textContent = `Original · ${(m.indices.length / 3).toLocaleString("en-US")} triangles`;
    } catch (e) {
        el.emptyA.textContent = `Could not preview this file: ${e.message}`;
    }
}

// ---- atlas view -------------------------------------------------------------------------------------------
// The map as a 2D image with the UV triangle edges on top: chart outlines (seams / boundaries) bold, interior
// edges faint. Row == v, so uv * size is the canvas position directly.
function drawAtlas(r) {
    const scale = Math.max(1, Math.ceil(1024 / r.map.width));
    const W = r.map.width * scale, H = r.map.height * scale;
    el.atlas.width = W;
    el.atlas.height = H;
    const ctx = el.atlas.getContext("2d");
    const img = new ImageData(r.map.width, r.map.height);
    for (let i = 0, n = r.map.width * r.map.height; i < n; i++) {
        img.data.set(r.map.rgb.subarray(3 * i, 3 * i + 3), 4 * i);
        img.data[4 * i + 3] = 255;
    }
    const tmp = new OffscreenCanvas(r.map.width, r.map.height);
    tmp.getContext("2d").putImageData(img, 0, 0);
    ctx.imageSmoothingEnabled = false;
    ctx.drawImage(tmp, 0, 0, W, H);

    const { inner, outline } = uvEdges(r);
    const stroke = (segs, style, width) => {
        ctx.beginPath();
        for (const [u0, v0, u1, v1] of segs) {
            ctx.moveTo(u0 * W, v0 * H);
            ctx.lineTo(u1 * W, v1 * H);
        }
        ctx.strokeStyle = style;
        ctx.lineWidth = width;
        ctx.stroke();
    };
    stroke(inner, "rgba(255,255,255,0.22)", 1);
    stroke(outline, "rgba(255,210,60,0.95)", 2 * scale / 2 + 1);
}

/** UV edges split into chart outlines (boundary, or uv differs across the 3D edge) and interior edges. */
function uvEdges(r) {
    const V = r.positions.length / 3;
    const seen = new Map();
    const inner = [], outline = new Set();
    const corner = (c) => [r.uv[2 * c], r.uv[2 * c + 1]];
    for (let f = 0; f < r.indices.length / 3; f++) {
        for (let e = 0; e < 3; e++) {
            const ca = 3 * f + e, cb = 3 * f + (e + 1) % 3;
            const a = r.indices[ca], b = r.indices[cb];
            const [lo, hi] = a < b ? [ca, cb] : [cb, ca];
            const seg = [...corner(lo), ...corner(hi)];
            const key = Math.min(a, b) * V + Math.max(a, b);
            const prev = seen.get(key);
            if (prev === undefined) {
                seen.set(key, seg);
                outline.add(seg); // boundary until a second face claims the edge
                continue;
            }
            const same = prev.every((x, i) => x === seg[i]);
            outline.delete(prev);
            if (same) {
                inner.push(seg);
            } else {
                outline.add(prev).add(seg);
            }
        }
    }
    return { inner, outline: [...outline] };
}

function updateAtlasVisibility() {
    const show = el.tAtlas.checked && current !== null && !panes[1].glbView;
    el.atlas.hidden = !show;
    panes[1].renderer.domElement.style.visibility = show ? "hidden" : "visible";
}

// ---- export -----------------------------------------------------------------------------------------------
const baseName = () => (current?.name ?? "mesh").replace(/\.[^.]+$/, "") + "_decimated";

function download(bytes, fileName, type) {
    const url = URL.createObjectURL(new Blob([bytes], { type }));
    const a = document.createElement("a");
    a.href = url;
    a.download = fileName;
    document.body.append(a);
    a.click();
    a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 10000);
}

async function exportGlb() {
    const bytes = await toGLB(current.result);
    download(bytes, `${baseName()}.glb`, "model/gltf-binary");
    return bytes;
}

async function exportObj() {
    const bytes = await toObjZip(current.result, { tangentSpace: current.options.tangentSpace, name: baseName() });
    download(bytes, `${baseName()}.zip`, "application/zip");
    return bytes;
}

// ---- GLB viewing (round-trip check of an exported file) --------------------------------------------------
function viewGlb(name, bytes) {
    new GLTFLoader().parse(bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength), "", (gltf) => {
        panes[1].glbView = true;
        current = null;
        setObject(panes[1], gltf.scene);
        el.tagB.textContent = `GLB · ${name}`;
        el.emptyB.textContent = "";
        let tris = 0;
        gltf.scene.traverse((o) => { if (o.isMesh) { tris += (o.geometry.index?.count ?? o.geometry.attributes.position.count) / 3; } });
        el.stats.textContent = `Viewing ${name} · ${tris.toLocaleString("en-US")} triangles`;
        el.dlGlb.disabled = el.dlObj.disabled = true;
        setState("glb");
    }, (e) => toast(`Cannot read ${name}: ${e.message ?? e}`));
}

// ---- input ------------------------------------------------------------------------------------------------
async function loadBytes(name, bytes) {
    if (/\.glb$/i.test(name)) {
        viewGlb(name, bytes);
        return;
    }
    panes[1].glbView = false;
    input = { name, bytes };
    el.run.disabled = false;
    await showOriginal(input);
    await process();
}

async function loadFile(file) {
    await loadBytes(file.name, new Uint8Array(await file.arrayBuffer()));
}

el.file.addEventListener("change", () => el.file.files[0] && loadFile(el.file.files[0]));
el.drop.addEventListener("keydown", (e) => { if (e.key === "Enter" || e.key === " ") { e.preventDefault(); el.file.click(); } });
el.spot.addEventListener("click", async () => {
    try {
        const res = await fetch("spot.obj");
        if (!res.ok) {
            throw new Error(`spot.obj: HTTP ${res.status}`);
        }
        await loadBytes("spot.obj", new Uint8Array(await res.arrayBuffer()));
    } catch (e) {
        toast(e.message);
    }
});
el.run.addEventListener("click", process);
el.dlGlb.addEventListener("click", () => exportGlb().catch((e) => toast(e.message)));
el.dlObj.addEventListener("click", () => exportObj().catch((e) => toast(e.message)));

let dragDepth = 0;
const hasFiles = (e) => [...(e.dataTransfer?.types ?? [])].includes("Files");
window.addEventListener("dragenter", (e) => { if (hasFiles(e)) { dragDepth++; document.body.classList.add("dragging"); } });
window.addEventListener("dragleave", () => { if (--dragDepth <= 0) { dragDepth = 0; document.body.classList.remove("dragging"); } });
window.addEventListener("dragover", (e) => { if (hasFiles(e)) { e.preventDefault(); } });
window.addEventListener("drop", (e) => {
    e.preventDefault();
    dragDepth = 0;
    document.body.classList.remove("dragging");
    const f = e.dataTransfer?.files?.[0];
    if (f) {
        loadFile(f);
    }
});

// ---- controls ---------------------------------------------------------------------------------------------
const markDirty = () => { if (input) { el.run.classList.add("dirty"); } };
el.ratio.addEventListener("input", () => { el.ratioOut.textContent = `${el.ratio.value} %`; markDirty(); });
el.charts.addEventListener("input", () => { el.chartsOut.textContent = el.charts.value; markDirty(); });
el.resolution.addEventListener("change", markDirty);
el.space.addEventListener("change", markDirty);
for (const t of [el.tNormal, el.tWire, el.tAtlas]) {
    t.addEventListener("change", applyViewToggles);
}
el.emptyB.textContent = "The result appears here.";
requestDraw();

// Test hook (web/test and the W5 driver): the same code paths the buttons use, minus the file dialogs.
window.__app = {
    loadBytes,
    glb: async () => exportGlb(),
    obj: async () => exportObj(),
    get state() { return document.body.dataset.state; },
};
