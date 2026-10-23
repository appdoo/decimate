// Runs the Wasm pipeline off the UI thread. The pipeline is one synchronous call, so progress is coarse
// (engine loading -> processing -> done); per-stage timings come back with the result.
// Protocol:  in  {type:'process', id, name, bytes: Uint8Array, options} | {type:'load', id, name, bytes}
//            out {type:'status', id, stage} | {type:'result', id, result} | {type:'error', id, message}
// An exception from process() (e.g. a malformed file) is reported and the worker keeps serving.
import createDecimate from "./decimate.js";

let modPromise = null;
const engine = () => (modPromise ??= createDecimate());

self.onmessage = async ({ data }) => {
    if (data.type !== "process" && data.type !== "load") {
        return;
    }
    const { id } = data;
    try {
        self.postMessage({ type: "status", id, stage: "loading engine" });
        const mod = await engine();
        if (data.type === "load") {
            const mesh = mod.loadMesh(data.name, data.bytes);
            self.postMessage({ type: "result", id, result: mesh }, [mesh.positions.buffer, mesh.indices.buffer]);
            return;
        }
        self.postMessage({ type: "status", id, stage: "processing" });
        const result = mod.process(data.name, data.bytes, data.options);
        const transfer = [result.positions, result.indices, result.uv, result.normals, result.tangents, result.map.rgb]
            .map((a) => a.buffer);
        self.postMessage({ type: "result", id, result }, transfer);
    } catch (e) {
        if (!modPromise || e instanceof TypeError || e instanceof WebAssembly.RuntimeError) {
            modPromise = null; // engine failed to start or trapped (heap state undefined): fresh instance next time
        }
        self.postMessage({ type: "error", id, message: String(e && e.message ? e.message : e) });
    }
};
