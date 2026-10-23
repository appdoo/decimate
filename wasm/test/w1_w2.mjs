// W1 + W2 (spec 34). usage: node w1_w2.mjs <decimate.js> <spot.obj> <wasm_ref binary>
// Exit 0 = pass, 1 = fail, 77 = skipped (native reference missing).
import { readFileSync, existsSync, mkdtempSync } from "node:fs";
import { execFileSync } from "node:child_process";
import { tmpdir } from "node:os";
import { join } from "node:path";
import { pathToFileURL } from "node:url";

const [, , glue, model, ref] = process.argv;
let failed = 0;
const check = (ok, msg) => { console.log((ok ? "PASS " : "FAIL ") + msg); if (!ok) failed++; };

const { default: createDecimate } = await import(pathToFileURL(glue));
const mod = await createDecimate();
const bytes = new Uint8Array(readFileSync(model));
const defaults = { targetRatio: 0.25, maxCharts: 16, resolution: 512, tangentSpace: false };
const out = mod.process("spot.obj", bytes, defaults);
const { stats, map } = out;

// W1
const ratio = stats.outVerts / stats.inVerts;
check(Math.abs(ratio - 0.25) <= 0.05, `W1 outVerts ${stats.outVerts}/${stats.inVerts} = ${ratio.toFixed(3)} within 25% +-5%`);
check(out.positions.length === 3 * stats.outVerts, "W1 positions length 3V");
check(out.indices.length % 3 === 0 && out.indices.every((i) => i < stats.outVerts), "W1 indices in range");
const corners = out.indices.length;
check(out.uv.length === 2 * corners && out.normals.length === 3 * corners && out.tangents.length === 4 * corners,
      "W1 per-corner uv/normals/tangents lengths");
check(out.uv.every((u) => u >= 0 && u <= 1), "W1 uv in [0,1]");
check(map.rgb.length === map.width * map.height * 3 && map.width === 512 && map.height === 512, `W1 map ${map.width}x${map.height} rgb bytes`);
check(map.rgb.some((b) => b !== map.rgb[0]), "W1 map not uniform");
console.log("wasm seconds", JSON.stringify(stats.seconds), "charts", stats.charts);

// W2
if (!existsSync(ref)) {
    console.log(`SKIP W2: native reference not built (${ref})`);
    process.exit(failed ? 1 : 77);
}
const dump = join(mkdtempSync(join(tmpdir(), "w2-")), "ref.rgb");
const nat = JSON.parse(execFileSync(ref, [model, dump]).toString());
check(nat.outVerts === stats.outVerts && nat.faces === out.indices.length / 3,
      `W2 counts native ${nat.outVerts}v/${nat.faces}f vs wasm ${stats.outVerts}v/${out.indices.length / 3}f`);
const nrgb = new Uint8Array(readFileSync(dump));
check(nrgb.length === map.rgb.length, "W2 map size equal");
let maxDiff = 0;
for (let i = 0; i < Math.min(nrgb.length, map.rgb.length); i++) maxDiff = Math.max(maxDiff, Math.abs(nrgb[i] - map.rgb[i]));
check(maxDiff <= 2, `W2 map max byte diff = ${maxDiff} (<= 2)`);
console.log("native seconds", JSON.stringify({ load: nat.load, reduce: nat.reduce, chartify: nat.chartify, bake: nat.bake }));
process.exit(failed ? 1 : 0);
