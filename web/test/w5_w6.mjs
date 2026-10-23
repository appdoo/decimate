// W5 + W6 (spec 34), headless Chrome. usage: node w5_w6.mjs <web dir> <out dir>
//   env PUPPETEER_CORE=<path to a puppeteer-core install>, CHROME=<chrome binary>
// Not part of ctest (needs a browser); run by hand, screenshots land in <out dir>.
import { createServer } from "node:http";
import { readFile, mkdir } from "node:fs/promises";
import { join, extname } from "node:path";
import { pathToFileURL } from "node:url";

const [, , webDir, outDir] = process.argv;
const puppeteer = (await import(pathToFileURL(join(process.env.PUPPETEER_CORE, "lib/esm/puppeteer/puppeteer-core.js")))).default;
await mkdir(outDir, { recursive: true });
const types = { ".html": "text/html", ".js": "text/javascript", ".wasm": "application/wasm", ".css": "text/css", ".obj": "text/plain" };
const server = createServer(async (req, res) => {
    try {
        const p = join(webDir, decodeURIComponent(new URL(req.url, "http://x").pathname).replace(/\/$/, "/index.html"));
        const body = await readFile(p);
        res.writeHead(200, { "content-type": types[extname(p)] ?? "application/octet-stream" }).end(body);
    } catch { res.writeHead(404).end(); }
}).listen(0);
const url = `http://127.0.0.1:${server.address().port}/`;

let failed = 0;
const check = (ok, msg) => { console.log((ok ? "PASS " : "FAIL ") + msg); if (!ok) failed++; };
const browser = await puppeteer.launch({
    executablePath: process.env.CHROME ?? "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
    headless: "new", args: ["--use-angle=swiftshader", "--enable-unsafe-swiftshader", "--ignore-gpu-blocklist"],
    defaultViewport: { width: 1440, height: 860 },
});
const open = async () => {
    const page = await browser.newPage();
    page.errors = [];
    page.on("pageerror", (e) => page.errors.push(String(e)));
    page.on("console", (m) => m.type() === "error" && page.errors.push(m.text()));
    await page.goto(url);
    return page;
};
const waitState = (page, s) => page.waitForFunction((s) => document.body.dataset.state === s, { timeout: 90000 }, s);
const shot = async (page, name) => { const f = join(outDir, name); await page.screenshot({ path: f }); console.log("screenshot", f); };

// ---- W5: drop spot.obj, process, export GLB, view it again
let page = await open();
await page.evaluate(async () => {
    const bytes = await (await fetch("spot.obj")).arrayBuffer();
    const dt = new DataTransfer();
    dt.items.add(new File([bytes], "spot.obj"));
    window.dispatchEvent(new DragEvent("dragenter", { dataTransfer: dt, bubbles: true }));
    window.dispatchEvent(new DragEvent("drop", { dataTransfer: dt, bubbles: true, cancelable: true }));
});
await waitState(page, "done");
check(true, "W5 drop spot.obj -> processed (object space)");
await shot(page, "w5_object_space.png");
await page.select("#space", "tangent");
await page.click("#run");
await waitState(page, "busy").catch(() => {});
await waitState(page, "done");
await shot(page, "w5_tangent_space.png");
check(await page.$eval("#dlGlb", (b) => !b.disabled), "W5 GLB button enabled for tangent space");
const stats = await page.$eval("#stats", (e) => e.textContent);
console.log("stats:", stats);
await page.click("#tWire");
await shot(page, "w5_wireframe.png");
await page.click("#tWire");
await page.click("#tAtlas");
await shot(page, "w5_atlas.png");
await page.click("#tAtlas");
const glb = await page.evaluate(async () => Array.from(await window.__app.glb()));
check(glb.length > 1000 && String.fromCharCode(...glb.slice(0, 4)) === "glTF", `W5 exported GLB ${glb.length} bytes`);
await page.evaluate((bytes) => window.__app.loadBytes("roundtrip.glb", new Uint8Array(bytes)), glb);
await waitState(page, "glb");
await new Promise((r) => setTimeout(r, 800));
await shot(page, "w5_glb_roundtrip.png");
check(page.errors.length === 0, `W5 no page errors ${JSON.stringify(page.errors)}`);
await page.close();

// ---- W6: malformed file -> C++ error text, worker survives a second file
page = await open();
await page.evaluate(() => window.__app.loadBytes("broken.obj", new TextEncoder().encode("v 1 2 3\nf 1 2 3\n")));
await waitState(page, "error");
const msg = await page.$eval("#toast", (e) => e.hidden ? "" : e.textContent);
console.log("toast:", msg);
check(msg.length > 0 && !/^undefined|\[object/.test(msg), `W6 error toast shows the C++ message: "${msg}"`);
await shot(page, "w6_error.png");
await page.evaluate(async () => window.__app.loadBytes("spot.obj", new Uint8Array(await (await fetch("spot.obj")).arrayBuffer())));
await waitState(page, "done");
check(true, "W6 worker survived: spot.obj processed after the error");
await page.close();

// dark scheme
page = await open();
await page.emulateMediaFeatures([{ name: "prefers-color-scheme", value: "dark" }]);
await page.evaluate(async () => window.__app.loadBytes("spot.obj", new Uint8Array(await (await fetch("spot.obj")).arrayBuffer())));
await waitState(page, "done");
await page.click("#tAtlas");
await shot(page, "w5_dark_atlas.png");
await page.close();

// mobile width
page = await open();
await page.setViewport({ width: 390, height: 844 });
await page.evaluate(async () => window.__app.loadBytes("spot.obj", new Uint8Array(await (await fetch("spot.obj")).arrayBuffer())));
await waitState(page, "done");
check(await page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth), "mobile: no horizontal scroll");
await page.screenshot({ path: join(outDir, "w5_mobile.png"), fullPage: true });
await page.close();
await browser.close();
server.close();
process.exit(failed ? 1 : 0);
