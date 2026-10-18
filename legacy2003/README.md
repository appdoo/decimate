# legacy2003 — original 2003 diploma-thesis source

This directory holds the **original 2003 source** of the diploma thesis that the
modern [`decimate`](../) library is a port of: an appearance-preserving
mesh-simplification + texture-atlas toolkit, written by Andor Greißl in 2003
(pre-AI, authored by hand).

It is kept here verbatim as the **source of truth** for the modernization: every
algorithm in `decimate` was ported from these files and characterized against
their behavior. It is intentionally **byte-stable** — not modernized, not
reformatted, not bug-fixed. The German comments, `$Author$`/`$Revision$`
SourceSafe keywords, raw `new`/`delete`, and Windows-isms are all original.

One exception to "verbatim": the vendor OpenGL extension headers
`engine/glnvidia.h`, `engine/glati.h` and `engine/wglati.h` (NVIDIA / ATI,
proprietary) were removed for licensing; `engine/ext_gl.h` still includes them,
so the legacy engine no longer builds as-is. The modern library never used them.

## License

MIT — Copyright (C) 2003 Andor Greißl. See [`LICENSE`](LICENSE). This is a
separate notice from the modern library at the repository root; the years and
scope differ.

## Building

This code builds **only under Visual Studio 2003 on Windows** (OpenGL + GLUT,
`windows.h`, `__int64`, etc.). The Visual Studio project files are kept as a
historical build reference:

- `polyreduction.sln` / `polyreduction.vcproj` — VS2003 solution/project
- `polyreduction.dsp` / `polyreduction.dsw` — older MSVC project (DevStudio)

It does **not** build on modern toolchains as-is — that is the whole point of the
`decimate` port at the repository root, which compiles cross-platform
(macOS/Linux/Windows + WebAssembly) on Clang/GCC. Do not try to "fix" this tree;
fix the port instead.

## Layout (the pipeline, in source form)

| Directory | What it contains |
|---|---|
| `math/` | vector/matrix math, half-edge mesh, sparse matrices (`smcr`, `sparsematrix`), bitsets |
| `loader/` | OBJ + Open Inventor (`.iv`) loaders, lexer |
| `reduce/` | QEM edge-collapse simplification (`edgecollapse`, `cquadricerrorsystem`, …) |
| `chartification/` | region-grow chart segmentation, LSCM conformal mapper, atlas packer (`catlaspacker`, `chorizon`) |
| `parametrization/` | parametrization support |
| `split/` | the `Splitter` (curvature feature-curve charting) + `texturesample` (normal bake) |
| `raytracing/` | ray/triangle intersection (used by the bake) |
| `engine/` | OpenGL/GLUT rendering engine (replaced by the modern `ui/` layer in `decimate`) |
| `image/` | image I/O |
| `main.cpp` | the interactive test driver |

## What is NOT tracked here

Build artifacts and binaries are excluded from git (see the repository
`.gitignore`): `glut32.lib`, `Debug/`, `Release/`, IDE caches
(`*.ncb/*.opt/*.suo/*.plg`), SourceSafe bindings (`*.scc/*.vspscc/*.vssscc`), and
the `*.vsd` design diagram. The `.cpp`/`.h` sources and the VS project files are
tracked.

## Encoding note

These files are **ISO-8859-1 / CRLF** (not UTF-8). German umlauts in comments
will look wrong under a naive UTF-8 read; use `iconv -f ISO-8859-1 -t UTF-8`
when inspecting them.
