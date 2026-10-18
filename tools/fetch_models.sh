#!/usr/bin/env bash
#
# fetch_models.sh — download the test meshes decimate's tests/demos use.
#
# These meshes are NOT bundled in this repository: they carry their own,
# non-MIT terms (the Stanford models are redistributable but NON-COMMERCIAL;
# the Cyberware horse has no explicit license grant). We do not redistribute
# them — this script pulls each from its canonical archive into the
# (gitignored) model/ directory, so you obtain them under their own terms.
# Every Georgia Tech download is sha256-verified; a mismatch deletes the file
# and fails.
#
# Default:  bunny.iv + horse.iv  (Georgia Tech Large Geometric Models Archive)
# --bench:  also Armadillo.obj + dragon.obj, the full-resolution Stanford
#           originals for tests/reduce_bench.cpp, converted from PLY by
#           tools/ply2obj.py (needs python3; vertex/face counts are asserted).
#
# NOT fetched: cow.obj. The only public mirror lacks the `vn` lines the cow
# characterization tests pin, so cow stays a local-only fixture and its tests
# skip when it is absent. The only mesh committed to the repo is
# model/spot.obj (CC0, Keenan Crane).
#
# Usage:  tools/fetch_models.sh [--bench]      (MODEL_DIR=... overrides model/)
# See:    model/README.md  (per-model license + provenance)

set -euo pipefail

GT="https://sites.cc.gatech.edu/data_files/large_models"
STANFORD="http://graphics.stanford.edu/pub/3Dscanrep"

BUNNY_SHA=f655154f643ef9bd38cc7fc307366fca339f83f80be38f188b267a1b472b91bb
HORSE_SHA=44caaf54477b2dffdfe4de3f5d32448844f7157ad6687abc4d3b9993cc83745f

TOOLS="$(cd "$(dirname "$0")" && pwd)"
DEST="${MODEL_DIR:-$TOOLS/../model}"

BENCH=0
for arg in "$@"; do
    case "$arg" in
        --bench) BENCH=1 ;;
        *) echo "usage: $0 [--bench]" >&2; exit 2 ;;
    esac
done

mkdir -p "$DEST"
DEST="$(cd "$DEST" && pwd)"
SCRATCH="$(mktemp -d)"
trap 'rm -rf "$SCRATCH"' EXIT

# sha256 of a file, portable (macOS shasum / Linux sha256sum)
sha256() {
    if command -v shasum >/dev/null 2>&1; then
        shasum -a 256 "$1" | cut -d' ' -f1
    else
        sha256sum "$1" | cut -d' ' -f1
    fi
}

# verify FILE EXPECTED_SHA — on mismatch delete FILE and fail loudly
verify() {
    local got
    got="$(sha256 "$1")"
    if [[ "$got" != "$2" ]]; then
        rm -f "$1"
        echo "ERROR: sha256 mismatch for $(basename "$1")" >&2
        echo "  expected $2" >&2
        echo "  got      $got" >&2
        echo "  (file deleted)" >&2
        exit 1
    fi
    echo "  ok     $(basename "$1") sha256 ${2:0:16}..."
}

# fetch_gt NAME SHA — Georgia Tech NAME.iv.gz -> gunzip -> DEST/NAME.iv
fetch_gt() {
    local out="$DEST/$1.iv"
    if [[ -f "$out" ]]; then
        echo "  have   $1.iv"
        return
    fi
    echo "  get    $1.iv"
    curl -fsSL "$GT/$1.iv.gz" -o "$SCRATCH/$1.iv.gz"
    gunzip -c "$SCRATCH/$1.iv.gz" > "$SCRATCH/$1.iv"
    verify "$SCRATCH/$1.iv" "$2"
    mv "$SCRATCH/$1.iv" "$out"
}

echo "Fetching test meshes into: $DEST"
echo "(spot.obj is already in the repo — CC0)"

fetch_gt bunny "$BUNNY_SHA"
fetch_gt horse "$HORSE_SHA"

if [[ "$BENCH" -eq 1 ]]; then
    command -v python3 >/dev/null 2>&1 || { echo "ERROR: --bench needs python3" >&2; exit 1; }

    if [[ -f "$DEST/Armadillo.obj" ]]; then
        echo "  have   Armadillo.obj"
    else
        echo "  get    Armadillo.obj (Stanford PLY -> OBJ)"
        curl -fsSL "$STANFORD/armadillo/Armadillo.ply.gz" -o "$SCRATCH/Armadillo.ply.gz"
        gunzip -c "$SCRATCH/Armadillo.ply.gz" > "$SCRATCH/Armadillo.ply"
        python3 -I "$TOOLS/ply2obj.py" "$SCRATCH/Armadillo.ply" "$SCRATCH/Armadillo.obj" \
            --verts 172974 --faces 345944
        mv "$SCRATCH/Armadillo.obj" "$DEST/Armadillo.obj"
    fi

    if [[ -f "$DEST/dragon.obj" ]]; then
        echo "  have   dragon.obj"
    else
        echo "  get    dragon.obj (Stanford PLY -> OBJ)"
        curl -fsSL "$STANFORD/dragon/dragon_recon.tar.gz" -o "$SCRATCH/dragon_recon.tar.gz"
        tar -xzf "$SCRATCH/dragon_recon.tar.gz" -C "$SCRATCH" dragon_recon/dragon_vrip.ply
        python3 -I "$TOOLS/ply2obj.py" "$SCRATCH/dragon_recon/dragon_vrip.ply" "$SCRATCH/dragon.obj" \
            --verts 437645 --faces 871414
        mv "$SCRATCH/dragon.obj" "$DEST/dragon.obj"
    fi
fi

cat <<'NOTE'

Done. NOTE on terms:
  - Stanford models (bunny, armadillo, dragon) are NON-COMMERCIAL use only and
    ask you to credit the Stanford Computer Graphics Laboratory. The dragon
    additionally asks: don't animate or morph it, no Boolean operators.
  - horse (Cyberware, Inc.) has no explicit redistribution license — test/
    research use.
  - cow.obj is not fetched; it is a local-only fixture and the cow-dependent
    tests skip when it is absent.
  - For commercial use of the Stanford models, email the Stanford lab.
Full provenance + terms: model/README.md
NOTE
