# model/ — test meshes

Only `spot.obj` is in the repository. Everything else is fetched on demand
under its own (non-MIT) terms. Tests whose model is missing are reported as
**Skipped**, not failed.

## Get the models

```sh
tools/fetch_models.sh           # bunny.iv + horse.iv   (sha256-verified)
tools/fetch_models.sh --bench   # + Armadillo.obj + dragon.obj for reduce_bench (needs python3)
```

Files land in this directory (gitignored except `spot.obj` and this README);
`MODEL_DIR=...` fetches elsewhere. A sha256 mismatch deletes the file and
fails. The `--bench` meshes are the full-resolution Stanford originals,
converted from PLY by `tools/ply2obj.py`, which asserts their vertex/face
counts.

## What's here

| File | How you get it | Terms | Source / credit |
|---|---|---|---|
| `spot.obj` | in the repo | **CC0 1.0** | Keenan Crane, [model repository](https://www.cs.cmu.edu/~kmcrane/Projects/ModelRepository/) |
| `bunny.iv` | `fetch_models.sh` | Stanford — non-commercial, credit required | Stanford bunny, via the [Georgia Tech Large Geometric Models Archive](https://sites.cc.gatech.edu/data_files/large_models/bunny.iv.gz) |
| `horse.iv` | `fetch_models.sh` | no explicit license — research/test use | Cyberware, Inc., via the [Georgia Tech archive](https://sites.cc.gatech.edu/data_files/large_models/horse.iv.gz) |
| `Armadillo.obj` | `--bench` | Stanford — non-commercial, credit required | [`Armadillo.ply.gz`](http://graphics.stanford.edu/pub/3Dscanrep/armadillo/Armadillo.ply.gz), 172 974 verts / 345 944 faces |
| `dragon.obj` | `--bench` | Stanford — non-commercial, credit required, no morphing / Booleans | [`dragon_recon.tar.gz`](http://graphics.stanford.edu/pub/3Dscanrep/dragon/dragon_recon.tar.gz) → `dragon_vrip.ply`, 437 645 verts / 871 414 faces |

## Terms

The Stanford models (bunny, Armadillo, dragon) come from the
[Stanford 3D Scanning Repository](https://graphics.stanford.edu/data/3Dscanrep/):

> You are welcome to use the data and models for research purposes. You are also welcome to mirror or redistribute them for free. … as long as credit is given to the Stanford Computer Graphics Laboratory. However, such models or images are not to be used for commercial purposes, nor should they appear in a product for sale … without our permission.

The repository adds, for the dragon, "Don't animate or morph them, don't
apply Boolean operators to them", and: "You can do anything you want to the
Stanford bunny or the armadillo." For commercial use, email the Stanford
Computer Graphics Laboratory.

Why nothing but `spot.obj` is committed: MIT covers this project's code, not
data placed next to it. A bundled bunny would still be non-commercial and
would make the tree mixed-license. So one CC0 mesh is bundled so that
`clone → build → run` always works, and the rest is yours to fetch.

This is a practical reading of the terms, not legal advice.
