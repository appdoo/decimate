#!/usr/bin/env python3
"""ply2obj.py - minimal PLY -> OBJ converter (python3 stdlib only).

Usage: ply2obj.py in.ply out.obj [--verts N --faces M]

Handles ascii, binary_little_endian and binary_big_endian. Reads the vertex
element's x/y/z (any extra vertex properties are skipped by their declared
type) and the face element's vertex_indices list. Faces are emitted as-is
(triangles stay triangles). Output is `v` lines then `f` lines, 1-based.
With --verts/--faces the converted counts are asserted (used by
fetch_models.sh to catch truncated or wrong downloads).
"""
import struct
import sys

TYPES = {
    "char": "b", "int8": "b", "uchar": "B", "uint8": "B",
    "short": "h", "int16": "h", "ushort": "H", "uint16": "H",
    "int": "i", "int32": "i", "uint": "I", "uint32": "I",
    "float": "f", "float32": "f", "double": "d", "float64": "d",
}


def parse_header(f):
    if f.readline().strip() != b"ply":
        sys.exit("not a PLY file")
    fmt, elements = None, []  # elements: [name, count, [(prop, kind, types)]]
    while True:
        line = f.readline()
        if not line:
            sys.exit("unterminated PLY header")
        tok = line.decode("ascii", "replace").split()
        if not tok or tok[0] in ("comment", "obj_info"):
            continue
        if tok[0] == "format":
            fmt = tok[1]
        elif tok[0] == "element":
            elements.append([tok[1], int(tok[2]), []])
        elif tok[0] == "property":
            if tok[1] == "list":
                elements[-1][2].append((tok[4], "list", (tok[2], tok[3])))
            else:
                elements[-1][2].append((tok[2], "scalar", (tok[1],)))
        elif tok[0] == "end_header":
            return fmt, elements


def read_ascii(f, elements):
    data = {}
    for name, count, props in elements:
        rows = []
        for _ in range(count):
            tok = f.readline().decode("ascii").split()
            pos, row = 0, {}
            for pname, kind, _types in props:
                if kind == "scalar":
                    row[pname] = tok[pos]
                    pos += 1
                else:
                    n = int(tok[pos])
                    row[pname] = [int(t) for t in tok[pos + 1:pos + 1 + n]]
                    pos += 1 + n
            rows.append(row)
        data[name] = rows
    return data


def read_binary(f, elements, endian):
    data = {}
    for name, count, props in elements:
        rows = []
        for _ in range(count):
            row = {}
            for pname, kind, types in props:
                if kind == "scalar":
                    c = TYPES[types[0]]
                    row[pname] = struct.unpack(endian + c, f.read(struct.calcsize(c)))[0]
                else:
                    cc, ic = TYPES[types[0]], TYPES[types[1]]
                    n = struct.unpack(endian + cc, f.read(struct.calcsize(cc)))[0]
                    row[pname] = list(struct.unpack(endian + "%d%s" % (n, ic),
                                                    f.read(n * struct.calcsize(ic))))
            rows.append(row)
        data[name] = rows
    return data


def main():
    args = sys.argv[1:]
    want_v = want_f = None
    if "--verts" in args:
        i = args.index("--verts")
        want_v = int(args[i + 1])
        del args[i:i + 2]
    if "--faces" in args:
        i = args.index("--faces")
        want_f = int(args[i + 1])
        del args[i:i + 2]
    if len(args) != 2:
        sys.exit(__doc__)

    with open(args[0], "rb") as f:
        fmt, elements = parse_header(f)
        if fmt == "ascii":
            data = read_ascii(f, elements)
        elif fmt == "binary_little_endian":
            data = read_binary(f, elements, "<")
        elif fmt == "binary_big_endian":
            data = read_binary(f, elements, ">")
        else:
            sys.exit("unsupported PLY format: %s" % fmt)

    verts, faces = data.get("vertex", []), data.get("face", [])
    with open(args[1], "w") as out:
        for v in verts:
            out.write("v %s %s %s\n" % (v["x"], v["y"], v["z"]))
        for fc in faces:
            idx = fc.get("vertex_indices", fc.get("vertex_index"))
            out.write("f " + " ".join(str(i + 1) for i in idx) + "\n")

    print("%s: %d vertices, %d faces" % (args[1], len(verts), len(faces)))
    if want_v is not None:
        assert len(verts) == want_v, "vertex count %d != expected %d" % (len(verts), want_v)
    if want_f is not None:
        assert len(faces) == want_f, "face count %d != expected %d" % (len(faces), want_f)


if __name__ == "__main__":
    main()
