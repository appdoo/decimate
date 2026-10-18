// This file is part of the decimate library; see decimate.h for version/license details.
/// @file inventor_loader.cpp
/// @brief Open Inventor ASCII (.iv) mesh loader — behavior-preserving port of
///        legacy InventorLdr::load.
///
/// Ported from legacy/loader/InventorLdr.cpp (Spec 06, Beat 2a).
///
/// Representation changes (no algorithmic or numeric effect):
///   - Output is MeshData (owned std::vector) instead of raw new[]/out-params.
///   - The vector → memcpy → new[] double-allocation in readCoordinate3 and
///     readIndexedFaceSet is dropped; we accumulate directly into result.
///   - Lexer is a stack value (not heap-allocated via new ifstream*).
///
/// Algorithm is UNCHANGED from legacy:
///   - load → readSeparator → readCoordinate3 / readIndexedFaceSet
///     recursive-descent over the Lexer token stream.
///   - Coordinate3: reads comma-separated XYZ float triples.
///   - IndexedFaceSet: reads exactly four ints per face (three indices + the
///     -1 terminator consumed and discarded at InventorLdr.cpp:210-211).
///   - Material / MaterialBinding: silently skipped.
///   - Triangle-only: hardcodes four-token reads per face (InventorLdr behavior).

#include <decimate/types/mesh_data.h>
#include <decimate/io/lexer.h>

#include <array>
#include <stdexcept>
#include <string>
#include <vector>

namespace decimate {

// ---------------------------------------------------------------------------
// IV keyword tokens
// ---------------------------------------------------------------------------
// Token values mirror the legacy InventorLdr enum: the first keyword is
// Token::UserToken (15) and each subsequent keyword adds 1.

namespace iv {

constexpr int kBase = static_cast<int>(Token::UserToken);

enum Tok : int {
    kInventor = kBase, kSeparator, kCoordinate3, kPoint, kMaterial, kDiffuseColor,
    kMaterialBinding, kValue, kPerVertexIndexed, kIndexedFaceSet, kCoordIndex,
    kMaterialIndex, kAscii,
};

// Keyword table (uppercased; the lexer matches case-insensitively). The order
// is irrelevant to correctness — each entry pairs its spelling with its token.
constexpr std::array<Lexer::Keyword, 13> kKeywords{{
    {"#INVENTOR", kInventor},          {"SEPARATOR", kSeparator},
    {"COORDINATE3", kCoordinate3},     {"POINT", kPoint},
    {"MATERIAL", kMaterial},           {"DIFFUSECOLOR", kDiffuseColor},
    {"MATERIALBINDING", kMaterialBinding}, {"VALUE", kValue},
    {"PER_VERTEX_INDEXED", kPerVertexIndexed}, {"INDEXEDFACESET", kIndexedFaceSet},
    {"COORDINDEX", kCoordIndex},       {"MATERIALINDEX", kMaterialIndex},
    {"ASCII", kAscii},
}};

}  // namespace iv

// ---------------------------------------------------------------------------
// readCoordinate3 — port of InventorLdr::readCoordinate3 (InventorLdr.cpp:108)
// ---------------------------------------------------------------------------
// Reads: { POINT [ <float triples separated by commas> ] }
// Populates coords with one Vec3f per vertex.

// Reads one "x y z" triple. The first component is already lexed (FLOAT or INT,
// read via getFloat()); the next two are matched as FLOAT. (Legacy InventorLdr
// only matched FLOAT; an INT-first component reads identically via getFloat().)
static decimate::math::Vec3f readVertexTriple(Lexer& lexer) {
    decimate::math::Vec3f v;
    v.x = lexer.getFloat();
    lexer.match(static_cast<int>(Token::FLOAT));
    v.y = lexer.getFloat();
    lexer.match(static_cast<int>(Token::FLOAT));
    v.z = lexer.getFloat();
    return v;
}

static void readCoordinate3(Lexer& lexer, std::vector<decimate::math::Vec3f>& coords) {
    lexer.match(static_cast<int>(Token::BROPEN));
    lexer.match(iv::kPoint);
    lexer.match(static_cast<int>(Token::SQUAREOPEN));

    int token = 0;
    while ((token = lexer.lex()) != EOF) {
        switch (token) {
            // FLOAT or INT both begin a vertex triple (see readVertexTriple).
            case static_cast<int>(Token::FLOAT):
            case static_cast<int>(Token::INT):
                coords.push_back(readVertexTriple(lexer));
                break;

            case static_cast<int>(Token::COMMA):
                break;  // separator between triples, skip

            case static_cast<int>(Token::SQUARECLOSE):
                lexer.match(static_cast<int>(Token::BRCLOSE));
                return;

            default:
                break;
        }
    }
    throw ParseError("Unexpected end of file in Coordinate3 block");
}

// ---------------------------------------------------------------------------
// readMaterial — port of InventorLdr::readMaterial (InventorLdr.cpp:149)
// ---------------------------------------------------------------------------
// Silently consumes everything until the matching closing brace.

static void readMaterial(Lexer& lexer) {
    lexer.match(static_cast<int>(Token::BROPEN));
    int token = 0;
    while ((token = lexer.lex()) != EOF) {
        if (token == static_cast<int>(Token::BRCLOSE)) {
            return;
        }
    }
    throw ParseError("Unexpected end of file in Material block");
}

// ---------------------------------------------------------------------------
// readMaterialBinding — port of InventorLdr::readMaterialBinding (line 170)
// ---------------------------------------------------------------------------

static void readMaterialBinding(Lexer& lexer) {
    lexer.match(static_cast<int>(Token::BROPEN));
    int token = 0;
    while ((token = lexer.lex()) != EOF) {
        if (token == static_cast<int>(Token::BRCLOSE)) {
            return;
        }
    }
    throw ParseError("Unexpected end of file in MaterialBinding block");
}

// ---------------------------------------------------------------------------
// readIndexedFaceSet — port of InventorLdr::readIndexedFaceSet (line 191)
// ---------------------------------------------------------------------------
// Reads: { COORDINDEX [ <v0, v1, v2, -1,>... ] }
// The -1 terminator is the fourth INT in each group — it is consumed and
// discarded (InventorLdr.cpp:210-211).  Triangle-only: exactly 3 indices + -1
// per face.

static void readIndexedFaceSet(Lexer& lexer, std::vector<int>& faces) {
    lexer.match(static_cast<int>(Token::BROPEN));
    lexer.match(iv::kCoordIndex);
    lexer.match(static_cast<int>(Token::SQUAREOPEN));

    int token = 0;
    while ((token = lexer.lex()) != EOF) {
        switch (token) {
            case static_cast<int>(Token::INT): {
                // First vertex index of a triangle.
                faces.push_back(lexer.getInt());
                lexer.match(static_cast<int>(Token::COMMA));

                lexer.match(static_cast<int>(Token::INT));
                faces.push_back(lexer.getInt());
                lexer.match(static_cast<int>(Token::COMMA));

                lexer.match(static_cast<int>(Token::INT));
                faces.push_back(lexer.getInt());
                lexer.match(static_cast<int>(Token::COMMA));

                // Consume the -1 terminator (InventorLdr.cpp:210-211).
                lexer.match(static_cast<int>(Token::INT));
                // -1 is read but not pushed; discard.
                lexer.match(static_cast<int>(Token::COMMA));
                break;
            }

            case static_cast<int>(Token::COMMA):
                break;  // extra comma between faces, skip

            case static_cast<int>(Token::SQUARECLOSE):
                return;

            default:
                break;
        }
    }
    throw ParseError("Unexpected end of file in IndexedFaceSet block");
}

// ---------------------------------------------------------------------------
// readSeparator — port of InventorLdr::readSeparator (InventorLdr.cpp:72)
// ---------------------------------------------------------------------------

static void readSeparator(Lexer& lexer, std::vector<decimate::math::Vec3f>& coords, std::vector<int>& faces) {
    lexer.match(static_cast<int>(Token::BROPEN));

    int token = 0;
    while ((token = lexer.lex()) != EOF) {
        switch (token) {
            case iv::kCoordinate3:
                readCoordinate3(lexer, coords);
                break;

            case iv::kMaterial:
                readMaterial(lexer);
                break;

            case iv::kMaterialBinding:
                readMaterialBinding(lexer);
                break;

            case iv::kIndexedFaceSet:
                readIndexedFaceSet(lexer, faces);
                break;

            case static_cast<int>(Token::BRCLOSE):
                return;

            default:
                break;  // ignore unrecognised tokens inside a Separator
        }
    }
    throw ParseError("Unexpected end of file in Separator block");
}

// ---------------------------------------------------------------------------
// loadInventor — port of InventorLdr::load (InventorLdr.cpp:48)
// ---------------------------------------------------------------------------

MeshData loadInventor(std::string_view path) {
    Lexer lexer{std::string(path)};
    if (!lexer) {
        throw std::runtime_error(
            std::string("loadInventor: unable to open file '") + std::string(path) + "'");
    }

    // Register the 13 Inventor keyword tokens.
    lexer.setKeywords(iv::kKeywords);

    MeshData result;

    int token = 0;
    while ((token = lexer.lex()) != EOF) {
        if (token == iv::kSeparator) {
            readSeparator(lexer, result.positions, result.indices);
        }
    }

    // Trust boundary: reject indices outside [0, positions.size()).
    for (const int idx : result.indices) {
        if (idx < 0 || static_cast<std::size_t>(idx) >= result.positions.size()) {
            throw std::runtime_error("loadInventor: coordIndex " + std::to_string(idx) +
                                     " out of range (vertex count " +
                                     std::to_string(result.positions.size()) + ")");
        }
    }

    return result;
}

}  // namespace decimate
