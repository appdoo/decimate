// This file is part of the decimate library; see decimate.h for version/license details.
/// @file lexer.cpp
/// @brief Character-by-character scanner for ASCII Open Inventor (.iv) files.
///
/// Ported from legacy/loader/lexer.cpp (Spec 06, Beat 2a).
///
/// Portability changes applied (representation, not algorithm):
///   - Removed #include <windows.h> (was dead: no Win32 type used anywhere).
///   - Replaced stricmp(a,b)==0 with portable case-insensitive compare.
///   - lexem_ is std::string instead of StringBuffer (byte-for-byte token
///     output identical — exponential growth vs additive-1024 is representation).
///   - German error string ("Fehler beim Einlesen...") → English.
///   - No printf/cout/cerr; parse failures throw ParseError.
///
/// Algorithm is unchanged from the 2009 legacy source:
///   - Per-character get() state machine (states 0–6).
///   - Case-insensitive uppercasing (c = toupper(c)).
///   - putback on non-consuming state transitions.
///   - Linear scan for keyword lookup.
///
/// NOTE (Beat 2b): The per-character get() path (hot edge 1) is NOT optimized
/// here.  Bulk-read replacement is deferred.

#include <decimate/io/lexer.h>

#include <cctype>
#include <cstdlib>
#include <sstream>
#include <string>

namespace decimate {

// ---------------------------------------------------------------------------
// Portable case-insensitive string compare
// ---------------------------------------------------------------------------
// Replaces legacy stricmp (MSVC extension).  The lexer already uppercases each
// character (c = toupper(c) in the state machine), so this is mostly defensive.
static bool caseInsensitiveEqual(std::string_view a, std::string_view b) {
    if (a.size() != b.size()) {
        return false;
    }
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (std::toupper(static_cast<unsigned char>(a[i])) !=
            std::toupper(static_cast<unsigned char>(b[i]))) {
            return false;
        }
    }
    return true;
}

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------

Lexer::Lexer(const std::string& fileName)
    : in_(fileName) {}

// Single-character punctuation tokens (state 0). Returns the token value, or
// -1 if `c` is not a standalone punctuation character.
static int singleCharToken(char c) {
    switch (c) {
        case '<': return static_cast<int>(Token::ANGLEOPEN);
        case '>': return static_cast<int>(Token::ANGLECLOSE);
        case '(': return static_cast<int>(Token::PARENOPEN);
        case ')': return static_cast<int>(Token::PARENCLOSE);
        case '[': return static_cast<int>(Token::SQUAREOPEN);
        case ']': return static_cast<int>(Token::SQUARECLOSE);
        case '{': return static_cast<int>(Token::BROPEN);
        case '}': return static_cast<int>(Token::BRCLOSE);
        case ',': return static_cast<int>(Token::COMMA);
        case '/': return static_cast<int>(Token::SLASH);
        default:  return -1;
    }
}

// ---------------------------------------------------------------------------
// lex() — tokenizer entry point
// ---------------------------------------------------------------------------
// Reads and uppercases one character at a time (legacy is case-insensitive),
// classifies it, and delegates to a per-kind sub-scanner.  This is a faithful
// restructuring of the legacy switch(status) DFA: each old state group becomes
// one scanner (number = states 1/2/5/6, keyword = 3, string = 4) so token
// output is byte-for-byte identical.  Helpers append to / pop from lexem_ and
// putback() the delimiter exactly as the legacy state machine did.

char Lexer::nextUpper(char& c) {
    if (!in_.get(c)) {
        return '\0';
    }
    c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return c;
}

// Finish scanning a number after its first character (digit or '-') is already
// in lexem_.  Returns INT or FLOAT; putbacks the first non-number character.
int Lexer::scanNumber() {
    int  status = 1;  // 1=int body, 2=frac, 5=exp-sign, 6=exp-digits
    char c      = 0;
    while (nextUpper(c) != '\0') {
        const bool digit = std::isdigit(static_cast<unsigned char>(c)) != 0;
        if (status == 1 && c == '.') {
            status = 2;
        } else if ((status == 1 || status == 2) && c == 'E') {
            status = 5;
        } else if (status == 5 && (c == '+' || c == '-' || digit)) {
            status = 6;
        } else if (!digit) {
            in_.putback(c);
            return static_cast<int>(status == 1 ? Token::INT : Token::FLOAT);
        }
        lexem_ += c;
    }
    return static_cast<int>(status == 1 ? Token::INT : Token::FLOAT);
}

// Finish scanning an alphanumeric keyword/identifier (first letter already in
// lexem_).  Returns the registered keyword token or Token::STRING.
int Lexer::scanKeyword() {
    char c = 0;
    while (nextUpper(c) != '\0') {
        if (std::isalnum(static_cast<unsigned char>(c)) == 0) {
            in_.putback(c);
            return checkLexem(lexem_);
        }
        lexem_ += c;
    }
    return checkLexem(lexem_);
}

// Finish scanning a double-quoted string (opening quote already consumed).
//
// DELIBERATE, SAFE deviation from legacy (flagged, not silent): the 2009 state
// machine double-appended each char and left the closing quote in the lexem
// (a buffer quirk).  This produces the clean string content instead.  The only
// observable difference is the text of a STRING token's lexem, which is used
// ONLY in error messages — never for data extraction — and no quoted strings
// occur in the data sections of the .iv assets.  Numeric output is unchanged
// (horse.iv golden confirms).
int Lexer::scanString() {
    char c = 0;
    while (nextUpper(c) != '\0') {
        if (c == '"') {
            return static_cast<int>(Token::STRING);
        }
        lexem_ += c;
    }
    return static_cast<int>(Token::STRING);
}

int Lexer::lex() {
    lexem_.clear();
    char c = 0;
    while (nextUpper(c) != '\0') {
        if (std::isalpha(static_cast<unsigned char>(c)) != 0) {
            lexem_ += c;
            return scanKeyword();
        }
        if (std::isdigit(static_cast<unsigned char>(c)) != 0 || c == '-') {
            lexem_ += c;
            return scanNumber();
        }
        if (int tok = singleCharToken(c); tok != -1) {
            return tok;
        }
        if (c == '"') {
            return scanString();
        }
        if (c == '\n') {
            ++line;  // newline counts as whitespace; keep scanning
        } else if (c != ' ' && c != '\t') {
            return static_cast<int>(Token::UNDEFINED);
        }
    }
    return EOF;
}

// ---------------------------------------------------------------------------
// match() — consume-and-assert
// ---------------------------------------------------------------------------

void Lexer::match(int expected) {
    if (lex() != expected) {
        std::ostringstream oss;
        oss << "Parse error reading file. Line: " << line
            << ", unexpected token '" << lexem_ << "'";
        throw ParseError(oss.str());
    }
}

// ---------------------------------------------------------------------------
// getInt / getFloat
// ---------------------------------------------------------------------------

int Lexer::getInt() const {
    // NOTE (Beat 2b): atoi is locale-independent for integers but is preserved
    // from legacy.  Replace with std::from_chars in Beat 2b.
    return std::atoi(lexem_.c_str());
}

float Lexer::getFloat() const {
    // NOTE (Beat 2b): atof is locale-dependent (LC_NUMERIC decimal separator).
    // This is the locale hazard from Finding B.  Replace with fast_float /
    // std::from_chars in Beat 2b.
    return static_cast<float>(std::atof(lexem_.c_str()));
}

// ---------------------------------------------------------------------------
// checkLexem() — portable keyword lookup
// ---------------------------------------------------------------------------

int Lexer::checkLexem(const std::string& lexem) const {
    // Linear scan (13 Inventor keywords — O(13), negligible).
    for (const Keyword& kw : keywords_) {
        if (caseInsensitiveEqual(kw.text, lexem)) {
            return kw.token;
        }
    }
    return static_cast<int>(Token::STRING);
}

}  // namespace decimate
