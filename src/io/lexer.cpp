// This file is part of the decimate library; see decimate.h for version/license details.
/// @file lexer.cpp
/// @brief Character-by-character scanner for ASCII Open Inventor (.iv) files.
///
/// Ported from legacy2003/loader/lexer.cpp (Spec 06, Beat 2a).
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

#include "parse_number.h"

#include <cctype>
#include <sstream>
#include <string>

namespace decimate {

/// Portable case-insensitive compare; replaces the MSVC-only stricmp. The lexer already
/// uppercases each character, so this is mostly defensive.
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


Lexer::Lexer(const std::string &fileName)
    : in_(fileName) {
}


/// Single-character punctuation tokens (state 0). Returns the token value, or
/// -1 if `c` is not a standalone punctuation character.
static int singleCharToken(char c) {
    switch (c) {
        case '<':
            return static_cast<int>(Token::ANGLEOPEN);
        case '>':
            return static_cast<int>(Token::ANGLECLOSE);
        case '(':
            return static_cast<int>(Token::PARENOPEN);
        case ')':
            return static_cast<int>(Token::PARENCLOSE);
        case '[':
            return static_cast<int>(Token::SQUAREOPEN);
        case ']':
            return static_cast<int>(Token::SQUARECLOSE);
        case '{':
            return static_cast<int>(Token::BROPEN);
        case '}':
            return static_cast<int>(Token::BRCLOSE);
        case ',':
            return static_cast<int>(Token::COMMA);
        case '/':
            return static_cast<int>(Token::SLASH);
        default:
            return -1;
    }
}


/// Reads one character, uppercased (the scanner is case-insensitive); returns '\0' at EOF
/// (a NUL byte in the file is indistinguishable from EOF).
char Lexer::nextUpper(char &c) {
    if (!in_.get(c)) {
        return '\0';
    }
    c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return c;
}


/// Token for a finished number lexem.  A lone '-' (no digit ever followed) is not a
/// number: it is UNDEFINED, so match(INT) reports a parse error instead of the 0 that
/// atoi used to read from it.
static int numberToken(const std::string &lexem, int status) {
    if (lexem == "-") {
        return static_cast<int>(Token::UNDEFINED);
    }
    return static_cast<int>(status == 1 ? Token::INT : Token::FLOAT);
}


/// Finishes scanning a number after its first character (digit or '-') is already
/// in lexem_.  Returns INT or FLOAT; puts back the first non-number character.
int Lexer::scanNumber() {
    int  status = 1; // 1=int body, 2=frac, 5=exp-sign, 6=exp-digits
    char c = 0;

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
            return numberToken(lexem_, status);
        }
        lexem_ += c;
    }
    return numberToken(lexem_, status);
}


/// Finishes scanning an alphanumeric keyword/identifier (first letter already in
/// lexem_).  Returns the registered keyword token or Token::STRING.
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


/// Finishes scanning a double-quoted string (opening quote already consumed).
///
/// DELIBERATE, SAFE deviation from legacy (flagged, not silent): the 2009 state
/// machine double-appended each char and left the closing quote in the lexem
/// (a buffer quirk).  This produces the clean string content instead.  The only
/// observable difference is the text of a STRING token's lexem, which is used
/// ONLY in error messages — never for data extraction — and no quoted strings
/// occur in the data sections of the .iv assets.  Numeric output is unchanged
/// (horse.iv golden confirms).
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


/// Reads and uppercases one character at a time, classifies it, and delegates to a per-kind
/// sub-scanner.  Restructures the legacy switch(status) DFA: each old state group becomes one
/// scanner (number = states 1/2/5/6, keyword = 3, string = 4), with the same token output.
/// Helpers append to lexem_ and put back the delimiter as the legacy state machine did.
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

        // Whitespace separates tokens; only newlines are counted, for error messages.
        if (c == '\n') {
            ++line;
        } else if (c != ' ' && c != '\t') {
            return static_cast<int>(Token::UNDEFINED);
        }
    }
    return EOF;
}


void Lexer::match(int expected) {
    if (lex() != expected) {
        std::ostringstream oss;
        oss << "Parse error reading file. Line: " << line
            << ", unexpected token '" << lexem_ << "'";
        throw ParseError(oss.str());
    }
}


int Lexer::getInt() const {
    return io_detail::parseInt(lexem_);
}


float Lexer::getFloat() const {
    return io_detail::parseFloat(lexem_);
}


/// Maps @p lexem to its registered keyword token, or Token::STRING if it is not a keyword.
int Lexer::checkLexem(const std::string &lexem) const {
    // Linear scan (13 Inventor keywords — O(13), negligible).
    for (const Keyword &kw : keywords_) {
        if (caseInsensitiveEqual(kw.text, lexem)) {
            return kw.token;
        }
    }
    return static_cast<int>(Token::STRING);
}

} // namespace decimate
