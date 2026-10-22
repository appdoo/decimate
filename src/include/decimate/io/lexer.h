// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

/// @file lexer.h
/// @brief Character-by-character scanner for ASCII Open Inventor (.iv) files.
///
/// Ported from legacy2003/loader/lexer.h + legacy2003/loader/lexer.cpp.  Portability fixes applied:
///   - Removed #include <windows.h> (was dead code in lexer.cpp:7).
///   - Replaced stricmp with a portable case-insensitive compare.
///   - German error message → English.
///   - No raw owning pointers in public API.
///   - No printf/cout/cerr; parse errors thrown as std::runtime_error.
///
/// This class is consumed ONLY by the .iv parser (InventorLdr).
/// The OBJ parser does not use it.

#include <fstream>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>

namespace decimate {

/// @brief Parse error thrown when the token stream is malformed.
///
/// Carries the file path, line number, and the unexpected token text.
struct ParseError : std::runtime_error {
    /// @param msg   Human-readable description of the parse failure.
    explicit ParseError(const std::string &msg)
        : std::runtime_error(msg) {
    }
};


/// @brief Token types returned by Lexer::lex().
enum class Token : int {
    INT = 0,
    FLOAT = 1,
    STRING = 2,
    BROPEN = 3,       // {
    BRCLOSE = 4,      // }
    PARENOPEN = 5,    // (
    PARENCLOSE = 6,   // )
    ANGLEOPEN = 7,    // <
    ANGLECLOSE = 8,   // >
    SQUAREOPEN = 9,   // [
    SQUARECLOSE = 10, // ]
    COMMA = 11,
    LEXEM = 12,
    UNDEFINED = 13,
    SLASH = 14,
    // User-defined keyword tokens start here.
    // InventorLdr registers its 13 keywords at UserToken + 0..12.
    UserToken = 15,
};


/// @brief Portable keyword-aware scanner for ASCII Open Inventor files.
///
/// Tokenizes a stream one character at a time (preserving legacy behavior —
/// bulk-read optimization is Beat 2b / Finding A).  Supports user-registered
/// keyword tokens returned as integer values >= Token::UserToken.
class Lexer {
public:
    /// @brief Open @p fileName for scanning.
    ///
    /// After construction, check operator bool() before using lex().
    explicit Lexer(const std::string &fileName);

    /// Destructor closes the file (RAII).
    ~Lexer() = default;

    // Non-copyable, movable.
    Lexer(const Lexer &) = delete;
    Lexer &operator=(const Lexer &) = delete;
    Lexer(Lexer &&) = default;
    Lexer &operator=(Lexer &&) = default;

    /// @brief Returns false if the file was not opened successfully.
    explicit operator bool() const {
        return static_cast<bool>(in_);
    }

    /// @brief A keyword and the token value lex() returns when it is matched.
    struct Keyword {
        std::string_view text;  ///< Keyword spelling (matched case-insensitively).
        int              token; ///< Token value returned for this keyword.
    };

    /// @brief Register the keyword table this scanner recognises.
    ///
    /// The span is borrowed (not owned); its storage must outlive this Lexer.
    void setKeywords(std::span<const Keyword> keywords) {
        keywords_ = keywords;
    }

    /// @brief Scan and return the next token.
    ///
    /// @returns The token as an int (cast from Token or a user-registered value),
    ///          or EOF (-1) at end of file.
    int lex();

    /// @brief Consume the next token and assert it equals @p expected.
    ///
    /// @throws ParseError if the next token does not match.
    void match(int expected);

    /// @brief Returns the current lexem string.
    [[nodiscard]] const std::string &getLexem() const {
        return lexem_;
    }

    /// @brief Returns the current lexem parsed as an integer.
    [[nodiscard]] int getInt() const;

    /// @brief Returns the current lexem parsed as a float.
    [[nodiscard]] float getFloat() const;

    /// Current source line (1-based), updated when '\n' is consumed.
    int line { 1 };

private:
    /// Returns the keyword token value for @p lexem (case-insensitive),
    /// or Token::STRING if not found.
    [[nodiscard]] int checkLexem(const std::string &lexem) const;

    /// Read the next character, uppercased, into @p c; returns it, or '\0' at EOF.
    char nextUpper(char &c);

    /// Sub-scanners invoked by lex() after the first character is classified.
    int scanNumber();  ///< INT/FLOAT body (first digit or '-' already in lexem_).
    int scanKeyword(); ///< keyword/identifier (first letter already in lexem_).
    int scanString();  ///< quoted string (opening quote already consumed).

    std::ifstream            in_;
    std::string              lexem_;
    std::span<const Keyword> keywords_;
};

} // namespace decimate
