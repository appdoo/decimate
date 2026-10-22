// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

/// @file parse_number.h
/// @brief Locale-independent, stateless number parsing shared by the OBJ loader and the
///        Inventor lexer (private to src/io/, not installed).
///
/// Replaces atof/atoi.  Those read the process-global LC_NUMERIC locale, so a "de_DE"
/// process turns "1.5" into 1.  These helpers take no locale and keep no state.  They
/// keep atof/atoi's leniency, which the 2003 readers relied on: an optional leading
/// '+', ".5", "1e5", trailing garbage ("1.0abc" -> 1.0, "251//255" -> 251) and 0 when
/// no number is there at all.

#include <charconv>
#include <fast_float/fast_float.h>
#include <string_view>
#include <system_error>

namespace decimate::io_detail {

/// Drops one leading '+', which atof/atoi accept and the converters below do not.
inline std::string_view skipPlus(std::string_view token) {
    if (!token.empty() && token.front() == '+') {
        token.remove_prefix(1);
    }
    return token;
}


/// atof replacement: the longest numeric prefix of @p token, 0 when there is none.
/// Parsed as double and then narrowed, exactly as `float(atof(...))` did, so the float
/// is bit-identical.  Out of range gives +-inf or 0 like strtod.
inline float parseFloat(std::string_view token) {
    token = skipPlus(token);
    double     value = 0.0;
    const auto result = fast_float::from_chars(token.data(), token.data() + token.size(), value,
                                               fast_float::chars_format::general);
    if (result.ec == std::errc::invalid_argument) {
        return 0.0F;
    }
    return static_cast<float>(value);
}


/// atoi replacement: the leading integer of @p token, 0 when there is none or it does
/// not fit an int (atoi is undefined there).  Parsing stops at the first non-digit, so
/// "251//255" yields 251 (FROZEN BUG 2 in obj_loader.cpp).
inline int parseInt(std::string_view token) {
    token = skipPlus(token);
    int        value = 0;
    const auto result = std::from_chars(token.data(), token.data() + token.size(), value);
    return result.ec == std::errc {} ? value : 0;
}

} // namespace decimate::io_detail
