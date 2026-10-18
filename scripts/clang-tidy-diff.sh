#!/usr/bin/env bash
# scripts/clang-tidy-diff.sh
#
# Run clang-tidy as ERRORS on files touched by the current diff vs. a base ref.
#
# Rollout contract (Spec 03-04 / tasks/spec_map.md Notes):
#   - Globally: .clang-tidy checks ride as warnings (WarningsAsErrors: '' means none).
#   - On diff:  this script promotes them to errors for touched/new files only.
#               Legacy code that hasn't been modernized yet is never force-linted.
#
# Usage:
#   # Errors on all files changed since main:
#   ./scripts/clang-tidy-diff.sh
#
#   # Errors on files changed in the last N commits:
#   ./scripts/clang-tidy-diff.sh HEAD~N
#
#   # Errors on staged changes only:
#   ./scripts/clang-tidy-diff.sh --staged
#
# Requirements:
#   - clang-tidy in PATH, or CLANG_TIDY env var pointing to the binary.
#   - build/compile_commands.json must exist (cmake -B build first).
#   - Python 3 in PATH (used by clang-tidy-diff.py if present in LLVM bin).
#
# Exit codes:
#   0  — no errors on touched files
#   1  — one or more error-level findings on touched files
#   2  — usage / environment error

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${REPO_ROOT}/build"
COMPILE_DB="${BUILD_DIR}/compile_commands.json"

# Honour override, fall back to PATH, then common Homebrew LLVM locations.
CLANG_TIDY="${CLANG_TIDY:-}"
if [[ -z "${CLANG_TIDY}" ]]; then
    for candidate in \
        clang-tidy \
        /opt/homebrew/opt/llvm/bin/clang-tidy \
        /usr/local/opt/llvm/bin/clang-tidy; do
        if command -v "${candidate}" &>/dev/null; then
            CLANG_TIDY="${candidate}"
            break
        fi
    done
fi

if [[ -z "${CLANG_TIDY}" ]]; then
    echo "ERROR: clang-tidy not found. Install via 'brew install llvm' or set CLANG_TIDY env var." >&2
    exit 2
fi

if [[ ! -f "${COMPILE_DB}" ]]; then
    echo "ERROR: ${COMPILE_DB} not found. Run 'cmake -B build' first." >&2
    exit 2
fi

# Determine the diff source.
STAGED=false
BASE_REF=""
if [[ "${1:-}" == "--staged" ]]; then
    STAGED=true
elif [[ -n "${1:-}" ]]; then
    BASE_REF="${1}"
else
    # Default: diff against the merge-base with main.
    BASE_REF="$(git -C "${REPO_ROOT}" merge-base HEAD main 2>/dev/null || git -C "${REPO_ROOT}" rev-parse main 2>/dev/null || echo "")"
    if [[ -z "${BASE_REF}" ]]; then
        echo "ERROR: cannot determine base ref — pass one explicitly, e.g. HEAD~1" >&2
        exit 2
    fi
fi

# Collect the unified diff.
if [[ "${STAGED}" == true ]]; then
    DIFF="$(git -C "${REPO_ROOT}" diff --cached -U0)"
else
    DIFF="$(git -C "${REPO_ROOT}" diff "${BASE_REF}" HEAD -U0)"
fi

if [[ -z "${DIFF}" ]]; then
    echo "No changes detected in diff — nothing to lint."
    exit 0
fi

# Prefer clang-tidy-diff.py from the same LLVM install as clang-tidy.
LLVM_BIN="$(dirname "${CLANG_TIDY}")"
TIDY_DIFF_PY="${LLVM_BIN}/clang-tidy-diff.py"

if [[ -f "${TIDY_DIFF_PY}" ]]; then
    # clang-tidy-diff.py runs clang-tidy on the line ranges that appear in the diff,
    # and we promote all findings to errors by redirecting the exit code.
    echo "Using clang-tidy-diff.py from ${LLVM_BIN}"
    echo "${DIFF}" | python3 "${TIDY_DIFF_PY}" \
        -clang-tidy-binary "${CLANG_TIDY}" \
        -p1 \
        -path "${BUILD_DIR}" \
        -config-file "${REPO_ROOT}/.clang-tidy" \
        -warnings-as-errors '*' \
        -j "$(sysctl -n hw.logicalcpu 2>/dev/null || nproc 2>/dev/null || echo 4)"
else
    # Fallback: collect touched .cpp/.h files and run clang-tidy directly with -warnings-as-errors.
    echo "clang-tidy-diff.py not found at ${TIDY_DIFF_PY} — falling back to whole-file mode."
    if [[ "${STAGED}" == true ]]; then
        TOUCHED_FILES="$(git -C "${REPO_ROOT}" diff --cached --name-only --diff-filter=ACM \
            | grep -E '\.(cpp|cxx|cc|h|hpp)$' | sed "s|^|${REPO_ROOT}/|" || true)"
    else
        TOUCHED_FILES="$(git -C "${REPO_ROOT}" diff "${BASE_REF}" HEAD --name-only --diff-filter=ACM \
            | grep -E '\.(cpp|cxx|cc|h|hpp)$' | sed "s|^|${REPO_ROOT}/|" || true)"
    fi

    if [[ -z "${TOUCHED_FILES}" ]]; then
        echo "No C++ source files in the diff — nothing to lint."
        exit 0
    fi

    echo "Linting (errors) the following touched files:"
    echo "${TOUCHED_FILES}"

    # Run clang-tidy with -warnings-as-errors='*' so every check becomes an error on these files.
    # shellcheck disable=SC2086
    "${CLANG_TIDY}" \
        -p "${BUILD_DIR}" \
        --config-file="${REPO_ROOT}/.clang-tidy" \
        --warnings-as-errors='*' \
        ${TOUCHED_FILES}
fi
