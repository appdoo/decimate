// This file is part of the decimate library; see decimate.h for version/license details.
#include <decimate.h>

// The umbrella target's one translation unit: defines the version symbol (and
// keeps `decimate` a real, linkable archive).

namespace decimate {

/// Packs the DECIMATE_VERSION_* macros of decimate.h into one comparable int.
int version() noexcept {
    return DECIMATE_VERSION_MAJOR * 100 + DECIMATE_VERSION_MINOR;
}

} // namespace decimate
