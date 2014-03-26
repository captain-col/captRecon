#ifndef ECaptRecon_hxx_seen
#define ECaptRecon_hxx_seen

#include <ECore.hxx>

namespace CP {
    /// Generate a base exception for the Analysis library
    EXCEPTION(ECaptRecon,ECore);

    /// An exception for when an invalid hit is found.
    EXCEPTION(EReconInvalidHit, ECaptRecon);
};
#endif
