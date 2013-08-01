#ifndef ostreamTVector3_hxx_seen
#define ostreamTVector3_hxx_seen
#include <ostream>
#if __INCLUDE_LEVEL__ > 1
#error ostreamTVector3.hxx included from an include file.
#else
#warning using ostreamTVector3.hxx in source.
#endif
inline std::ostream& operator<<(std::ostream& s, const TVector3& v) {
    return s << '(' << v.X() << ',' << v.Y() << ',' << v.Z() << ')';
}
#endif
