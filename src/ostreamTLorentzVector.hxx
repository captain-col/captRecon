#ifndef ostreamTLorentzVector_hxx_seen
#define ostreamTLorentzVector_hxx_seen
#include <ostream>

#if __INCLUDE_LEVEL__ > 1
#error ostreamTLorentzVector.hxx included from an include file.
#else
#warning using ostreamTLorentzVector.hxx in source.
#endif
inline std::ostream& operator<<(std::ostream& s, const TLorentzVector& v) {
    return s << '('<<v.X()<<','<<v.Y()<<','<<v.Z()<<','<<v.T()<<')';
}
#endif
