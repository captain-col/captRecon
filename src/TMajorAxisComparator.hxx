#ifndef TMajorAxisComparator_hxx_seen
#define TMajorAxisComparator_hxx_seen

#include <TReconBase.hxx>
#include <TReconCluster.hxx>

#include <TPrincipal.h>


namespace CP {
    class TMajorAxisComparator;
};

/// A class to take a TReconObjectContainer and compare it's objects by their
/// position along the major axis.
class CP::TMajorAxisComparator {
public:
    /// Create a sorter that will orient the hits along the principal axis 
    explicit TMajorAxisComparator(const CP::TReconObjectContainer& container) {
        fPrincipal = new TPrincipal(3,"");
        for (CP::TReconObjectContainer::const_iterator c = container.begin();
             c != container.end(); ++c) {
            CP::THandle<CP::TReconCluster> cluster = *c;
            if (!cluster) continue;
            double row[3];
            cluster->GetPosition().Vect().GetXYZ(row);
            fPrincipal->AddRow(row);
        }
        fPrincipal->MakePrincipals();
    }
    
    bool operator() (const CP::THandle<CP::TReconBase> lhs,
                     const CP::THandle<CP::TReconBase> rhs) {
        CP::THandle<CP::TReconCluster> lc = lhs;
        CP::THandle<CP::TReconCluster> rc = rhs;
        double lX[3] = {lc->GetPosition().X(),
                        lc->GetPosition().Y(),
                        lc->GetPosition().Z()};
        double lP[3] = {0,0,0};
        fPrincipal->X2P(lX,lP);
        double rX[3] = {rc->GetPosition().X(),
                        rc->GetPosition().Y(),
                        rc->GetPosition().Z()};
        double rP[3] = {0,0,0};
        fPrincipal->X2P(rX,rP);
        return lP[0] < rP[0];
    }
    
private:
    TPrincipal* fPrincipal;
};

#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// End:
