#ifndef CreateCluster_hxx_seen
#define CreateCluster_hxx_seen

#include "HitUtilities.hxx"
#include "ECaptRecon.hxx"

#include <TReconCluster.hxx>
#include <TUnitsTable.hxx>
#include <THandle.hxx>
#include <TReconHit.hxx>
#include <THit.hxx>
#include <CaptGeomId.hxx>

namespace CP {

    /// A base exception for the create track template.
    EXCEPTION(ECreateCluster,ECaptRecon);

    /// An exception that there are repeated objects in the track.
    EXCEPTION(EClusterRepeatedObject, ECreateCluster);

    /// An exception that there are repeated objects in the track.
    EXCEPTION(EClusterNonHit, ECreateCluster);

    /// Take iterators to a collection of CP::THandle<THit> objects and
    /// construct a cluster.  A typical example would be to create a cluster
    /// from a THitSelection.
    ///
    /// \code
    /// CP::THitSelection hits;
    /// CP::THandle<CP::TReconCluster> cluster
    ///      = CreateCluster("cluster", hits.begin(), hits.end());
    /// \endcode
    ///
    /// If the optional recalculateUncertainty argument is true, then the
    /// uncertainty is based on the size of the cluster.  This is the
    /// appropriate measure of position uncertainty when the hit represents a
    /// measured charge distribution.
    template<typename hitIterator>
    CP::THandle<CP::TReconCluster> 
    CreateCluster(const char* name, hitIterator begin, hitIterator end,
                  bool recalculateUncertainty = true);
};

//////////////////////////////////////////////////////////////////
// IMPLEMENTATION
//////////////////////////////////////////////////////////////////

template<typename hitIterator>
CP::THandle<CP::TReconCluster> 
CP::CreateCluster(const char* name, hitIterator begin, hitIterator end,
                  bool recalculateUncertainty) {
    
#ifdef DEBUG_CREATE_CLUSTER
    int count1 = 0;
    for (hitIterator i = begin; i!=end; ++i) {
        CP::THandle< CP::THit > h = *i;
        if (!h) {
            CaptError("Non-hit at position " << count1);
            throw CP::EClusterNonHit();
        }
        hitIterator j = i;
        int count2 = 0;
        while ((++j) != end) {
            ++count2;
            if (CP::GetPointer(*i) != CP::GetPointer(*j)) continue;
            CaptError("Invalid cluster: multiple copies of object ("
                      << count1 << " matchs " << count1 + count2 << ")");
            throw CP::EClusterRepeatedObject();
        }
        ++count1;
    }
#endif

    CaptNamedInfo("CreateCluster", name );

    CP::THandle<CP::TReconCluster> cluster(new CP::TReconCluster);
    cluster->FillFromHits(name,begin,end);

#ifdef CREATE_CLUSTER_SUMMARY_INFORMATION
    // Collect the unique X, V and U hits.
    double summedCharge = 0.0;
    double summedVar = 0.0;
    std::set< CP::THandle<CP::THit> > xHits;
    std::set< CP::THandle<CP::THit> > vHits;
    std::set< CP::THandle<CP::THit> > uHits;
    hitIterator p = begin;
    while (p != end) {
        CP::THandle<CP::TReconHit> rHit = *p;
        if (!rHit) throw CP::EClusterNonHit();
        double sigma = rHit->GetChargeUncertainty();
        summedCharge += rHit->GetCharge();
        summedVar +=  sigma*sigma;
        for (int i=0; i<rHit->GetConstituentCount(); ++i) {
            CP::THandle<CP::THit> c = rHit->GetConstituent(i);
            if (CP::GeomId::Captain::IsUWire(c->GetGeomId())) {
                uHits.insert(rHit->GetConstituent(i));
            }
            else if (CP::GeomId::Captain::IsVWire(c->GetGeomId())) {
                vHits.insert(rHit->GetConstituent(i));
            }
            else if (CP::GeomId::Captain::IsXWire(c->GetGeomId())) {
                xHits.insert(rHit->GetConstituent(i));
            }
        }
        ++p;
    }

    // Find the X charge and it's variance.
    double xCharge = 0.0;
    double xVariance = 0.0;
    for (std::set< CP::THandle<CP::THit> >::iterator h = xHits.begin();
         h != xHits.end(); ++h) {
        xCharge += (*h)->GetCharge();
        double v = (*h)->GetChargeUncertainty();
        xVariance += v*v;
    }

    // Find the V charge and it's variance.
    double vCharge = 0.0;
    double vVariance = 0.0;
    for (std::set< CP::THandle<CP::THit> >::iterator h = vHits.begin();
         h != vHits.end(); ++h) {
        vCharge += (*h)->GetCharge();
        double v = (*h)->GetChargeUncertainty();
        vVariance += v*v;
    }

    // Find the U charge and it's variance.
    double uCharge = 0.0;
    double uVariance = 0.0;
    for (std::set< CP::THandle<CP::THit> >::iterator h = uHits.begin();
         h != uHits.end(); ++h) {
        uCharge += (*h)->GetCharge();
        double v = (*h)->GetChargeUncertainty();
        uVariance += v*v;
    }

    CP::THandle<CP::TClusterState> state = cluster->GetState();

    CaptNamedInfo(
        "CreateCluster", name << " Q: "
        << unit::AsString(state->GetEDeposit(), 
                          std::sqrt(state->GetEDepositVariance()), 
                          "charge")
        << " X: " << unit::AsString(xCharge, 
                                    std::sqrt(xVariance), "charge")
        << " V: " << unit::AsString(vCharge, 
                                    std::sqrt(vVariance), "charge")
        << " U: " << unit::AsString(uCharge, 
                                    std::sqrt(uVariance), "charge"));
#endif
    
    if (recalculateUncertainty) {
        // Adjust the position covariance to say that the position uncertainty
        // is not the weighted average position of the hits, but the moments
        // of the charge distribution.  Empirically, this seems to not give
        // the right chi-squared values for track fits, but gives a better
        // estimate in a shower where the clusters are not done by time
        // sliced.  The default value is set above in the declaration of the
        // CreateCluster template.
        const CP::TReconCluster::MomentMatrix& moments 
            = cluster->GetMoments();
        CP::THandle<CP::TClusterState> covState = cluster->GetState();
        for (int i=0; i<3; ++i) {
            for (int j=0; j<3; ++j) {
                covState->SetPositionCovariance(i,j,moments(i,j));
            }
        }
    }
        
    return cluster;
}
#endif
