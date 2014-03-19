#ifndef CreateCluster_hxx_seen
#define CreateCluster_hxx_seen

#include "HitUtilities.hxx"
#include "ECaptRecon.hxx"

#include <TReconCluster.hxx>
#include <TUnitsTable.hxx>
#include <THandle.hxx>
#include <TReconHit.hxx>
#include <THit.hxx>

namespace CP {

    /// A base exception for the create track template.
    EXCEPTION(ECreateCluster,ECaptRecon);

    /// An exception that there are repeated objects in the track.
    EXCEPTION(EClusterRepeatedObject, ECreateCluster);

    /// An exception that there are repeated objects in the track.
    EXCEPTION(EClusterNonHit, ECreateCluster);

    /// Take iterators from a container holding THandle<THit> objects and
    /// construct a cluster.  The track will be "fitted" using the
    /// TSegmentClusterFit which actually just connects the clusters with
    /// straight line segments.  The resulting track can be refit by
    /// TClusterFit.
    template<typename iterator>
    CP::THandle<CP::TReconCluster> 
    CreateCluster(const char* name, iterator begin, iterator end,
                  bool momentsAsUncertainty = false) {

#ifdef DEBUG_CREATE_CLUSTER
        for (iterator i = begin; i!=end; ++i) {
            CP::THandle< CP::THit > h = *i;
            if (!h) throw CP::EClusterNonHit();
            iterator j = i;
            while ((++j) != end) {
                if (CP::GetPointer(*i) != CP::GetPointer(*j)) continue;
                CaptError("Invalid cluster: multiple copies of object");
                throw CP::EClusterRepeatedObject();
            }
        }
#endif

        CaptNamedInfo("CreateCluster", name );

        CP::THandle<CP::TReconCluster> cluster(new CP::TReconCluster);
        cluster->FillFromHits(name,begin,end);

        // Collect the unique X, V and U hits.
        double summedCharge = 0.0;
        double summedVar = 0.0;
        std::set< CP::THandle<CP::THit> > xHits;
        std::set< CP::THandle<CP::THit> > vHits;
        std::set< CP::THandle<CP::THit> > uHits;
        iterator p = begin;
        while (p != end) {
            CP::THandle<CP::TReconHit> rHit = *p;
            if (!rHit) throw CP::EClusterNonHit();
            double sigma = rHit->GetChargeUncertainty();
            summedCharge += rHit->GetCharge();
            summedVar +=  sigma*sigma;
            xHits.insert(rHit->GetConstituent(0));
            vHits.insert(rHit->GetConstituent(1));
            uHits.insert(rHit->GetConstituent(2));
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

        // Construct the cluster charge from the average of X, V, and U.
        double wireCharge = xCharge/xVariance 
            + vCharge/vVariance
            + uCharge/uVariance;
        double wireVar = 1.0/xVariance + 1.0/vVariance + 1.0/uVariance;
        wireCharge /= wireVar;
        wireVar = 1.0/wireVar; 

        CP::THandle<CP::TClusterState> state = cluster->GetState();
#ifdef RESET_CLUSTER_CHARGE_WITH_WIRE_AVERAGE
        // Update the cluster state with the new charge and it's variance.
        state->SetEDeposit(wireCharge);
        state->SetEDepositVariance(wireVar);
#endif

#define RESET_CLUSTER_CHARGE_WITH_WIRE_SUM
#ifdef RESET_CLUSTER_CHARGE_WITH_WIRE_SUM
        // Update the cluster state with the new charge and it's variance.
        state->SetEDeposit(summedCharge);
        state->SetEDepositVariance(summedVar);
#endif

        CaptNamedInfo(
            "CreateCluster", name << " Q: "
            << unit::AsString(state->GetEDeposit(), 
                              std::sqrt(state->GetEDepositVariance()), 
                              "charge")
            << " on wires: " << unit::AsString(wireCharge, 
                              std::sqrt(wireVar), "charge")
            << " X: " << unit::AsString(xCharge, 
                                        std::sqrt(xVariance), "charge")
            << " V: " << unit::AsString(vCharge, 
                                        std::sqrt(vVariance), "charge")
            << " U: " << unit::AsString(uCharge, 
                                        std::sqrt(uVariance), "charge"));

        if (momentsAsUncertainty) {
            // Adjust the position covariance to say that the position
            // uncertainty is not the weighted average position of the hits,
            // but the moments of the charge distribution.  Empirically, this
            // seems to not give the right chi-squared values for track fits,
            // but gives a better estimate in a shower where the clusters are
            // not done by time sliced.  The default is with
            // momentsAsUncertainty set to be false.
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
};

#endif
