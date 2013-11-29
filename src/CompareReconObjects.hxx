#ifndef CompareReconObjects_hxx_seen
#define CompareReconObjects_hxx_seen

#include "ECaptRecon.hxx"

#include <TReconCluster.hxx>
#include <TReconTrack.hxx>
#include <TReconPID.hxx>
#include <TReconVertex.hxx>
#include <THandle.hxx>

namespace CP {

    /// This is a predicate to order recon objects.  This puts the vertices,
    /// followed by pids, tracks and clusters.  Within any type the objects
    /// are ordered by size.
    struct CompareReconObjects {
        bool operator () (CP::THandle<CP::TReconBase> lhs, 
                          CP::THandle<CP::TReconBase> rhs) {
            CP::THandle<CP::TReconVertex> lv = lhs;
            CP::THandle<CP::TReconVertex> rv = rhs;
            CP::THandle<CP::TReconPID> lp = lhs;
            CP::THandle<CP::TReconPID> rp = rhs;
            CP::THandle<CP::TReconTrack> lt = lhs;
            CP::THandle<CP::TReconTrack> rt = rhs;
            CP::THandle<CP::TReconCluster> lc = lhs;
            CP::THandle<CP::TReconCluster> rc = rhs;
            
            if (lv && rv) return lv->GetNodes().size() > rv->GetNodes().size();
            if (lv && !rv) return true;
            if (!lv && rv) return false;

            if (lp && rp) return lp->GetNodes().size() > rp->GetNodes().size();
            if (lp && !rp) return true;
            if (!lp && rp) return false;

            if (lt && rt) return lt->GetNodes().size() > rt->GetNodes().size();
            if (lt && !rt) return true;
            if (!lt && rt) return false;

            if (lc && rc) return lc->GetEDeposit() > rc->GetEDeposit();
            if (lc && !rc) return true;
            if (!lc && rc) return false;

            return CP::GetPointer(lhs) < CP::GetPointer(rhs);
        }
    };
};
#endif
