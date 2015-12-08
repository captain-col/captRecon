#ifndef TRemoveOutliers_hxx_seen
#define TRemoveOutliers_hxx_seen

#include <THitSelection.hxx>
#include <TReconHit.hxx>

#include <map>
#include <set>

namespace CP {
    class TRemoveOutliers;
};

/// Remove 3D hits that are "outliers" for the constituient 2D hits.  This
/// takes the hypothesis that when a wire is hit at a particular time (ie a 2D
/// hit), the charge deposition on that wire is actually fairly well localized
/// along the wire.  Additionally, it's fairly unlikely that there will be two
/// clumps of charge deposition along the wire in the same time slice (ie
/// charge from two particles that are well separated in space don't usually
/// arrive at the wire at the same time).  Obviously, this isn't true in all
/// cases, but does cover most "reasonable" event geometries.  In one time
/// slice and for one wire, this uses the other wire planes to estimate where
/// the charge is located along that wire.  
class CP::TRemoveOutliers {
public:

    struct HitInfo {
        // The 3D hits that contain the hit.
        CP::THitSelection fContainedBy;
    };

    typedef std::map< CP::THandle<CP::THit>, HitInfo> HitMap;

    /// Start outlier removal using an input hit selection that contains
    /// TReconHit objects.  This hits in the hit selection will be modified.
    /// See the fMaxOutlierSize, fSmallGroupSize, and fBigGroupSize for the
    /// parameter documentation.
    TRemoveOutliers(int maxOutlier=4, double small=0.1, double big=0.2);
    
    // Mark the 3D hits that are outliers and remove them from the input THit
    // selection.
    void Apply(CP::THitSelection& hits); 

private:

    /// Determine if a particular TReconHit is an outlier.  The index is the
    /// index of the constituent to check.
    bool IsOutlier(CP::THandle<CP::THit> hit, int index);

    /// Find the wires that are in the same group as the wire.
    int LocalGroup(int wire, std::set<int>& input);

    /// A map from the 2D hits back to the 3D hits.
    HitMap fHitMap;

    /// The maximum size for an outlier.  Any group bigger than this isn't an
    /// outlier.
    int fMaxOutlierSize;

    /// If the size of the group that a wire is in is less than this fraction
    /// of the total number of cross wires, then the wire is in an outlier.
    double fSmallGroupFraction;

    /// If the size of both groups that contain the 3D recon hit are less than
    /// this fraction of the total number of cross wires, then the wire is an
    /// outlier.
    double fBigGroupFraction;
};

#endif
