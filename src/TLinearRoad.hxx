#ifndef TLinearRoad_hxx_seen
#define TLinearRoad_hxx_seen

#include <TReconBase.hxx>
#include <TReconCluster.hxx>
#include <TReconTrack.hxx>

#include <TPrincipal.h>

#include <deque>
#include <list>

namespace CP {
    class TLinearRoad;
};

/// A concrete class that implements a simple track following algorithm.
/// Given an object collection of clusters, and a seed track, the algorithm
/// tries to extend the track at both ends.
/// 
/// To use this class
/// \code
/// std::unique_ptr<TLinearRoad> obj(new TLinearRoad(clusters, seed));
/// obj->Process();
/// CP::THandle<CP::TReconTrack> foundTrack = obj->GetTrack();
/// \endcode
class CP::TLinearRoad {
public:
    typedef std::deque< CP::THandle<CP::TReconCluster> > SeedContainer;
    typedef std::list< CP::THandle<CP::TReconCluster> > RemainingContainer;

    /// Construct a new TLinearRoad objects that takes a seed and clusters to
    /// check if it should be expanded.  The seed must be sorted by position
    /// along the candidate track.  A new TLinearRoad object must be
    /// constructed for each track seed.  There is an optional parameter to
    /// limit the number of clusters added to the ends of the track during
    /// road following.
    explicit TLinearRoad(int maxClusters = 100000);
    virtual ~TLinearRoad() { }
    
    /// Expand the track seed in both directions using the clusters.  The
    /// output of this can be gotten using the GetTrack, FillRemains and
    /// FillSeed methods.
    void Process(const CP::TReconObjectContainer& seed,
                 const CP::TReconObjectContainer& clusters);
   
    /// Add the clusters found by Process into a (partially filled)
    /// TReconTrack object.  The track will have the clusters assigned to
    /// nodes and in the "proper" order, but it needs to be fit.
    CP::THandle<CP::TReconTrack> GetTrack();

    /// Get any clusters that were not used in the road.  This mutates the
    /// input TReconObjectContainer (first clears, and the fills).
    void  FillRemains(CP::TReconObjectContainer& remains) const;

    /// Get any clusters that were used in the road.  This mutates the input
    /// TReconObjectContainer (first clears, and the fills).  The FillSeed
    /// method is provided so that diagnostics can be run.  Usually, the
    /// GetTrack method should be used to get the output of the road
    /// following.
    void  FillSeed(CP::TReconObjectContainer& seed) const;

    /// Set the full width (not the half width) of the road to be used in the
    /// road following.
    void SetRoadWidth(double r) {fRoadWidth = r;}

    /// Set the maximum step along the road direction when searching for the
    /// next cluster.  This is in units of the "cluster sizes" and references
    /// the "gap" between the clusters.
    void SetRoadStep(int l) {fRoadStep = l;}

    /// Set the opening angle of the road.  The road opening angle should be a
    /// small value to take into account multiple scattering, and geometric
    /// effects (e.g. wire "stepping" in almost vertical tracks).  This is
    /// actually the tangent of the opening angle (which is the opening angle
    /// in the small angle approximation).
    void SetOpeningAngle(double a) {fOpeningAngle = a;}

    /// Set the maximum number of clusters in the seed used to follow the
    /// road.  A seed must have more than this number of clusters, and be
    /// longer than the seed length.
    void SetSeedSize(int h) {fSeedSize=h;}

    /// Set the maximum distance between the first and last cluster in the
    /// seed used to follow the road.  A seed must be longer than this
    /// length, and have more clusters than the seed size.
    void SetSeedLength(double len) {fSeedLength=len;}

private:
    /// Find the cluster to expand the track.  This assumes that the
    /// neighboring hits have already been found and just selects the best
    /// neighbor.  The clusters in the seed are assumed to be in order from
    /// the front of the SeedContainer to the back.  If extendBack is true
    /// (false), then the next cluster is found relative to the back (front)
    /// of the seed.  The upstream end of the seed is defined as the front,
    /// and the downstream end is defined as the back.
    CP::THandle<CP::TReconCluster> 
    NextCluster(const SeedContainer& seed,
                const RemainingContainer& neighbors,
                bool extendBack = true);


    /// Take a TVector3 and find it's "position" along the principal axis.  
    double FindPositionPrincipal(TPrincipal& pca, const TVector3& pos);

    /// Take a position along the principal axis and find the associated
    /// TVector3.
    TVector3 FindPrincipalPosition(TPrincipal& pca, double principal);

    CP::THandle<CP::TTrackState> 
    CreateTrackState(CP::THandle<CP::TReconCluster> object,
                     const TVector3& direction);

    /// Clusters that belong to this track.
    SeedContainer fTrackClusters;
    
    /// Clusters not currently associated with any track.  This is the source
    /// of clusters to be added to the track.
    RemainingContainer fRemainingClusters;
    
    /// The original clusters in the seed.  This makes sure that in the worst
    /// case, the seed comes out unscathed.
    CP::TReconObjectContainer fOriginalClusters;

    /// The maximum number of clusters to add to the ends of the track during
    /// road following.  The default is unlimited.
    int fMaxClusters;

    /// The road width to be used.
    double fRoadWidth;

    /// The maximum step between clusters as a fraction of the cluster size.
    double fRoadStep;

    /// The opening angle of the road (should be small).
    double fOpeningAngle;

    /// The number of hits in the road following segment. 
    unsigned int fSeedSize;

    /// The length of the road following segment.
    double fSeedLength;
};
#endif
