#ifndef TSplitTracks_hxx_seen
#define TSplitTracks_hxx_seen
#include <TAlgorithm.hxx>
#include <TAlgorithmResult.hxx>

namespace CP {
    class TSplitTracks;
    class TReconCluster;
};

/// This takes a algorithm result with a TReconObjectContainer of with tracks
/// and first exams the tracks to see if they should be split because they
/// have sharp corners, or big gaps.
class CP::TSplitTracks
    : public CP::TAlgorithm {
public:
    TSplitTracks();
    virtual ~TSplitTracks();

    /// Apply the algorithm.
    CP::THandle<CP::TAlgorithmResult> 
    Process(const CP::TAlgorithmResult& input,
            const CP::TAlgorithmResult& input1 = CP::TAlgorithmResult::Empty,
            const CP::TAlgorithmResult& input2 = CP::TAlgorithmResult::Empty);

    /// The container type for the clusters being built into tracks.  This
    /// needs to have a random acces iterator.
    typedef std::vector< CP::THandle<CP::TReconCluster> > ClusterContainer;

private:

    /// Save a track to the final object container, but first check that it
    /// makes since.
    void SaveTrack(CP::TReconObjectContainer& container,
                   ClusterContainer::iterator begin, 
                   ClusterContainer::iterator end);

    /// Return the chi2 for the three clusters falling in a line.  The chi2
    /// will have 3 d.o.f. since there are 9 measurements, and 6 parameters.
    double ThreeInLine(CP::THandle<CP::TReconCluster> a,
                       CP::THandle<CP::TReconCluster> b, 
                       CP::THandle<CP::TReconCluster> c) const;

    /// Return the radius of a circle passing through the three clusters.  The
    /// clusters must be oriented with b between a and c.  If the angle abc is
    /// more than 90 degrees, the radius is returned as zero (this is
    /// equivalent to b not being between a and c).
    double RadiusOfCurvature(CP::THandle<CP::TReconCluster> a,
                             CP::THandle<CP::TReconCluster> b, 
                             CP::THandle<CP::TReconCluster> c) const;


    /// Return the size of a cluster transverse to a particular direction.
    double Transversity(CP::THandle<CP::TReconCluster> cluster,
                        const TVector3& dir) const;

    // Return the maximum width of a track.
    double  MaxTransversity(const CP::TReconNodeContainer& nodes) const;

    // Return the maximum width of a track.
    double  MaxTransversity(ClusterContainer::iterator begin,
                            ClusterContainer::iterator end) const;

    /// Return the length of a group of clusters.  The clusters are assumed to
    /// represent a track.
    double TrackLength(const CP::TReconNodeContainer& nodes) const;
    double TrackLength(ClusterContainer::iterator begin,
                       ClusterContainer::iterator end) const;
                       
    /// Return the kink angle at an iterator.  The begin and end iterator
    /// define the entire cluster container.  The kink angle will be zero for
    /// a perfectly straight track, and is calculated based on a minimum of 3
    /// clusters to each side.
    double KinkAngle(ClusterContainer::iterator here, 
                     ClusterContainer::iterator begin,
                     ClusterContainer::iterator end) const;

    /// The cut value for the maximum distance between clusters in at the end
    /// of a track that are in the same track.  Tracks that have gaps bigger
    /// than this may be remerged later, but are split here.
    double fEndDistanceCut;

    /// The maximum ratio between the transverse radius of two clusters at the
    /// end of a track.  If the end cluster is fatter than this, it is
    /// probably part of a horizontal track that is very fat.
    double fTransversityCut;
    
    /// The cut value for the chi2 for three clusters to be considered as from
    /// the same line.
    double fThreeInLineCut;

    /// The minimum radius of curvature for a track.  If the track has a
    /// radius of curvature bigger than this, then it is OK.
    double fRadiusOfCurvature;

    /// The cut value for the maximum distance between clusters in the middle
    /// of a track that are in the same track.  Tracks that have gaps bigger
    /// than this may be remerged later, but are split here.
    double fSplitDistanceCut;

    /// The maximum angle between different segments at a kink.
    double fKinkAngleCut;

};
#endif
