#ifndef TMergeTracks_hxx_seen
#define TMergeTracks_hxx_seen

#include <TAlgorithm.hxx>
#include <TAlgorithmResult.hxx>
#include <TReconTrack.hxx>

namespace CP {
    class TMergeTracks;
    class TReconCluster;
};

/// This takes a algorithm result with a TReconObjectContainer of with tracks
/// and first exams the tracks to see if they should be split because they
/// have sharp corners, or big gaps.
class CP::TMergeTracks
    : public CP::TAlgorithm {
public:
    TMergeTracks();
    virtual ~TMergeTracks();

    /// Apply the algorithm.
    CP::THandle<CP::TAlgorithmResult> 
    Process(const CP::TAlgorithmResult& input,
            const CP::TAlgorithmResult& input1 = CP::TAlgorithmResult::Empty,
            const CP::TAlgorithmResult& input2 = CP::TAlgorithmResult::Empty);

    /// An enum specifiying how the tracks are oriented.
    typedef enum {kNotClose, kFrontFront, kFrontBack, kBackFront, kBackBack} 
        Orientation;

private:
    /// Return how the tracks are oriented relative to each other.
    Orientation TrackOrientation(CP::THandle<CP::TReconTrack> t1,
                                 CP::THandle<CP::TReconTrack> t2);

    /// Return how the positions are oriented relative to each other.
    Orientation TrackOrientation(const TVector3& firstFront,
                                 const TVector3& firstBack,
                                 const TVector3& secondFront,
                                 const TVector3& secondBack);

    /// Return a goodness for how well the tracks are matched.  The goodness
    /// is based on the likelihood that the ends of the tracks line up.  A
    /// higher value is a worse match.
    double MatchGoodness(CP::THandle<CP::TReconTrack> t1,
                         CP::THandle<CP::TReconTrack> t2);

    /// Take two tracks and merge the clusters in the right order to build a
    /// third track which is returned.  If the tracks can't be merged this
    /// will return an empty handle.
    CP::THandle<CP::TReconTrack> MergeTracks(CP::THandle<CP::TReconTrack> t1, 
                                             CP::THandle<CP::TReconTrack> t2);

    /// Return the chi2 for the three clusters falling in a line.  The chi2
    /// will have 3 d.o.f. since there are 9 measurements, and 6 parameters.
    double ThreeInLine(CP::THandle<CP::TReconCluster> a,
                       CP::THandle<CP::TReconCluster> b, 
                       CP::THandle<CP::TReconCluster> c);

    /// The cut value for the chi2 for clusters to be considered as from
    /// the same line.
    double fGoodnessCut;

    /// The cut value for the maximum distance between clusters at the end of
    /// a track to be in the same track.
    double fMergeDistanceCut;

};
#endif
