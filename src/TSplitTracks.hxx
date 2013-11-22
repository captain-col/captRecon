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
                       CP::THandle<CP::TReconCluster> c);

    /// The cut value for the chi2 for three clusters to be considered as from
    /// the same line.
    double fThreeInLineCut;

    /// The cut value for the maximum distance between clusters in the middle
    /// of a track that are in the same track.  Tracks that have gaps bigger
    /// than this may be remerged later, but are split here.
    double fSplitDistanceCut;

    /// The cut value for the maximum distance between clusters in at the end
    /// of a track that are in the same track.  Tracks that have gaps bigger
    /// than this may be remerged later, but are split here.
    double fEndDistanceCut;

};
#endif
