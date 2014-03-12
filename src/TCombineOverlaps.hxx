#ifndef TCombineOverlaps_hxx_seen
#define TCombineOverlaps_hxx_seen

#include <TAlgorithm.hxx>
#include <TAlgorithmResult.hxx>

namespace CP {
    class TCombineOverlaps;
};

/// This takes a algorithm result with a TReconObjectContainer of with objects
/// and combines objects that have a high degree of overlaps between their 2D
/// hits.  The individual 3D hits are then density clustered for further
/// processing.  Any object that is not combined with another object is passed
/// through to the output.  The assumption is that this algorithm will be
/// followed by a more detailed fitting algorithm.
class CP::TCombineOverlaps
    : public CP::TAlgorithm {
public:
    TCombineOverlaps();
    virtual ~TCombineOverlaps();

    /// Apply the algorithm.
    CP::THandle<CP::TAlgorithmResult> 
    Process(const CP::TAlgorithmResult& input,
            const CP::TAlgorithmResult& input1 = CP::TAlgorithmResult::Empty,
            const CP::TAlgorithmResult& input2 = CP::TAlgorithmResult::Empty);

    /// Set the overlapping hit cut.
    void SetOverlapCut(double value) {fOverlapCut = value;}

    /// Get the overlapping hit cut.
    double GetOverlapCut() const {return fOverlapCut;}
    
private:

    /// Check the amount of overlap between the hits in the two objectgs.
    double CheckOverlap(CP::THandle<CP::TReconBase> object1,
                        CP::THandle<CP::TReconBase> object2);

    /// Merge the to objects into a single object.
    CP::THandle<CP::TReconBase> 
    MergeObjects(CP::THandle<CP::TReconBase> object1,
                 CP::THandle<CP::TReconBase> object2);

    /// Objects with more than this amount of overlap in the 2D hits will be
    /// combined.
    double fOverlapCut;

};
#endif
