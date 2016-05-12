#ifndef TCaptainRecon_hxx_seen
#define TCaptainRecon_hxx_seen

#include <TAlgorithm.hxx>
#include <TAlgorithmResult.hxx>

#include <memory>

namespace CP {
    class TCaptainRecon;
};

/// Apply the CAPTAIN reconstruction to a TAlgorithmResult handle
/// containing 2D wire hits.  Since a THitSelection can be converted into
/// a TAlgorithmResult this can also be called with a THandle to a
/// THitSelection.  The output TAlgorithmResult will contain:
///
///   * used -- A hit selection of the hits that are used in the final
///                  object.
/// 
///   * unused -- A hit selection of the hits that were not used in the
///                  final object.
///
///   * final -- The final reconstruction objectgs for this algorithm
/// 
/// The hit selection in the first input algorithm result is expected to be
/// the wire hits.  The hit selection in the second input algorithm result is
/// expected to be the PMT hits.
class CP::TCaptainRecon: public CP::TAlgorithm {
public:
    TCaptainRecon();
    virtual ~TCaptainRecon();

    /// Actually apply the reconstruction algorithms.
    CP::THandle<CP::TAlgorithmResult> 
    Process(const CP::TAlgorithmResult& input,
            const CP::TAlgorithmResult& input1 = CP::TAlgorithmResult::Empty,
            const CP::TAlgorithmResult& input2 = CP::TAlgorithmResult::Empty);


private:

    /// A template to simplify calling sub-algorithms.  This handles the
    /// TAlgorithm memory management.
    template<typename T>
    CP::THandle<CP::TAlgorithmResult> Run(const CP::TAlgorithmResult& in) {
        std::unique_ptr<T> ptr(new T);
        return ptr->Process(in);
    }

    /// A template to simplify calling sub-algorithms.  This handles the
    /// TAlgorithm memory management.
    template<typename T>
    CP::THandle<CP::TAlgorithmResult> Run(const CP::TAlgorithmResult& in1,
                                          const CP::TAlgorithmResult& in2) {
        std::unique_ptr<T> ptr(new T);
        return ptr->Process(in1,in2);
    }
};
#endif
