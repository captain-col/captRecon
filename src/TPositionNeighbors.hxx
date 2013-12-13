#ifndef TPositionNeighbors_hxx_seen
#define TPositionNeighbors_hxx_seen

#include "TIterativeNeighbors.hxx"

#include <THandle.hxx>
#include <THit.hxx>

namespace CP {
    template <typename PositionHandle> class TPositionNeighbors;
}

/// This is a light wrapper around the TIterativeNeighbors template that
/// provides quick lookup to find the neighbors of a hit, or cluster.  This
/// works for any class that implements the PositionHandle concept.  This
/// concept means that the object looks like a pointer to a class that
/// implements a GetPosition() method.  The GetPosition() method must then
/// have sub-methods for X(), Y(), and Z().  Examples of objects that
/// implement the PositionHandle concept are THandle<THit> and
/// THandle<TReconCluster>.  This inherits all of the TIterativeNeighbors
/// public methods.
template <typename PositionHandle>
class CP::TPositionNeighbors 
    : public CP::TIterativeNeighbors< PositionHandle > {
public:
    using CP::TIterativeNeighbors< PositionHandle >::begin;

    typedef 
    typename CP::TIterativeNeighbors< PositionHandle >::iterator  iterator;

    virtual ~TPositionNeighbors();

    /// Construct an empty object (Handles are added using AddHandle(), or
    /// AddPoint()).
    TPositionNeighbors() {}

    /// A templated constructor to take iterators to the PositionHandle
    /// and add them the neighbor tree.
    template<typename HandleIterator>
    TPositionNeighbors(HandleIterator begin, HandleIterator end) {
        while (begin != end) {
            AddHandle(*begin);
            ++begin;
        }
    }

    /// Add a specialization of TIterativeNeighbors::AddPoint() to make it
    /// easy to add a handle.  The handle will usually be a THandle<THit>, or
    /// a THandle<TReconCluster>
    virtual void AddHandle(PositionHandle handle);

    /// Start the search for a neighbor based on a handle (a specialization of
    /// TIterativeNeighbors::begin()).  The handle will usually be a
    /// THandle<THit> or a THandle<TReconCluster>.  If you want to find the
    /// closest handle to a position, you can still use the
    /// TIterativeNeighbors::begin(double,double,double) method.
    virtual iterator begin(PositionHandle handle);
};

template <typename PositionHandle>
CP::TPositionNeighbors<PositionHandle>::~TPositionNeighbors() {}

template <typename PositionHandle>
void CP::TPositionNeighbors<PositionHandle>::AddHandle(PositionHandle handle) {
    CP::TIterativeNeighbors<PositionHandle>
        ::AddPoint(handle,
                   handle->GetPosition().X(), 
                   handle->GetPosition().Y(), 
                   handle->GetPosition().Z());
}

template <typename PositionHandle>
typename CP::TPositionNeighbors<PositionHandle>::iterator 
CP::TPositionNeighbors<PositionHandle>::begin(PositionHandle handle) {
    return begin(handle->GetPosition().X(),
                 handle->GetPosition().Y(),
                 handle->GetPosition().Z());
}
#endif
