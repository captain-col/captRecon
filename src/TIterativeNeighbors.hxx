#ifndef TIterativeNeighbors_hxx_seen
#define TIterativeNeighbors_hxx_seen

#include <CGAL/Search_traits.h>
#include <CGAL/Orthogonal_incremental_neighbor_search.h>
namespace CP {
    template<typename ValueType> class TIterativeNeighbors;
}


/// A template to quickly find points that are near to a position.  This is
/// typically used to find the neigbors of a hit, but can also lookup a
/// generate xyz position.  This is used like this:
///
/// \code
/// typedef CP::TIterativeNeighbors< CP::THandle<CP::THit> > Neighbors
/// Neighbors neighbors
///
/// for (CP::THitSelection::iterator h = hits.begin(); h!= hits.end(); ++h) {
///    neighbors.AddPoint((*h),
///                       (*h)->GetPosition().X(),
///                       (*h)->GetPosition().Y(),
///                       (*h)->GetPosition().Z());
/// }
/// 
/// for (Neighbors::iterator h 
///          = neighbors.begin(0.0, 0.0, 0.0);
///      h != neighbors.end(); ++h) { 
///     std::cout << "The neighbor " << h->first
///               << " is " << std::sqrt(h->second) << " away"
///               << " std::endl;
/// }
/// \endcode
///
/// The template argument is the type of a value to be attached to each point.
/// Typically, this is a hit or a cluster.
template<typename ValueType>
class CP::TIterativeNeighbors {
public:
    /// An internal representation of the point that is not intended to be
    /// directly used by the caller.
    struct Point {
        /// The type used to represent the point position for quick lookups.
        /// It doesn't need much precision.
        typedef float coord_t;

        /// The type of value being carried by the point.
        typedef ValueType value_t;

        /// A default constructor
        Point() {fVec[0]=fVec[1]=fVec[2]=0.0;}

        /// A constructor for a point
        Point(const value_t& v, double x, double y, double z) : fValue(v)
            {fVec[0]=x; fVec[1]=y; fVec[2]=z;}
        
        /// A constructor from a position.
        Point(double x, double y, double z) 
            {fVec[0]=x; fVec[1]=y; fVec[2]=z;}
        
        /// Used by CGAL to check for exact equality for numeric calculations. 
        bool operator ==(const Point& rhs) {
            for (int i=0; i<3; ++i) if (fVec[i] != rhs.fVec[i]) return false;
            return true;
        }
        bool operator !=(const Point& rhs) {return !(*this == rhs);}

        /// A constructor for iterator used by CGAL to access
        /// the coordinates of this point.
        struct coord_iterator {
            /// The type returned by the coordinate iterator.
            typedef const typename Point::coord_t* result_type;
            result_type operator() (const Point& p) const 
                {return static_cast<result_type>(p.fVec);}
            result_type operator() (const Point& p, int) const 
                {return static_cast<result_type>(p.fVec+3);}
        };

        /// The position information for this point.
        coord_t fVec[3];

        /// The value associated with this point.
        value_t fValue;
    };

    /// A typedef required by CGAL.
    typedef CGAL::Search_traits<
        typename Point::coord_t, 
        Point,
        typename Point::coord_iterator::result_type,
        typename Point::coord_iterator > TreeTraits;
    /// A typedef required by CGAL.
    typedef CGAL::Orthogonal_incremental_neighbor_search<TreeTraits> Neighbors;

    /// An input iterator that returns the values in order of closest to
    /// furthest.  The iterator points to a std::pair<ValueType,float> where
    /// the first element is the value at the point, and the second element is
    /// the distance to the search point provided in the begin method.  The
    /// iterator holds a copy of the ValueType so the original object cannot
    /// be modified (even if you use const_cast.
    class iterator {
    public:
        typedef typename Neighbors::iterator base_iterator;
        typedef std::pair<ValueType, 
                          typename base_iterator::value_type::second_type> 
            value_type;
        typedef std::input_iterator_tag iterator_category;
        iterator() {}
        iterator(const base_iterator& i) : fBase(i) {}
        iterator& operator ++() {++fBase; return *this;} 
        const value_type& operator *() {
            fValue.first = fBase->first.fValue;
            fValue.second = fBase->second;
            return fValue;
        }
        const value_type* operator ->() {
            fValue.first = fBase->first.fValue;
            fValue.second = fBase->second;
            return &fValue;
        }
        bool operator == (const iterator& rhs) const {return fBase==rhs.fBase;}
        bool operator !=(const iterator& rhs) const {return !((*this)==rhs);}
    private:
        base_iterator fBase;
        value_type fValue;
    };
        
    /// An input iterator over all of the values in the neighbor tree.  This
    /// does not interfere with the neighbor iterator.  The pointed to value
    /// is a constant.
    class value_iterator {
    public:
        typedef typename Neighbors::Tree::iterator base_iterator;
        typedef ValueType value_type;
        typedef std::input_iterator_tag iterator_category;
        value_iterator() {}
        value_iterator(const base_iterator& i) : fBase(i) {}
        value_iterator& operator ++() {++fBase; return *this;} 
        const value_type& operator *() {
            return fBase->fValue;
        }
        const value_type* operator ->() {
            return &fBase->fValue;
        }
        bool operator == (const value_iterator& rhs) const {
            return fBase==rhs.fBase;}
        bool operator !=(const value_iterator& rhs) const {
            return !((*this)==rhs);}
    private:
        base_iterator fBase;
    };
        
    virtual ~TIterativeNeighbors() {if (fNeighbors) delete fNeighbors;}

    /// Create an empty set of points that will be filled using the AddPoint()
    /// method.
    TIterativeNeighbors() : fNeighbors(NULL) {}

    /// Add a new point to the neighbor search.  The first argument is the
    /// value associated with the point, and the remaining three are the
    /// cartesian coordinates of the point.  The cartesian coordinates are
    /// used to define the neighbors of this point.
    virtual void AddPoint(const ValueType& v, 
                          double x, double y, double z) {
        fNeighborTree.insert(Point(v,x,y,z));
    }

    /// Get the neighbors.  This returns the neighbors in order of closest to
    /// the furthest.  The returned iterator points to a
    /// std::pair<ValueType,float> where the first element is the value at the
    /// point, and the second element is the distance to the search point
    /// provided in the begin method.  Only one search for neighbors is
    /// allowed at any given time and any iterators are invalidated by the
    /// next call to begin().  A typical loop might look like
    /// \code
    /// CP::TIterativeNeighbors<int>::iterator begin = neighbors.begin(0,0,0);
    /// CP::TIterativeNeighbors<int>::iterator end = neighbors.end();
    /// while (begin != end) {
    ///    std::cout << "Value: " << begin->first 
    ///              << " is " << std::sqrt(begin->second) << " away"
    ///              << std::endl;
    ///    if (begin->second > threshold) break;
    ///    ++begin;
    /// }
    /// \endcode
    virtual iterator begin(double x, double y, double z) {
        if (fNeighbors) delete fNeighbors;
        fNeighbors = new Neighbors(fNeighborTree, Point(x,y,z));
        iterator it(fNeighbors->begin());
        fCurrentEnd = iterator(fNeighbors->end());
        return it;
    }

    /// One past the last neighbor of the point.  This is invalidated by the
    /// next call to begin.
    virtual const iterator& end() const {return fCurrentEnd;}

    /// Return an interator to all of the points in the neighbors search.
    /// This does not invalidate the neighbors iterator (returned by begin()).
    virtual value_iterator begin_values() {
        return fNeighborTree.begin();
    }

    /// Return an interator to all of the points in the neighbors search.
    /// This does not invalidate the neighbors iterator (returned by begin()).
    virtual value_iterator end_values() {
        return fNeighborTree.end();
    }

private:
    /// A tree allocated to hold the mosts recent nearest neighbors.
    Neighbors* fNeighbors;

    /// A k-d tree representation of the point positions.  This allows very
    /// quick lookup of neighbors.
    typename Neighbors::Tree fNeighborTree;

    /// The current end of the neighbor points.  This is actually one past the
    /// furthest point.  This is set by the Neighbors method.
    iterator fCurrentEnd;
};
#endif
