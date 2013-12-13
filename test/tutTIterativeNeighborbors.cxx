///
/// An unbelievably ugly hack to test private methods.  This needs to be done
/// before the class being tested,
///
// #define private public
// #define protected public
#include <TIterativeNeighbors.hxx>
#undef private
#undef protected

#include <TCaptLog.hxx>

#include <TVector3.h>

#include <tut.h>

namespace tut {
    struct baseNeighbors {
        baseNeighbors() {
            // Run before each test.
        }
        ~baseNeighbors() {
            // Run after each test.
        }
    };

    // Declare the test
    typedef test_group<baseNeighbors>::object testNeighbors;
    test_group<baseNeighbors> groupNeighbors("TIterativeNeighbors");

    // Test the nearest neighbor search using an int as a payload.
    template<> template<> void testNeighbors::test<1> () {
        typedef CP::TIterativeNeighbors<int> Neighbors;
        Neighbors neighbors;

        double sign = -1;
        for (int i=0; i<5; ++i) {
            sign = -1.0*sign;
            neighbors.AddPoint(i,0.5*i,sign,0.0);
        }

        Neighbors::iterator begin = neighbors.begin(0,0,0);
        Neighbors::iterator end = neighbors.end();
        int expected = 0;
        while (begin != end) {
            ensure_equals("Payload is in order",begin->first,expected);
            ++begin;
            ++expected;
        }
    }

    // Test the iterator over all of the values in the search.
    template<> template<> void testNeighbors::test<2> () {
        typedef CP::TIterativeNeighbors<int> Neighbors;
        Neighbors neighbors;

        double sign = -1;
        for (int i=0; i<5; ++i) {
            sign = -1.0*sign;
            neighbors.AddPoint(i,0.5*i,sign,0.0);
        }

        Neighbors::value_iterator begin = neighbors.begin_values();
        Neighbors::value_iterator end = neighbors.end_values();
        int count = 0;
        while (begin != end) {
            ++count;
            ++begin;
        }
        ensure_equals("Iteration over values completes",count,5);
    }

    // Test the iterator and value_iterator don't interfere with each other.
    template<> template<> void testNeighbors::test<3> () {
        typedef CP::TIterativeNeighbors<TVector3> Neighbors;
        Neighbors neighbors;

        int numberOfValues = 5;

        double sign = -1;
        for (int i=0; i<numberOfValues; ++i) {
            sign = -1.0*sign;
            TVector3 vec(0.5*i,sign,0.0);
            neighbors.AddPoint(vec,vec.X(), vec.Y(), vec.Z());
        }

        Neighbors::value_iterator value = neighbors.begin_values();
        Neighbors::value_iterator end_value = neighbors.end_values();
        while (value != end_value) {
            Neighbors::iterator neighbor = neighbors.begin(value->X(), 
                                                           value->Y(),
                                                           value->Z());
            Neighbors::iterator end_neighbor = neighbors.end();
            int count = 0;
            while (neighbor != end_neighbor) {
                if (!count) {
                    ensure_distance("Neighbor X coordinate agrees",
                                    value->X(), neighbor->first.X(), 0.001);
                    ensure_distance("Neighbor Y coordinate agrees",
                                    value->X(), neighbor->first.X(), 0.001);
                    ensure_distance("Neighbor Z coordinate agrees",
                                    value->X(), neighbor->first.X(), 0.001);
                }
                ++count;
                ++neighbor;
            }
            ensure_equals("Same number of neighbors as values",
                          count, numberOfValues);
            ++value;
        }
    }

};

// Local Variables:
// mode:c++
// c-basic-offset:4
// End:
