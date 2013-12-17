//
// An unbelievably ugly hack to test private methods.  This needs to be done
// before the class being tested, and is only required if you're testing
// private methods or members.
//
// #define private public
// #define protected public
#include <TPositionNeighbors.hxx>
#undef private
#undef protected

#include <TCaptLog.hxx>

#include <TVector3.h>

#include <tut.h>

namespace tut {
    struct basePositionNeighbors {
        basePositionNeighbors() {
            // Run before each test.
        }
        ~basePositionNeighbors() {
            // Run after each test.
        }
    };

    struct TVector3Handle {
        TVector3Handle(const TVector3& v) : fVector(v) {}
        TVector3& GetPosition() {return fVector;}
        TVector3 fVector;
    };

    // Declare the test
    typedef test_group<basePositionNeighbors>::object testPositionNeighbors;
    test_group<basePositionNeighbors> groupPositionNeighbors(
        "TPositionNeighbors");

    // Test the nearest neighbor search using an int as a payload.
    template<> template<> void testPositionNeighbors::test<1> () {
        typedef CP::TPositionNeighbors<TVector3Handle*> PositionNeighbors;
        PositionNeighbors neighbors;

        double sign = -1;
        for (int i=0; i<5; ++i) {
            sign = -1.0*sign;
            neighbors.AddPoint(NULL,0.5*i,sign,0.0);
        }

        PositionNeighbors::iterator begin = neighbors.begin(0,0,0);
        PositionNeighbors::iterator end = neighbors.end();
        int count = 0;
        while (begin != end) {
            ++count;
            ++begin;
        }
        ensure_equals("All values are seen",count,5);
    }

    // Test the iterator over all of the values in the search.
    template<> template<> void testPositionNeighbors::test<2> () {
        typedef CP::TPositionNeighbors<TVector3Handle*> Neighbors;
        Neighbors neighbors;

        double sign = -1;
        for (int i=0; i<5; ++i) {
            sign = -1.0*sign;
            neighbors.AddPoint(NULL,0.5*i,sign,0.0);
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
    template<> template<> void testPositionNeighbors::test<3> () {
        typedef CP::TPositionNeighbors<TVector3Handle*> Neighbors;
        Neighbors neighbors;

        int numberOfValues = 5;

        double sign = -1;
        for (int i=0; i<numberOfValues; ++i) {
            sign = -1.0*sign;
            TVector3Handle* vec = new TVector3Handle(TVector3(0.5*i,sign,0.0));
            neighbors.AddHandle(vec);
        }

        Neighbors::value_iterator value = neighbors.begin_values();
        Neighbors::value_iterator end_value = neighbors.end_values();
        while (value != end_value) {
            Neighbors::iterator neighbor = neighbors.begin(
                (*value)->GetPosition().X(), 
                (*value)->GetPosition().Y(), 
                (*value)->GetPosition().Z());
            Neighbors::iterator end_neighbor = neighbors.end();
            int count = 0;
            while (neighbor != end_neighbor) {
                if (!count) {
                    ensure_distance("Neighbor X coordinate agrees",
                                    (*value)->GetPosition().X(), 
                                    neighbor->first->GetPosition().X(),
                                    0.001);
                    ensure_distance("Neighbor Y coordinate agrees",
                                    (*value)->GetPosition().Y(), 
                                    neighbor->first->GetPosition().Y(),
                                    0.001);
                    ensure_distance("Neighbor Z coordinate agrees",
                                    (*value)->GetPosition().Z(), 
                                    neighbor->first->GetPosition().Z(),
                                    0.001);

                }
                ++count;
                ++neighbor;
            }
            ensure_equals("Same number of neighbors as values",
                          count, numberOfValues);
            ++value;
        }
    }

    // Test the iterator and value_iterator don't interfere with each other.
    template<> template<> void testPositionNeighbors::test<4> () {
        typedef CP::TPositionNeighbors<TVector3Handle*> Neighbors;
        Neighbors neighbors;
        
        int numberOfValues = 5;
        
        double sign = -1;
        for (int i=0; i<numberOfValues; ++i) {
            sign = -1.0*sign;
            TVector3Handle* vec = new TVector3Handle(TVector3(0.5*i,sign,0.0));
            neighbors.AddHandle(vec);
        }

        Neighbors::value_iterator value = neighbors.begin_values();
        Neighbors::value_iterator end_value = neighbors.end_values();
        while (value != end_value) {
            Neighbors::iterator neighbor = neighbors.begin(*value);
            Neighbors::iterator end_neighbor = neighbors.end();
            int count = 0;
            while (neighbor != end_neighbor) {
                if (!count) {
                    ensure_distance("Neighbor X coordinate agrees",
                                    (*value)->GetPosition().X(), 
                                    neighbor->first->GetPosition().X(), 0.001);
                    ensure_distance("Neighbor Y coordinate agrees",
                                    (*value)->GetPosition().Y(), 
                                    neighbor->first->GetPosition().Y(), 0.001);
                    ensure_distance("Neighbor Z coordinate agrees",
                                    (*value)->GetPosition().Z(), 
                                    neighbor->first->GetPosition().Z(), 0.001);

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
