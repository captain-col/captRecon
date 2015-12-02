#include <TPositionDensityCluster.hxx>
#include <TCaptLog.hxx>
#include <TVector3.h>
#include <tut.h>

#include <memory>

namespace tut {
    struct basePDCluster {
        basePDCluster() {
            // Run before each test.
        }
        ~basePDCluster() {
            // Run after each test.
        }
    };

    struct TVector3Object {
        TVector3Object(const TVector3& v) : fVector(v) {}
        TVector3& GetPosition() {return fVector;}
        TVector3 fVector;
    };
    typedef struct TVector3Object* TVector3Handle;

    // Declare the test
    typedef test_group<basePDCluster>::object testPDCluster;
    test_group<basePDCluster> groupPDCluster(
        "TPositionDensityCluster");

    // Test that the cluster can be constructed.
    template<> template<> void testPDCluster::test<1> () {
        const int minPoints = 5;
        const double maxDist = 1.0;
        CP::TPositionDensityCluster<TVector3Handle> dCluster(minPoints,maxDist);
        ensure_equals("No clusters in when empty",
                      dCluster.GetClusterCount(), 0);
        std::vector<TVector3Handle> vectors;
        dCluster.Cluster(vectors.begin(), vectors.end());
        ensure_equals("No clusters no input objects",
                      dCluster.GetClusterCount(), 0);        
    }

    // Test that the cluster can be constructed.
    template<> template<> void testPDCluster::test<2> () {
        // Make a vector of handles.
        std::vector<TVector3Handle> vectors;
        vectors.push_back(new TVector3Object(TVector3(3,3,3)));
        vectors.push_back(new TVector3Object(TVector3(3.5,3.5,3.5)));
        vectors.push_back(new TVector3Object(TVector3(2.5,3.5,3.5)));
        vectors.push_back(new TVector3Object(TVector3(2.5,2.5,3.5)));
        vectors.push_back(new TVector3Object(TVector3(3.5,2.5,3.5)));
        vectors.push_back(new TVector3Object(TVector3(3.5,3.5,2.5)));
        vectors.push_back(new TVector3Object(TVector3(2.5,3.5,2.5)));
        vectors.push_back(new TVector3Object(TVector3(2.5,2.5,2.5)));
        vectors.push_back(new TVector3Object(TVector3(3.5,2.5,2.5)));
        vectors.push_back(new TVector3Object(TVector3(-3,-3,-3)));
        vectors.push_back(new TVector3Object(TVector3(-3.5,-3.5,-3.5)));
        vectors.push_back(new TVector3Object(TVector3(2.5,-3.5,-3.5)));
        vectors.push_back(new TVector3Object(TVector3(2.5,2.5,-3.5)));
        vectors.push_back(new TVector3Object(TVector3(-3.5,2.5,-3.5)));
        vectors.push_back(new TVector3Object(TVector3(-3.5,-3.5,2.5)));
        vectors.push_back(new TVector3Object(TVector3(2.5,-3.5,2.5)));
        vectors.push_back(new TVector3Object(TVector3(2.5,2.5,2.5)));
        vectors.push_back(new TVector3Object(TVector3(-3.5,2.5,2.5)));
                          
#ifdef PRINT_IT
        for (std::vector<TVector3Handle>::iterator v = vectors.begin();
             v != vectors.end(); ++v) {
            std::cout << " " << (*v)->GetPosition().X()
                      << " " << (*v)->GetPosition().Y()
                      << " " << (*v)->GetPosition().Z()
                      << std::endl;
        }
#endif

        const int minPoints = 5;
        const double maxDist = 1.0;
        CP::TPositionDensityCluster<TVector3Handle> dCluster(minPoints,maxDist);

        dCluster.Cluster(vectors.begin(), vectors.end());

    }

};

// Local Variables:
// mode:c++
// c-basic-offset:4
// End:
