#include "TDestroyShortTracks.hxx"


#include <THandle.hxx>
#include <TReconTrack.hxx>
#include <TReconCluster.hxx>
#include <TCaptLog.hxx>
#include <HEPUnits.hxx>
#include <TUnitsTable.hxx>
#include <TRuntimeParameters.hxx>
 
#include <memory>
#include <set>
#include <cmath>

CP::TDestroyShortTracks::TDestroyShortTracks()
    : TAlgorithm("TDestroyShortTracks", 
                 "Break up short tracks into separate hits") {
}

CP::TDestroyShortTracks::~TDestroyShortTracks() { }

CP::THandle<CP::TAlgorithmResult>
CP::TDestroyShortTracks::Process(const CP::TAlgorithmResult& input,
                               const CP::TAlgorithmResult&,
                               const CP::TAlgorithmResult&) {
    
    CP::THandle<CP::TReconObjectContainer> inputObjects 
        = input.GetResultsContainer();

    CaptLog("TDestroyShortTracks Process " << GetEvent().GetContext());

    if (!inputObjects) {
        CaptError("No input objects");
        return CP::THandle<CP::TAlgorithmResult>();
    }

    // Create the output containers.
    CP::THandle<CP::TAlgorithmResult> result = CreateResult();
    std::unique_ptr<CP::TReconObjectContainer> 
        final(new CP::TReconObjectContainer("final"));

    // Divide into objects to save and objects to disintegrate.
    for (CP::TReconObjectContainer::iterator it = inputObjects->begin();
         it != inputObjects->end(); ++it) {
        CP::THandle<CP::TReconTrack> track = *it;
        if (track) {
            	double minLength=15*unit::mm;
		    double length = (track->GetFront()->GetPosition().Vect()-track->GetBack()->GetPosition().Vect()).Mag();
		    if(length > minLength )
		      {
		        final->push_back(*it);	
		      }  
		  }
        }
    

  
    result->AddResultsContainer(final.release());

    return result;
}
