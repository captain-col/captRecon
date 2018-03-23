#include "TFitting3D.hxx"
#include "TTrackFit.hxx"

#include <THandle.hxx>
#include <TReconTrack.hxx>
#include <TReconCluster.hxx>
#include <TCaptLog.hxx>
#include <TRuntimeParameters.hxx>

#include <memory>
#include <cmath>
#include <set>
#include <list>
#include <algorithm>

/// Fit 3D tracks.
CP::TFitting3D::TFitting3D()
    : TAlgorithm("TFitting3D", 
                 "Fiting 3D") {
   
}

CP::TFitting3D::~TFitting3D() { }


CP::THandle<CP::TAlgorithmResult>
CP::TFitting3D::Process(const CP::TAlgorithmResult& input,
                               const CP::TAlgorithmResult&,
                               const CP::TAlgorithmResult&) {

    CP::THandle<CP::TReconObjectContainer> inputObjects 
        = input.GetResultsContainer();
    CaptLog("TFitting3D Process " << GetEvent().GetContext());

    if (!inputObjects) {
        CaptError("No input objects");
        return CP::THandle<CP::TAlgorithmResult>();
    }

    // Create the output containers.
    CP::THandle<CP::TAlgorithmResult> result = CreateResult();
    std::unique_ptr<CP::TReconObjectContainer> 
        final(new CP::TReconObjectContainer("final"));

    CP::TReconObjectContainer tracks;
    
    for (CP::TReconObjectContainer::iterator tr = inputObjects->begin();
	 tr != inputObjects->end(); ++tr){
      CP::THandle<CP::TReconTrack> track = *tr;
      if (!track) {
	final->push_back(*tr);
	continue;
      }
      tracks.push_back(*tr);
    }

   
for(CP::TReconObjectContainer::iterator tr=tracks.begin();tr!=tracks.end();++tr){
 
  CP::THandle<CP::TReconTrack> trackS=*tr;
                    TTrackFit fitter;
		    CP::THandle<CP::TReconTrack> trackfit= fitter(trackS);

            final->push_back(trackfit);
 }

    result->AddResultsContainer(final.release());
    return result;
}
