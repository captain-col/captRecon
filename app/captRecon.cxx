#include "TCaptainRecon.hxx"

#include <eventLoop.hxx>

#include <TH1F.h>
#include <TPad.h>

#include <memory>

/// Run the cluster calibration on a file.
class TCaptReconLoop: public CP::TEventLoopFunction {
public:
  TCaptReconLoop() {}

  virtual ~TCaptReconLoop() {};

  void Usage(void) {     }

  virtual bool SetOption(std::string option,std::string value="") {
    return true;
  }

  void Initialize() {
    XUV_tracks = new TH1F("XUV_tracks","",10,0,10);    
  }

  
  bool operator () (CP::TEvent& event) {
    // Make sure the electronics simulated is created.
    CP::THandle<CP::THitSelection> drift(event.GetHits("drift"));
    CP::THandle<CP::THitSelection> pmt(event.GetHits("pmt"));

    /// An event without any drift hits.  Just pass it through.
    if (!drift) return true;

    // Run the event reconstruction on the event.
    std::auto_ptr<CP::TAlgorithm> captRecon(new CP::TCaptainRecon());
    CP::THandle<CP::TAlgorithmResult> result;
    if (pmt) result = captRecon->Process(*drift,*pmt);
    else result = captRecon->Process(*drift);
    if (result) {
      event.AddFit(result);
      if(result->GetResultsContainer("match3Tr")){
	//std::cout<<"TEST2="<<(result->GetResultsContainer("match3Tr"))->size()<<std::endl;
	XUV_tracks->Fill((result->GetResultsContainer("match3Tr"))->size());
    }
    }
    else {
      CaptError("No reconstruction result");
    }

    // Save everything.
    return true;
  }

  void Finalize(CP::TRootOutput * const output) {    
    XUV_tracks->Draw();
    gPad->Print("plots/XUV_tracks.root");
  }

private:
  TH1F* XUV_tracks;// = new TH1F("XUV_tracks","",50,0,50);
  
};



int main(int argc, char **argv) {
  TCaptReconLoop userCode;
  CP::eventLoop(argc,argv,userCode);
}
