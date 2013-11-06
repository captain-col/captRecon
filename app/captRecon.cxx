#include "TCaptainRecon.hxx"

#include <eventLoop.hxx>

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

    bool operator () (CP::TEvent& event) {
        // Make sure the electronics simulated is created.
        CP::THandle<CP::THitSelection> drift(event.GetHits("drift"));
        CP::THandle<CP::THitSelection> pmt(event.GetHits("pmt"));

        /// An event without any drift hits.  Just pass it through.
        if (!drift) return true;

        /// An event without any PMT hits.  Just pass it through.
        if (!pmt) return true;

        // Run the event reconstruction on the event.
        std::auto_ptr<CP::TAlgorithm> captRecon(new CP::TCaptainRecon());
        CP::THandle<CP::TAlgorithmResult> result 
            = captRecon->Process(*drift,*pmt);
        if (result) event.AddFit(result);
        else {
            CaptError("No reconstruction result");
        }

        // Save everything.
        return true;
    }

};

int main(int argc, char **argv) {
    TCaptReconLoop userCode;
    CP::eventLoop(argc,argv,userCode);
}
