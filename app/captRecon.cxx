#include "TCluster3D.hxx"

#include <eventLoop.hxx>

/// Run the cluster calibration on a file.
class TCaptReconLoop: public CP::TEventLoopFunction {
public:
    TCaptReconLoop() {
        fCaptRecon = NULL;
    }

    virtual ~TCaptReconLoop() {};

    void Usage(void) {     }

    virtual bool SetOption(std::string option,std::string value="") {
        return true;
    }

    bool operator () (CP::TEvent& event) {
        // Make sure the electronics simulated is created.
        if (!fCaptRecon) fCaptRecon = new CP::TCluster3D();

        CP::THandle<CP::THitSelection> drift(event.GetHitSelection("drift"));

        // Run the simulation on the event.
        CP::THandle<CP::TAlgorithmResult> cluster3D
            = fCaptRecon->Process(*drift);
        event.AddFit(cluster3D);

        // Save everything.
        return true;
    }

private:
    CP::TCluster3D* fCaptRecon;
};

int main(int argc, char **argv) {
    TCaptReconLoop userCode;
    CP::eventLoop(argc,argv,userCode);
}
