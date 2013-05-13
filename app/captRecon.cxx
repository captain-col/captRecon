#include "TCaptRecon.hxx"

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
        if (!fCaptRecon) fCaptRecon = new CP::TCaptRecon();

        // Run the simulation on the event.
        (*fCaptRecon)(event);

        // Save everything.
        return true;
    }

private:
    CP::TCaptRecon* fCaptRecon;
};

int main(int argc, char **argv) {
    TCaptReconLoop userCode;
    CP::eventLoop(argc,argv,userCode,1);
}
