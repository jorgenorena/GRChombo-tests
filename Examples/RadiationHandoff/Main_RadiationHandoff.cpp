#include "DefaultLevelFactory.hpp"
#include "GRAMR.hpp"
#include "RadiationHandoffLevel.hpp"
#include "SetupFunctions.hpp"

class HandoffAMR : public GRAMR
{
  public:
    void write_inspected_slice() const { writeCheckpointFile(); writePlotFile(); }
};
int main(int argc, char *argv[])
{
    mainSetup(argc,argv);
    GRParmParse pp(argc-2,argv+2,NULL,argv[1]);
    SimulationParameters params(pp);
    if (!params.just_check_params)
    {
        HandoffAMR amr;
        DefaultLevelFactory<RadiationHandoffLevel> factory(amr,params);
        setupAMRObject(amr,factory);
        amr.get_gramrlevels()[0]->specificPostTimeStep();
        // The generic restart setup assumes the output directory already exists.
        // A handoff reads another code's file and writes to a fresh directory.
        if (!FilesystemTools::directory_exists(params.hdf5_path))
            FilesystemTools::mkdir_recursive(params.hdf5_path);
        amr.write_inspected_slice(); // no run(), no evolution or field reset
        pout() << "Pass D handoff inspection finished; zero evolution steps." << std::endl;
    }
    mainFinalize();
    return 0;
}
