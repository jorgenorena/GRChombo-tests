#include "DefaultLevelFactory.hpp"
#include "GRAMR.hpp"
#include "RadiationFLRWLevel.hpp"
#include "SetupFunctions.hpp"

int main(int argc, char *argv[])
{
    mainSetup(argc, argv);
    GRParmParse pp(argc - 2, argv + 2, NULL, argv[1]);
    SimulationParameters params(pp);
    if (!params.just_check_params)
    {
        GRAMR amr;
        DefaultLevelFactory<RadiationFLRWLevel> factory(amr, params);
        setupAMRObject(amr, factory);
        // Write the initial or restored diagnostics before any evolution.
        for (auto *level : amr.get_gramrlevels())
            level->specificPostTimeStep();
        amr.run(params.stop_time, params.max_steps);
        amr.conclude();
        pout() << "RadiationFLRW finished." << std::endl;
    }
    mainFinalize();
    return 0;
}
