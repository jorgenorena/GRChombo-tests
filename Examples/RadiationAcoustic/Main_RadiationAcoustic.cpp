#include "DefaultLevelFactory.hpp"
#include "GRAMR.hpp"
#include "RadiationAcousticLevel.hpp"
#include "SetupFunctions.hpp"

int main(int argc, char *argv[])
{
    mainSetup(argc, argv);
    GRParmParse pp(argc - 2, argv + 2, NULL, argv[1]);
    SimulationParameters params(pp);
    if (!params.just_check_params)
    {
        GRAMR amr;
        DefaultLevelFactory<RadiationAcousticLevel> factory(amr, params);
        setupAMRObject(amr, factory);
        for (auto *level : amr.get_gramrlevels())
            level->specificPostTimeStep();
        amr.run(params.stop_time, params.max_steps);
        amr.conclude();
        pout() << "RadiationAcoustic finished." << std::endl;
    }
    mainFinalize();
    return 0;
}
