#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_
#include "GRParmParse.hpp"
#include "Potential.hpp"
#include "RadiationFLRWBackground.hpp"
#include "SimulationParametersBase.hpp"
#include <map>

class SimulationParameters : public SimulationParametersBase
{
  public:
    RadiationFLRWBackground::params_t background;
    Potential::params_t potential_params;
    bool radiation_enabled;
    double phi_true;

    explicit SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        pp.load("radiation_enabled", radiation_enabled, true);
        pp.load("radiation_rho0", background.rho0);
        pp.load("radiation_kappa", background.kappa, 1.);
        pp.load("radiation_t0", background.t0, 0.);
        pp.load("G_Newton", background.G, 1.);
        check_parameter("radiation_rho0", background.rho0,
                        radiation_enabled ? background.rho0 > 0. : background.rho0 == 0.,
                        "must be positive when enabled, or zero when explicitly disabled");
        const RadiationFLRWBackground bg(background);
        bg.at(0.);
        bg.at(stop_time);

        double thickness, compactness_o_Ro;
        pp.load("thickness", thickness, 0.5);
        pp.load("compactness_o_Ro", compactness_o_Ro, 0.25);
        pp.load("q", potential_params.q, 0.);
        check_parameter("thickness", thickness, std::isfinite(thickness) && thickness > 0.,
                        "must be positive and finite");
        check_parameter("compactness_o_Ro", compactness_o_Ro,
                        std::isfinite(compactness_o_Ro) && compactness_o_Ro > 0.,
                        "must be positive and finite");
        potential_params.v = std::sqrt(3. / (32. * M_PI) * thickness * compactness_o_Ro);
        potential_params.lambda = 2. / (thickness * thickness *
                                      potential_params.v * potential_params.v);
        phi_true = DomainWallPotential::true_vacuum(
            potential_params.v, potential_params.lambda, potential_params.q);

        // This executable deliberately implements Pass A only. Additional
        // validation is required before relaxing these restrictions.
        check_parameter("max_level", max_level, max_level == 0,
                        "RadiationFLRW Pass A requires a uniform grid");
        for (int i = 0; i < CH_SPACEDIM; ++i)
            check_parameter("isPeriodic", boundary_params.is_periodic[i],
                            boundary_params.is_periodic[i], "Pass A requires periodic boundaries");
        check_parameter("max_spatial_derivative_order", max_spatial_derivative_order,
                        max_spatial_derivative_order == 4, "Pass A uses fourth-order derivatives");
        check_parameter("ignore_checkpoint_name_mismatch", ignore_checkpoint_name_mismatch,
                        !ignore_checkpoint_name_mismatch, "Pass A requires exact checkpoint names");
        check_parameter("activate_extraction", activate_extraction, !activate_extraction,
                        "extraction is not part of Pass A");
#ifdef USE_AHFINDER
        check_parameter("AH_activate", AH_activate, !AH_activate,
                        "apparent-horizon searches are not part of Pass A");
#endif
        pout() << "Radiation FLRW reference: t0=" << background.t0
               << " rho0=" << background.rho0 << " H0=" << bg.H0()
               << " Pi0=" << bg.Pi0() << " kappa=" << background.kappa << std::endl;
    }

    std::map<std::string, double> checkpoint_metadata() const
    {
        const RadiationFLRWBackground bg(background);
        return {{"radiation_t0", background.t0}, {"radiation_rho0", background.rho0},
                {"radiation_kappa", background.kappa}, {"G_Newton", background.G},
                {"radiation_H0", bg.H0()}, {"radiation_Pi0", bg.Pi0()},
                {"wall_v", potential_params.v}, {"wall_lambda", potential_params.lambda},
                {"wall_q", potential_params.q}};
    }
};
#endif
