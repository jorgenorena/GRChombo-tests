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
    bool radiation_enabled, wall_enabled, handoff_refill_boundaries;
    double thickness, Ro, phi_true, initial_data_residual_tolerance;
    explicit SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        pp.load("radiation_enabled", radiation_enabled, true);
        pp.load("wall_enabled", wall_enabled, true);
        pp.load("handoff_refill_boundaries", handoff_refill_boundaries, false);
        pp.get("radiation_rho0", background.rho0);
        pp.load("radiation_kappa", background.kappa, 1.);
        pp.load("radiation_t0", background.t0, 0.);
        pp.load("G_Newton", background.G, 1.);
        check_parameter("radiation_enabled", radiation_enabled, radiation_enabled && background.rho0 > 0.,
                        "Pass D requires nonzero dynamical radiation");
        const RadiationFLRWBackground bg(background);
        check_parameter("radiation_t0", background.t0, background.t0 == 0., "fresh handoff slice is at t=t0=0");
        double compactness;
        pp.load("thickness", thickness, 0.5);
        pp.load("Ro", Ro, 2.);
        pp.load("compactness_o_Ro", compactness, 0.1);
        pp.load("q", potential_params.q, 0.);
        pp.load("initial_data_residual_tolerance", initial_data_residual_tolerance, 1.e-5);
        check_parameter("wall_parameters", thickness,
            std::isfinite(thickness) && thickness > 0. && std::isfinite(Ro) && Ro > 0. &&
            std::isfinite(compactness) && compactness > 0., "wall scales must be positive and finite");
        check_parameter("initial_data_residual_tolerance", initial_data_residual_tolerance,
            std::isfinite(initial_data_residual_tolerance) && initial_data_residual_tolerance > 0.,
            "must be positive and finite (percent)");
        potential_params.v = std::sqrt(3. / (32. * M_PI) * thickness * compactness);
        potential_params.lambda = 2. / (thickness * thickness * potential_params.v * potential_params.v);
        phi_true = DomainWallPotential::true_vacuum(potential_params.v, potential_params.lambda, potential_params.q);
        check_parameter("restart_file", pp.contains("restart_file"), pp.contains("restart_file"), "Pass D requires fresh GRTresna data");
        check_parameter("max_level", max_level, max_level == 0, "Pass D checks uniform initial-data handoff");
        check_parameter("stop_time", stop_time, stop_time == 0. && max_steps == 0, "Pass D inspects the initial slice without evolving");
        check_parameter("max_spatial_derivative_order", max_spatial_derivative_order,
                        max_spatial_derivative_order == 4, "use fourth-order imported constraints");
        check_parameter("ignore_checkpoint_name_mismatch", ignore_checkpoint_name_mismatch,
                        !ignore_checkpoint_name_mismatch, "exact component names are required");
        check_parameter("activate_extraction", activate_extraction, !activate_extraction, "not part of Pass D");
        for (int d = 0; d < CH_SPACEDIM; ++d)
            if (!boundary_params.is_periodic[d] &&
                (boundary_params.lo_boundary[d] == BoundaryConditions::REFLECTIVE_BC ||
                 boundary_params.hi_boundary[d] == BoundaryConditions::REFLECTIVE_BC))
                check_parameter("radiation_parity", boundary_params.vars_parity[c_psi_rad],
                    boundary_params.vars_parity[c_psi_rad] == BoundaryConditions::EVEN &&
                    boundary_params.vars_parity[c_Pi_rad] == BoundaryConditions::EVEN,
                    "radiation fields must have even parity at reflecting boundaries");
#ifdef USE_AHFINDER
        check_parameter("AH_activate", AH_activate, !AH_activate, "not part of Pass D");
#endif
    }
    std::map<std::string, double> checkpoint_metadata() const
    {
        const RadiationFLRWBackground bg(background);
        return {{"radiation_t0", background.t0}, {"radiation_rho0", background.rho0},
                {"radiation_kappa", background.kappa}, {"G_Newton", background.G},
                {"radiation_H0", bg.H0()}, {"radiation_Pi0", bg.Pi0()},
                {"wall_v", potential_params.v}, {"wall_lambda", potential_params.lambda},
                {"wall_q", potential_params.q}, {"wall_Ro", Ro}, {"wall_thickness", thickness},
                {"initial_sign_of_K", -1.}, {"initial_psi_asymptotic", 1.}};
    }
};
#endif
