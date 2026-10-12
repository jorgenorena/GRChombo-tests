#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_
#include "AcousticGeometry.hpp"
#include "GRParmParse.hpp"
#include "Potential.hpp"
#include "SimulationParametersBase.hpp"
#include <map>

class SimulationParameters : public SimulationParametersBase
{
  public:
    double G_Newton, radiation_rho0, radiation_kappa, Pi0;
    double wave_amplitude, phi_true;
    double wavelength, static_refine_lo, static_refine_hi;
    int tagging_mode; // 0: radiation gradient; 1: fixed slab
    bool standing_wave, mirror_tags;
    AcousticGeometry geometry;
    Potential::params_t potential_params;

    explicit SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        pp.load("G_Newton", G_Newton, 0.);
        pp.load("radiation_rho0", radiation_rho0, 0.75);
        pp.load("radiation_kappa", radiation_kappa, 1.);
        pp.load("wave_amplitude", wave_amplitude, 1.e-5);
        pp.load("map_amplitude", geometry.amplitude, 0.);
        pp.load("wavelength", wavelength, 1.);
        pp.load("standing_wave", standing_wave, false);
        pp.load("tagging_mode", tagging_mode, 0);
        pp.load("static_refine_lo", static_refine_lo, 0.25);
        pp.load("static_refine_hi", static_refine_hi, 0.75);
        pp.load("mirror_tags", mirror_tags, false);
        check_parameter("wavelength", wavelength, std::isfinite(wavelength) && wavelength > 0., "must be positive");
        geometry.k = 2. * M_PI / wavelength;
        check_parameter("map_amplitude", geometry.amplitude, geometry.amplitude == 0., "Pass C uses Cartesian flat geometry");
        check_parameter("tagging_mode", tagging_mode, tagging_mode == 0 || tagging_mode == 1, "0=gradient, 1=fixed slab");
        check_parameter("static_refine_hi", static_refine_hi, static_refine_hi > static_refine_lo, "empty fixed refinement interval");
        check_parameter("regrid_threshold", regrid_thresholds[0], regrid_thresholds[0] > 0., "must be positive so a homogeneous state is not tagged");
        check_parameter("G_Newton", G_Newton, G_Newton == 0.,
                        "Pass C uses a fixed Minkowski metric and requires G_Newton=0");
        check_parameter("radiation_rho0", radiation_rho0,
                        std::isfinite(radiation_rho0) && radiation_rho0 > 0.,
                        "must be positive and finite");
        check_parameter("radiation_kappa", radiation_kappa,
                        std::isfinite(radiation_kappa) && radiation_kappa > 0.,
                        "must be positive and finite");
        Pi0 = std::pow(4. * radiation_rho0 / (3. * radiation_kappa), 0.25);
        check_parameter("Pi0", Pi0, std::isfinite(Pi0) && Pi0 > 0.,
                        "derived radiation momentum must be positive and finite");
        check_parameter("wave_amplitude", wave_amplitude,
                        std::isfinite(wave_amplitude) && std::abs(wave_amplitude) < 0.5,
                        "must have magnitude <0.5 for timelike initial data");
        check_parameter("map_amplitude", geometry.amplitude,
                        std::isfinite(geometry.amplitude) && std::abs(geometry.amplitude) < 0.5,
                        "must have magnitude <0.5 for a regular map and resolved coordinate speeds");
        check_parameter("dt_multiplier", dt_multiplier,
                        dt_multiplier <= 0.25 * (1. - std::abs(geometry.amplitude)),
                        "must resolve the smallest mapped physical cell (<=0.25*(1-|map_amplitude|))");

        double thickness, compactness;
        pp.load("thickness", thickness, 0.5);
        pp.load("compactness_o_Ro", compactness, 0.25);
        pp.load("q", potential_params.q, 0.);
        check_parameter("thickness", thickness, std::isfinite(thickness) && thickness > 0.,
                        "must be positive and finite");
        check_parameter("compactness_o_Ro", compactness,
                        std::isfinite(compactness) && compactness > 0.,
                        "must be positive and finite");
        potential_params.v = std::sqrt(3. / (32. * M_PI) * thickness * compactness);
        potential_params.lambda = 2. / (thickness * thickness *
                                      potential_params.v * potential_params.v);
        phi_true = DomainWallPotential::true_vacuum(
            potential_params.v, potential_params.lambda, potential_params.q);
        check_parameter("max_level", max_level, max_level >= 0 && max_level <= 1,
                        "Pass C permits at most one refinement level");
        for (int i = 0; i < CH_SPACEDIM; ++i)
        {
            const bool reflected = boundary_params.lo_boundary[i] == BoundaryConditions::REFLECTIVE_BC &&
                                   boundary_params.hi_boundary[i] == BoundaryConditions::REFLECTIVE_BC;
            check_parameter("isPeriodic", boundary_params.is_periodic[i],
                            boundary_params.is_periodic[i] || reflected, "use periodic or two reflective faces");
            if (!boundary_params.is_periodic[i])
            {
                check_parameter("radiation_parity", boundary_params.vars_parity[c_psi_rad],
                                boundary_params.vars_parity[c_psi_rad] == BoundaryConditions::EVEN &&
                                boundary_params.vars_parity[c_Pi_rad] == BoundaryConditions::EVEN,
                                "psi_rad and Pi_rad must have even reflection parity");
                if (i == 0)
                {
                    check_parameter("standing_wave", standing_wave, standing_wave, "reflecting x requires a standing wave");
                    const double half_waves = 2. * (ivN[0] + 1) * coarsest_dx / wavelength;
                    check_parameter("wavelength", wavelength, std::abs(half_waves-std::round(half_waves)) < 1.e-12,
                                    "reflection faces must lie at standing-wave extrema");
                }
            }
        }
        check_parameter("max_spatial_derivative_order", max_spatial_derivative_order,
                        max_spatial_derivative_order == 4, "Pass C uses fourth-order derivatives");
        check_parameter("ignore_checkpoint_name_mismatch", ignore_checkpoint_name_mismatch,
                        !ignore_checkpoint_name_mismatch, "exact component names are required");
        check_parameter("activate_extraction", activate_extraction, !activate_extraction,
                        "wave extraction is not used in this test driver");
#ifdef USE_AHFINDER
        check_parameter("AH_activate", AH_activate, !AH_activate,
                        "Pass C does not perform horizon searches");
#endif
        pout() << "Pass C: fixed flat geometry, G=0, Pi0=" << Pi0
               << " wave amplitude=" << wave_amplitude
               << " map amplitude=" << geometry.amplitude << std::endl;
    }

    std::map<std::string, double> checkpoint_metadata() const
    {
        return {{"G_Newton", G_Newton}, {"radiation_rho0", radiation_rho0},
                {"radiation_kappa", radiation_kappa}, {"radiation_Pi0", Pi0},
                {"wave_amplitude", wave_amplitude}, {"wave_k", geometry.k},
                {"map_amplitude", geometry.amplitude}, {"standing_wave", standing_wave},
                {"tagging_mode", tagging_mode}, {"tag_threshold", regrid_thresholds[0]},
                {"static_refine_lo", static_refine_lo}, {"static_refine_hi", static_refine_hi}, {"mirror_tags", mirror_tags}, {"wall_v", potential_params.v},
                {"wall_lambda", potential_params.lambda}, {"wall_q", potential_params.q}};
    }
};
#endif
