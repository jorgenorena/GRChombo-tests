#include "RadiationFLRWLevel.hpp"
#include "AMRReductions.hpp"
#include "BoxLoops.hpp"
#include "ComputePack.hpp"
#include "CosmoMovingPunctureGauge.hpp"
#include "NanCheck.hpp"
#include "NewMatterConstraints.hpp"
#include "SetValue.hpp"
#include "SmallDataIO.hpp"
#include "TraceARemoval.hpp"
#include <fstream>

namespace
{
class InitialFLRW
{
    RadiationFLRWBackground::state_t m_bg;
    double m_phi;
  public:
    InitialFLRW(RadiationFLRWBackground::state_t bg, double phi) : m_bg(bg), m_phi(phi) {}
    template <class T> void compute(Cell<T> cell) const
    {
        cell.store_vars(T(m_bg.chi), c_chi);
        cell.store_vars(T(1.), c_h11);
        cell.store_vars(T(1.), c_h22);
        cell.store_vars(T(1.), c_h33);
        cell.store_vars(T(m_bg.K), c_K);
        cell.store_vars(T(1.), c_lapse);
        cell.store_vars(T(m_phi), c_phi);
        cell.store_vars(T(m_bg.psi), c_psi_rad);
        cell.store_vars(T(m_bg.Pi), c_Pi_rad);
    }
};
}

void RadiationFLRWLevel::initialData()
{
    BoxLoops::loop(make_compute_pack(SetValue(0.),
        InitialFLRW(RadiationFLRWBackground(m_p.background).at(m_time), m_p.phi_true)),
        m_state_new, m_state_new, INCLUDE_GHOST_CELLS);
    refresh_diagnostics();
}

void RadiationFLRWLevel::refresh_diagnostics()
{
    fillAllGhosts();
    BoxLoops::loop(RadiationStateCheck(m_dx, m_p.radiation_enabled),
                  m_state_new, m_state_new, EXCLUDE_GHOST_CELLS, disable_simd());
    const auto model = matter();
    BoxLoops::loop(make_compute_pack(
        MatterConstraints<WallRadiation>(model, m_dx, m_p.background.G,
                                         c_Ham, Interval(c_Mom, c_Mom)),
        FLRWDiagnostics(model, m_dx, m_p.background.kappa, m_p.radiation_enabled)),
        m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
}

void RadiationFLRWLevel::postRestart()
{
    // Every restart time, including positive times: read the restored fields,
    // without replacing any of them by the analytic homogeneous reference.
    refresh_diagnostics();
    pout() << "Reconstructed radiation/total constraints at restart t=" << m_time
           << "; original reference t0=" << m_p.background.t0 << std::endl;
}

void RadiationFLRWLevel::specificEvalRHS(GRLevelData &soln, GRLevelData &rhs,
                                       double time)
{
    BoxLoops::loop(TraceARemoval(), soln, soln, INCLUDE_GHOST_CELLS);
    BoxLoops::loop(RadiationStateCheck(m_dx, m_p.radiation_enabled),
                  soln, soln, EXCLUDE_GHOST_CELLS, disable_simd());
    CosmoMovingPunctureGauge gauge(m_p.ccz4_params);
    gauge.set_K_mean(RadiationFLRWBackground(m_p.background).at(time).K);
    MatterCCZ4RHS<WallRadiation, CosmoMovingPunctureGauge> evolution(
        matter(), m_p.ccz4_params, m_dx, m_p.sigma, m_p.formulation, m_p.background.G);
    evolution.set_gauge(gauge);
    BoxLoops::loop(evolution, soln, rhs, EXCLUDE_GHOST_CELLS);
}

void RadiationFLRWLevel::specificAdvance()
{
    BoxLoops::loop(TraceARemoval(), m_state_new, m_state_new, INCLUDE_GHOST_CELLS);
    if (m_p.nan_check)
        BoxLoops::loop(NanCheck(m_dx, m_p.center, "RadiationFLRW advance"),
                      m_state_new, m_state_new, EXCLUDE_GHOST_CELLS, disable_simd());
}

void RadiationFLRWLevel::specificUpdateODE(GRLevelData &soln,
                                          const GRLevelData &, Real)
{
    BoxLoops::loop(TraceARemoval(), soln, soln, INCLUDE_GHOST_CELLS);
}

void RadiationFLRWLevel::computeTaggingCriterion(FArrayBox &criterion,
                                                 const FArrayBox &,
                                                 const FArrayBox &)
{
    criterion.setVal(0.); // Pass A enforces max_level=0.
}

void RadiationFLRWLevel::specificPostTimeStep()
{
    refresh_diagnostics();
    AMRReductions<VariableType::diagnostic> diag(m_gr_amr);
    AMRReductions<VariableType::evolution> evo(m_gr_amr);
    const double coordinate_volume = evo.get_domain_volume();
    const double proper_volume = diag.sum(c_sqrt_gamma);
    const auto bg = RadiationFLRWBackground(m_p.background).at(m_time);
    // Both averages use a consistent measure. In homogeneous FLRW both give
    // chi=a^-2, unlike sum(chi)/proper_volume which incorrectly gives a^-5.
    const double chi_coordinate = evo.sum(c_chi) / coordinate_volume;
    const double chi_proper = diag.sum(c_chi_proper) / proper_volume;
    const std::string file = m_p.data_path + "flrw";
    const bool first = (m_time == 0.) || !std::ifstream(file + ".dat").good();
    SmallDataIO output(file, m_dt, m_time, m_restart_time, SmallDataIO::APPEND,
                       first, ".dat", 16, 10);
    if (!first)
        output.remove_duplicate_time_data();
    if (first)
        output.write_header_line({"a_bg", "chi_mean", "chi_proper_mean", "Pi_rad_mean",
            "psi_rad_mean", "rho_rad_mean", "rho_wall_mean", "K_mean", "lapse_mean",
            "Ham_rms", "Mom_rms", "Ham_max", "Mom_max", "S_minus_E_max", "current",
            "timelike_min", "proper_volume", "Theta_max", "chi_min", "chi_max", "lapse_min"});
    output.write_time_data_line({bg.a, chi_coordinate, chi_proper,
        evo.sum(c_Pi_rad) / coordinate_volume, evo.sum(c_psi_rad) / coordinate_volume,
        diag.sum(c_rho_rad) / coordinate_volume, diag.sum(c_rho_wall) / coordinate_volume,
        evo.sum(c_K) / coordinate_volume, evo.sum(c_lapse) / coordinate_volume,
        diag.norm(c_Ham, 2, true), diag.norm(c_Mom, 2, true),
        diag.norm(c_Ham, 0), diag.norm(c_Mom, 0), diag.norm(c_S_minus_E_rad, 0),
        diag.sum(c_current_proper), diag.min(c_timelike_margin), proper_volume,
        evo.norm(c_Theta, 0), evo.min(c_chi), evo.max(c_chi), evo.min(c_lapse)});
}

#ifdef CH_USE_HDF5
void RadiationFLRWLevel::writeCheckpointHeader(HDF5Handle &handle) const
{
    GRAMRLevel::writeCheckpointHeader(handle);
    HDF5HeaderData header;
    header.m_int["radiation_flrw_schema"] = 1;
    header.m_int["radiation_enabled"] = m_p.radiation_enabled;
    for (const auto &entry : m_p.checkpoint_metadata())
        header.m_real[entry.first] = entry.second;
    header.writeToFile(handle);
}

void RadiationFLRWLevel::readCheckpointHeader(HDF5Handle &handle)
{
    GRAMRLevel::readCheckpointHeader(handle);
    HDF5HeaderData header;
    header.readFromFile(handle);
    if (header.m_int.count("radiation_flrw_schema") == 0 ||
        header.m_int.at("radiation_flrw_schema") != 1 ||
        header.m_int.count("radiation_enabled") == 0 ||
        header.m_int.at("radiation_enabled") != int(m_p.radiation_enabled))
        MayDay::Error("RadiationFLRW checkpoint schema/enabled flag missing or incompatible");
    for (const auto &entry : m_p.checkpoint_metadata())
    {
        if (header.m_real.count(entry.first) == 0 ||
            header.m_real.at(entry.first) != entry.second)
            MayDay::Error(("RadiationFLRW restart normalization mismatch: " + entry.first).c_str());
    }
}
#endif
