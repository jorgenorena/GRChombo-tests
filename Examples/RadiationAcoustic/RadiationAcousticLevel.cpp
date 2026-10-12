#include "RadiationAcousticLevel.hpp"
#include "AMRReductions.hpp"
#include "BoxLoops.hpp"
#include "ComputePack.hpp"
#include "NanCheck.hpp"
#include "NewMatterConstraints.hpp"
#include "SetValue.hpp"
#include "SmallDataIO.hpp"
#include <fstream>

namespace
{
class InitialAcoustic
{
    AcousticGeometry m_geometry;
    double m_dx, m_Pi0, m_amplitude, m_phi;
  public:
    InitialAcoustic(const SimulationParameters &p, double dx)
        : m_geometry(p.geometry), m_dx(dx), m_Pi0(p.Pi0),
          m_amplitude(p.wave_amplitude), m_phi(p.phi_true) {}
    template <class T> void compute(Cell<T> cell) const
    {
        const Coordinates<T> coords(cell.get_int_vect(), m_dx);
        const T f = m_geometry.jacobian(coords.x);
        const T fp = m_geometry.jacobian_derivative(coords.x);
        const T phase = m_geometry.k * m_geometry.X(coords.x);
        const T chi = pow(f, -2. / 3.);
        cell.store_vars(chi, c_chi);
        cell.store_vars(chi * f * f, c_h11);
        cell.store_vars(chi, c_h22);
        cell.store_vars(chi, c_h33);
        cell.store_vars((4. / 3.) * pow(f, -7. / 3.) * fp, c_Gamma1);
        cell.store_vars(T(1.), c_lapse);
        cell.store_vars(T(m_phi), c_phi);
        cell.store_vars(m_amplitude * m_Pi0 / m_geometry.k * sin(phase), c_psi_rad);
        cell.store_vars(m_Pi0 * (1. - m_amplitude / std::sqrt(3.) * cos(phase)), c_Pi_rad);
    }
};
}

void RadiationAcousticLevel::initialData()
{
    BoxLoops::loop(make_compute_pack(SetValue(0.), InitialAcoustic(m_p, m_dx)),
                  m_state_new, m_state_new, INCLUDE_GHOST_CELLS);
    refresh_diagnostics();
}

void RadiationAcousticLevel::refresh_diagnostics()
{
    fillAllGhosts();
    BoxLoops::loop(RadiationStateCheck(m_dx, true), m_state_new, m_state_new,
                  EXCLUDE_GHOST_CELLS, disable_simd());
    const auto model = matter();
    BoxLoops::loop(make_compute_pack(
        MatterConstraints<WallRadiation>(model, m_dx, m_p.G_Newton, c_Ham, Interval(c_Mom, c_Mom)),
        FLRWDiagnostics(model, m_dx, m_p.radiation_kappa, true),
        AcousticDiagnostics(model, m_dx, m_p.geometry, m_p.Pi0, m_p.radiation_kappa, m_time)),
        m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
}

void RadiationAcousticLevel::specificEvalRHS(GRLevelData &soln, GRLevelData &rhs,
                                            double)
{
    BoxLoops::loop(RadiationStateCheck(m_dx, true), soln, soln,
                  EXCLUDE_GHOST_CELLS, disable_simd());
    // Prescribed stationary flat geometry isolates the matter PDE. The full
    // composite matter RHS is unchanged; neither field is reset to a wave.
    BoxLoops::loop(make_compute_pack(SetValue(0.),
        FixedMetricMatterRHS(matter(), m_dx, m_p.sigma)), soln, rhs, EXCLUDE_GHOST_CELLS);
}

void RadiationAcousticLevel::specificAdvance()
{
    if (m_p.nan_check)
        BoxLoops::loop(NanCheck(m_dx, m_p.center, "RadiationAcoustic advance"),
                      m_state_new, m_state_new, EXCLUDE_GHOST_CELLS, disable_simd());
}

void RadiationAcousticLevel::computeTaggingCriterion(FArrayBox &criterion,
                                                     const FArrayBox &,
                                                     const FArrayBox &)
{
    criterion.setVal(0.); // max_level=0 is required for Pass B.
}

void RadiationAcousticLevel::specificPostTimeStep()
{
    refresh_diagnostics();
    AMRReductions<VariableType::diagnostic> diag(m_gr_amr);
    AMRReductions<VariableType::evolution> evo(m_gr_amr);
    const double volume = diag.sum(c_sqrt_gamma);
    const std::string file = m_p.data_path + "acoustic";
    const bool first = m_time == 0. || !std::ifstream(file + ".dat").good();
    SmallDataIO output(file, m_dt, m_time, m_restart_time, SmallDataIO::APPEND,
                       first, ".dat", 16, 10);
    if (!first)
        output.remove_duplicate_time_data();
    if (first)
        output.write_header_line({"current", "energy", "proper_volume", "timelike_min",
            "rho_min", "rho_max", "trace_max", "Pi_min", "Pi_max", "grad_Pi_max",
            "grad_velocity_max", "Pi_cos", "Pi_sin", "psi_cos", "psi_sin", "Ham_max", "Mom_max"});
    output.write_time_data_line({diag.sum(c_current_proper), diag.sum(c_energy_proper), volume,
        diag.min(c_timelike_margin), diag.min(c_rho_rad), diag.max(c_rho_rad),
        diag.norm(c_S_minus_E_rad, 0), evo.min(c_Pi_rad), evo.max(c_Pi_rad),
        diag.max(c_grad_Pi), diag.max(c_grad_velocity),
        2. * diag.sum(c_Pi_cos) / volume, 2. * diag.sum(c_Pi_sin) / volume,
        2. * diag.sum(c_psi_cos) / volume, 2. * diag.sum(c_psi_sin) / volume,
        diag.norm(c_Ham, 0), diag.norm(c_Mom, 0)});
}

#ifdef CH_USE_HDF5
void RadiationAcousticLevel::writeCheckpointHeader(HDF5Handle &handle) const
{
    GRAMRLevel::writeCheckpointHeader(handle);
    HDF5HeaderData header;
    header.m_int["radiation_acoustic_schema"] = 1;
    for (const auto &entry : m_p.checkpoint_metadata())
        header.m_real[entry.first] = entry.second;
    header.writeToFile(handle);
}

void RadiationAcousticLevel::readCheckpointHeader(HDF5Handle &handle)
{
    GRAMRLevel::readCheckpointHeader(handle);
    HDF5HeaderData header;
    header.readFromFile(handle);
    if (header.m_int.count("radiation_acoustic_schema") == 0 ||
        header.m_int.at("radiation_acoustic_schema") != 1)
        MayDay::Error("RadiationAcoustic checkpoint schema missing or incompatible");
    for (const auto &entry : m_p.checkpoint_metadata())
        if (header.m_real.count(entry.first) == 0 || header.m_real.at(entry.first) != entry.second)
            MayDay::Error(("RadiationAcoustic restart mismatch: " + entry.first).c_str());
}
#endif
