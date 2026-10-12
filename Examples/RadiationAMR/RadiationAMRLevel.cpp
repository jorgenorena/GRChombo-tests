#include "RadiationAMRLevel.hpp"
#include "AMRReductions.hpp"
#include "BoxIterator.H"
#include "BoxLoops.hpp"
#include "ComputePack.hpp"
#include "NanCheck.hpp"
#include "NewMatterConstraints.hpp"
#include "RadiationGradientTagging.hpp"
#include "SetValue.hpp"
#include "SmallDataIO.hpp"
#include <fstream>
#include <iomanip>

namespace
{
class InitialWave
{
    double m_dx, m_Pi0, m_amplitude, m_k, m_phi;
    bool m_standing;
  public:
    InitialWave(const SimulationParameters &p, double dx)
        : m_dx(dx), m_Pi0(p.Pi0), m_amplitude(p.wave_amplitude),
          m_k(p.geometry.k), m_phi(p.phi_true), m_standing(p.standing_wave) {}
    template <class T> void compute(Cell<T> cell) const
    {
        const Coordinates<T> coords(cell.get_int_vect(), m_dx);
        const T phase = m_k * coords.x;
        cell.store_vars(T(1.), c_chi);
        cell.store_vars(T(1.), c_h11);
        cell.store_vars(T(1.), c_h22);
        cell.store_vars(T(1.), c_h33);
        cell.store_vars(T(1.), c_lapse);
        cell.store_vars(T(m_phi), c_phi);
        cell.store_vars(m_amplitude * m_Pi0 / m_k *
                        (m_standing ? cos(phase) : sin(phase)), c_psi_rad);
        cell.store_vars(m_Pi0 * (1. - (m_standing ? T(0.) :
                        m_amplitude / std::sqrt(3.) * cos(phase))), c_Pi_rad);
    }
};
}

void RadiationAMRLevel::initialData()
{
    BoxLoops::loop(make_compute_pack(SetValue(0.), InitialWave(m_p, m_dx)),
                  m_state_new, m_state_new, INCLUDE_GHOST_CELLS);
    refresh_diagnostics();
}

void RadiationAMRLevel::check_radiation(GRLevelData &state)
{
    const auto &layout = state.disjointBoxLayout();
    for (DataIterator dit = state.dataIterator(); dit.ok(); ++dit)
    {
        auto &fab = state[dit];
        // A first-derivative stencil reaches two cells: with three ghosts,
        // check valid cells AND the first interpolated ghost layer safely.
        Box checked = grow(layout[dit], 1);
        BoxLoops::loop(RadiationStateCheck(m_dx, true), fab, fab, checked, disable_simd());
        for (BoxIterator bit(fab.box()); bit.ok(); ++bit)
            if (!(std::isfinite(fab(bit(), c_psi_rad)) &&
                  std::isfinite(fab(bit(), c_Pi_rad)) && fab(bit(), c_Pi_rad) != 0.))
                MayDay::Error("Invalid radiation interpolation: nonfinite field or zero Pi_rad");
    }
}

void RadiationAMRLevel::refresh_diagnostics()
{
    if (m_grids.size() == 0) return;
    fillAllGhosts();
    check_radiation(m_state_new);
    const auto model = matter();
    BoxLoops::loop(make_compute_pack(
        MatterConstraints<WallRadiation>(model, m_dx, 0., c_Ham, Interval(c_Mom, c_Mom)),
        FLRWDiagnostics(model, m_dx, m_p.radiation_kappa, true),
        AcousticDiagnostics(model, m_dx, m_p.geometry, m_p.Pi0, m_p.radiation_kappa, m_time),
        RadiationGradientTagging(m_dx, m_p.Pi0, m_p.geometry.k)),
        m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
}

void RadiationAMRLevel::refresh_hierarchy()
{
    for (auto *level : m_gr_amr.get_gramrlevels())
        static_cast<RadiationAMRLevel *>(level)->refresh_diagnostics();
}

void RadiationAMRLevel::specificEvalRHS(GRLevelData &soln, GRLevelData &rhs, double)
{
    // The core has already interpolated the fine ghosts at this RK-stage time.
    check_radiation(soln);
    BoxLoops::loop(make_compute_pack(SetValue(0.), FixedMetricMatterRHS(matter(), m_dx, m_p.sigma)),
                  soln, rhs, EXCLUDE_GHOST_CELLS);
}

void RadiationAMRLevel::specificAdvance()
{
    if (m_p.nan_check)
        BoxLoops::loop(NanCheck(m_dx, m_p.center, "RadiationAMR advance"),
                      m_state_new, m_state_new, EXCLUDE_GHOST_CELLS, disable_simd());
}

void RadiationAMRLevel::computeTaggingCriterion(FArrayBox &criterion,
    const FArrayBox &, const FArrayBox &diagnostics)
{
    for (BoxIterator bit(criterion.box()); bit.ok(); ++bit)
    {
        if (m_p.tagging_mode == 0)
            criterion(bit(), 0) = diagnostics(bit(), c_radiation_tag);
        else
        {
            double x = (bit()[0] + 0.5) * m_dx;
            if (m_p.mirror_tags && x > 0.5 * m_p.wavelength) x = m_p.wavelength - x;
            criterion(bit(), 0) = (x >= m_p.static_refine_lo && x < m_p.static_refine_hi) ?
                                 2. * m_p.regrid_thresholds[m_level] : 0.;
        }
    }
}

std::array<double, 5> RadiationAMRLevel::transfer_state()
{
    AMRReductions<VariableType::diagnostic> d(m_gr_amr);
    double levels = 0., fine_cells = 0.;
    for (auto *base : m_gr_amr.get_gramrlevels())
    {
        const auto *level = static_cast<RadiationAMRLevel *>(base);
        if (level->m_grids.size() == 0) continue;
        ++levels;
        if (level->m_level > 0)
            for (LayoutIterator lit = level->m_grids.layoutIterator(); lit.ok(); ++lit)
                fine_cells += level->m_grids[lit].numPts();
    }
    return {{d.sum(c_current_proper), d.sum(c_energy_proper), d.min(c_timelike_margin), levels, fine_cells}};
}

void RadiationAMRLevel::write_layout()
{
    if (procID() != 0) return;
    std::ofstream out(m_p.data_path + "layouts.dat", std::ios::app);
    out << std::setprecision(17);
    for (auto *base : m_gr_amr.get_gramrlevels())
    {
        const auto *level = static_cast<RadiationAMRLevel *>(base);
        for (LayoutIterator lit = level->m_grids.layoutIterator(); lit.ok(); ++lit)
        {
            const auto &box = level->m_grids[lit];
            out << m_time << " " << level->m_level;
            for (int d = 0; d < CH_SPACEDIM; ++d) out << " " << box.smallEnd(d);
            for (int d = 0; d < CH_SPACEDIM; ++d) out << " " << box.bigEnd(d);
            out << "\n";
        }
    }
}

void RadiationAMRLevel::preRegrid(int base, const Vector<Vector<Box>> &)
{
    if (m_level != base) return;
    refresh_hierarchy();
    m_before_transfer = transfer_state();
}

void RadiationAMRLevel::specificPostRegrid(int base)
{
    if (m_level != base) return;
    refresh_hierarchy(); // guard newly interpolated cells before any evolution
    const auto after = transfer_state();
    const std::string file = m_p.data_path + "transfers";
    const bool first = !std::ifstream(file + ".dat").good();
    SmallDataIO out(file, m_dt, m_time, m_restart_time, SmallDataIO::APPEND, first, ".dat", 16, 12);
    if (first) out.write_header_line({"Q_before", "Q_after", "energy_before", "energy_after",
        "timelike_before", "timelike_after", "levels_before", "levels_after", "fine_cells_before", "fine_cells_after"});
    if (!first) out.remove_duplicate_time_data();
    out.write_time_data_line({m_before_transfer[0], after[0], m_before_transfer[1], after[1],
        m_before_transfer[2], after[2], m_before_transfer[3], after[3], m_before_transfer[4], after[4]});
    write_layout();
}

void RadiationAMRLevel::specificPostTimeStep()
{
    if (m_level != 0) return;
    // Coarse postTimeStep follows fine subcycling and averaging down. Refresh
    // EVERY level here, so reductions never combine asynchronous diagnostics.
    refresh_hierarchy();
    AMRReductions<VariableType::diagnostic> d(m_gr_amr);
    AMRReductions<VariableType::evolution> e(m_gr_amr);
    const auto state = transfer_state();
    const double volume = d.sum(c_sqrt_gamma);
    const std::string file = m_p.data_path + "amr";
    const bool first = !std::ifstream(file + ".dat").good();
    SmallDataIO out(file, m_dt, m_time, m_restart_time, SmallDataIO::APPEND, first, ".dat", 16, 12);
    if (first) out.write_header_line({"current", "energy", "proper_volume", "timelike_min",
        "trace_max", "Pi_min", "Pi_max", "grad_Pi_max", "grad_velocity_max", "tag_max",
        "Pi_cos", "Pi_sin", "psi_cos", "psi_sin", "Ham_max", "Mom_max", "levels", "fine_cells"});
    if (!first) out.remove_duplicate_time_data();
    out.write_time_data_line({state[0], state[1], volume, state[2], d.norm(c_S_minus_E_rad, 0),
        e.min(c_Pi_rad), e.max(c_Pi_rad), d.max(c_grad_Pi), d.max(c_grad_velocity), d.max(c_radiation_tag),
        2*d.sum(c_Pi_cos)/volume, 2*d.sum(c_Pi_sin)/volume, 2*d.sum(c_psi_cos)/volume,
        2*d.sum(c_psi_sin)/volume, d.norm(c_Ham, 0), d.norm(c_Mom, 0), state[3], state[4]});
    if (first) write_layout();
}

#ifdef CH_USE_HDF5
void RadiationAMRLevel::writeCheckpointHeader(HDF5Handle &handle) const
{
    GRAMRLevel::writeCheckpointHeader(handle);
    HDF5HeaderData h;
    h.m_int["radiation_amr_schema"] = 1;
    for (const auto &entry : m_p.checkpoint_metadata()) h.m_real[entry.first] = entry.second;
    h.writeToFile(handle);
}
void RadiationAMRLevel::readCheckpointHeader(HDF5Handle &handle)
{
    GRAMRLevel::readCheckpointHeader(handle);
    HDF5HeaderData h;
    h.readFromFile(handle);
    if (h.m_int.count("radiation_amr_schema") == 0 || h.m_int.at("radiation_amr_schema") != 1)
        MayDay::Error("RadiationAMR checkpoint schema missing or incompatible");
    for (const auto &entry : m_p.checkpoint_metadata())
        if (h.m_real.count(entry.first) == 0 || h.m_real.at(entry.first) != entry.second)
            MayDay::Error(("RadiationAMR restart mismatch: " + entry.first).c_str());
}
#endif
