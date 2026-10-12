#include "RadiationHandoffLevel.hpp"
#include "AMRReductions.hpp"
#include "BoxLoops.hpp"
#include "ComputePack.hpp"
#include "NanCheck.hpp"
#include "NewMatterConstraints.hpp"
#include "SmallDataIO.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace
{
bool close(double a, double b)
{
    return std::isfinite(a) && std::isfinite(b) &&
           std::abs(a-b) <= 64.*std::numeric_limits<double>::epsilon()*
                            std::max({std::abs(a), std::abs(b), 1.e-300});
}
bool provenance_key(const std::string &s)
{
    return s.find("radiation_") == 0 || s.find("wall_") == 0 ||
           s.find("initial_") == 0 || s == "G_Newton";
}
class InitialMatterCheck
{
    double m_Pi;
  public:
    explicit InitialMatterCheck(double Pi) : m_Pi(Pi) {}
    void compute(Cell<double> cell) const
    {
        if (!close(cell.load_vars(c_Pi_rad), m_Pi) || cell.load_vars(c_psi_rad) != 0. ||
            cell.load_vars(c_Pi) != 0.)
            MayDay::Error("RadiationHandoff fields do not match the initial radiation normalization");
    }
};
class TotalDiagnostics
{
    WallRadiation m_matter;
    FourthOrderDerivatives m_deriv;
    double m_G;
  public:
    TotalDiagnostics(WallRadiation matter, double dx, double G)
        : m_matter(matter), m_deriv(dx), m_G(G) {}
    template <class T> void compute(Cell<T> cell) const
    {
        const auto v = cell.template load_vars<FLRWDiagnostics::Vars>();
        const auto d = m_deriv.template diff1<FLRWDiagnostics::Vars>(cell);
        const auto inv = TensorAlgebra::compute_inverse_sym(v.h);
        const auto chris = TensorAlgebra::compute_christoffel(d.h, inv);
        const auto em = m_matter.compute_emtensor(v, d, inv, chris.ULL);
        cell.store_vars(em.rho, c_rho_total);
        cell.store_vars(v.K*v.K - 24.*M_PI*m_G*em.rho, c_hybrid_defect);
    }
};
}

void RadiationHandoffLevel::refresh_diagnostics()
{
    // This zero-step inspector retains the solver's cubic physical ghosts.
    // A separate control can apply the evolution code's linear refill to
    // measure boundary incompatibility before any time evolution is attempted.
    if (m_p.handoff_refill_boundaries)
        fillAllGhosts();
    else
        m_state_new.exchange(Interval(0, NUM_VARS-1), m_exchange_copier);
    BoxLoops::loop(make_compute_pack(
        NanCheck(m_dx, m_p.center, "RadiationHandoff import"),
        RadiationStateCheck(m_dx, true),
        InitialMatterCheck(RadiationFLRWBackground(m_p.background).Pi0())),
        m_state_new, m_state_new, EXCLUDE_GHOST_CELLS, disable_simd());
    const auto model = matter();
    BoxLoops::loop(make_compute_pack(
        MatterConstraints<WallRadiation>(model, m_dx, m_p.background.G, c_Ham,
            Interval(c_Mom, c_Mom), c_Ham_abs, Interval(c_Mom_abs, c_Mom_abs)),
        MatterConstraints<WallRadiation>(model, m_dx, m_p.background.G, -1, Interval(c_Mom1, c_Mom3)),
        MatterConstraints<ScalarField<Potential>>(model.wall(), m_dx, m_p.background.G, c_Ham_wall, Interval()),
        FLRWDiagnostics(model, m_dx, m_p.background.kappa, true),
        TotalDiagnostics(model, m_dx, m_p.background.G)),
        m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
}
void RadiationHandoffLevel::postRestart()
{
    if (m_time != 0.) MayDay::Error("Pass D requires an initial slice at t=0");
    refresh_diagnostics();
    pout() << "Pass D reconstructed wall+radiation sources and constraints at t=0, without changing valid fields." << std::endl;
}
void RadiationHandoffLevel::specificPostTimeStep()
{
    refresh_diagnostics();
    AMRReductions<VariableType::diagnostic> d(m_gr_amr);
    AMRReductions<VariableType::evolution> e(m_gr_amr);
    const double volume = e.get_domain_volume();
    // sum() includes reflected copies; domain volume and normalized norms do
    // not. Report integrals over the actual computational domain consistently.
    const double symmetry = get_symmetry_factor();
    const auto integral = [&](int var) { return d.sum(var)/symmetry; };
    const auto maximum_abs = [&](int var) { return std::max(std::abs(d.min(var)), std::abs(d.max(var))); };
    SmallDataIO out(m_p.data_path + "handoff", m_dt, m_time, m_restart_time,
                    SmallDataIO::NEW, true, ".dat", 16, 12);
    out.write_header_line({"rho_rad_mean", "rho_wall_mean", "rho_total_mean", "K_min", "K_max",
        "chi_min", "chi_max", "Ham_rms", "Ham_max", "Mom_rms", "Mom_max", "Ham_terms_rms",
        "Mom_terms_rms", "trace_max", "timelike_min", "current", "proper_volume"});
    out.write_time_data_line({integral(c_rho_rad)/volume, integral(c_rho_wall)/volume, integral(c_rho_total)/volume,
        e.min(c_K), e.max(c_K), e.min(c_chi), e.max(c_chi), d.norm(c_Ham,2,true), maximum_abs(c_Ham),
        d.norm(c_Mom,2,true), maximum_abs(c_Mom), d.norm(c_Ham_abs,2,true), d.norm(c_Mom_abs,2,true),
        maximum_abs(c_S_minus_E_rad), d.min(c_timelike_margin), integral(c_current_proper), integral(c_sqrt_gamma)});
}
void RadiationHandoffLevel::readCheckpointHeader(HDF5Handle &handle)
{
    GRAMRLevel::readCheckpointHeader(handle);
    HDF5HeaderData h; h.readFromFile(handle);
    auto fail = [](const std::string &name) { MayDay::Error(("RadiationHandoff metadata mismatch: " + name).c_str()); };
    if (h.m_int["radiation_handoff_schema"] != 1 || h.m_int["radiation_enabled"] != 1 ||
        h.m_int.count("initial_wall_enabled") == 0 || h.m_int["initial_wall_enabled"] != int(m_p.wall_enabled)) fail("schema/wall flag");
    if (h.m_string["initial_data_method"] != "CTTKHybrid-local-K-expanding" ||
        h.m_string["initial_data_source"] != "GRTresna-RadiationID") fail("initial-data prescription");
    for (const auto &entry : m_p.checkpoint_metadata())
        if (!h.m_real.count(entry.first) || !close(h.m_real[entry.first],entry.second)) fail(entry.first);
    for (const auto &key : {"initial_Ham_error_percent", "initial_Mom_error_percent"})
        if (!h.m_real.count(key) || !std::isfinite(h.m_real[key]) || h.m_real[key] < 0. ||
            h.m_real[key] > m_p.initial_data_residual_tolerance) fail(key);
    if (!h.m_real.count("initial_residual_tolerance_percent") ||
        !std::isfinite(h.m_real["initial_residual_tolerance_percent"]) ||
        h.m_real["initial_residual_tolerance_percent"] <= 0.) fail("export tolerance");
    for (int d = 0; d < CH_SPACEDIM; ++d)
    {
        const auto suffix = std::to_string(d);
        const int n = m_p.ivN[d]+1;
        for (const auto &key : {"initial_N_", "initial_periodic_", "initial_lo_boundary_", "initial_hi_boundary_"})
            if (!h.m_int.count(key+suffix)) fail(key+suffix);
        for (const auto &key : {"initial_L_", "initial_center_"})
            if (!h.m_real.count(key+suffix)) fail(key+suffix);
        if (h.m_int["initial_N_"+suffix] != n ||
            h.m_int["initial_periodic_"+suffix] != int(m_p.boundary_params.is_periodic[d]) ||
            !close(h.m_real["initial_L_"+suffix],n*m_p.coarsest_dx) ||
            !close(h.m_real["initial_center_"+suffix],m_p.center[d])) fail("domain "+suffix);
        if (!m_p.boundary_params.is_periodic[d])
        {
            const int lo = h.m_int["initial_lo_boundary_"+suffix] == 1 ? BoundaryConditions::REFLECTIVE_BC : BoundaryConditions::EXTRAPOLATING_BC;
            const int hi = h.m_int["initial_hi_boundary_"+suffix] == 1 ? BoundaryConditions::REFLECTIVE_BC : BoundaryConditions::EXTRAPOLATING_BC;
            if (lo != m_p.boundary_params.lo_boundary[d] || hi != m_p.boundary_params.hi_boundary[d]) fail("boundaries "+suffix);
        }
    }
    for (const auto &v : h.m_real) if (provenance_key(v.first)) m_initial_header.m_real[v.first]=v.second;
    for (const auto &v : h.m_int) if (provenance_key(v.first)) m_initial_header.m_int[v.first]=v.second;
    for (const auto &v : h.m_string) if (provenance_key(v.first)) m_initial_header.m_string[v.first]=v.second;
}
void RadiationHandoffLevel::writeCheckpointHeader(HDF5Handle &handle) const
{
    GRAMRLevel::writeCheckpointHeader(handle);
    m_initial_header.writeToFile(handle);
    HDF5HeaderData inspection;
    inspection.m_int["handoff_refilled_boundaries"] = m_p.handoff_refill_boundaries;
    inspection.writeToFile(handle);
}
