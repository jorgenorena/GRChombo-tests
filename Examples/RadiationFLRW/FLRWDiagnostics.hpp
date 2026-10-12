#ifndef FLRWDIAGNOSTICS_HPP_
#define FLRWDIAGNOSTICS_HPP_
#include "MatterCCZ4RHS.hpp"
#include "Potential.hpp"
#include "ScalarRadiation.hpp"

using WallRadiation = ScalarRadiation<Potential>;

class FLRWDiagnostics
{
    WallRadiation m_matter;
    FourthOrderDerivatives m_deriv;
    double m_kappa;
    bool m_enabled;

  public:
    template <class T> using Vars = MatterCCZ4RHS<WallRadiation>::Vars<T>;
    FLRWDiagnostics(WallRadiation matter, double dx, double kappa, bool enabled)
        : m_matter(matter), m_deriv(dx), m_kappa(kappa), m_enabled(enabled) {}

    template <class data_t> void compute(Cell<data_t> cell) const
    {
        const auto vars = cell.template load_vars<Vars>();
        const auto d1 = m_deriv.template diff1<Vars>(cell);
        const auto h_UU = TensorAlgebra::compute_inverse_sym(vars.h);
        const auto chris = TensorAlgebra::compute_christoffel(d1.h, h_UU);
        const auto rad = m_matter.radiation().compute_emtensor(vars, d1, h_UU, chris.ULL);
        const auto wall = m_matter.wall().compute_emtensor(vars, d1, h_UU, chris.ULL);
        const data_t volume = sqrt(TensorAlgebra::compute_determinant_sym(vars.h)) *
                              pow(vars.chi, -1.5);
        data_t margin = 1., J = 0.;
        if (m_enabled)
        {
            data_t w2 = 0.;
            FOR(i, j)
            {
                w2 += vars.chi * h_UU[i][j] * (d1.psi_rad[i] / vars.Pi_rad) *
                      (d1.psi_rad[j] / vars.Pi_rad);
            }
            margin -= w2;
            J = m_kappa * vars.Pi_rad * vars.Pi_rad * vars.Pi_rad * margin;
        }
        cell.store_vars(rad.rho, c_rho_rad);
        cell.store_vars(wall.rho, c_rho_wall);
        cell.store_vars(rad.S - rad.rho, c_S_minus_E_rad);
        cell.store_vars(volume, c_sqrt_gamma);
        cell.store_vars(volume * vars.chi, c_chi_proper);
        cell.store_vars(volume * J, c_current_proper);
        cell.store_vars(margin, c_timelike_margin);
    }
};

// Fail explicitly before dividing by Pi_rad. This check changes no fields.
// Use with disable_simd() at every RK stage, initial slice and restart.
class RadiationStateCheck
{
    FourthOrderDerivatives m_deriv;
    bool m_enabled;
  public:
    RadiationStateCheck(double dx, bool enabled) : m_deriv(dx), m_enabled(enabled) {}
    void compute(Cell<double> cell) const
    {
        if (!m_enabled)
            return;
        const auto vars = cell.load_vars<FLRWDiagnostics::Vars>();
        const auto d1 = m_deriv.diff1<FLRWDiagnostics::Vars>(cell);
        if (!(std::isfinite(vars.Pi_rad) && vars.Pi_rad != 0. &&
              std::isfinite(vars.chi) && vars.chi > 0.))
            MayDay::Error("Invalid radiation state: Pi_rad zero/nonfinite or chi nonpositive");
        const auto h_UU = TensorAlgebra::compute_inverse_sym(vars.h);
        double w2 = 0.;
        FOR(i, j)
        {
            w2 += vars.chi * h_UU[i][j] * (d1.psi_rad[i] / vars.Pi_rad) *
                  (d1.psi_rad[j] / vars.Pi_rad);
        }
        if (!(std::isfinite(w2) && w2 >= 0. && w2 < 1.))
            MayDay::Error("Invalid radiation state: timelike margin is nonpositive/nonfinite");
    }
};
#endif
