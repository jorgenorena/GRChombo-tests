#ifndef ACOUSTICMATTER_HPP_
#define ACOUSTICMATTER_HPP_
#include "AcousticGeometry.hpp"
#include "Coordinates.hpp"
// Reuse the composite matter diagnostics and the timelike guard from Pass A.
#include "../RadiationFLRW/FLRWDiagnostics.hpp"

class FixedMetricMatterRHS
{
    WallRadiation m_matter;
    FourthOrderDerivatives m_deriv;
    double m_sigma;
  public:
    template <class T> using Vars = FLRWDiagnostics::Vars<T>;
    FixedMetricMatterRHS(WallRadiation matter, double dx, double sigma)
        : m_matter(matter), m_deriv(dx), m_sigma(sigma) {}
    template <class T> void compute(Cell<T> cell) const
    {
        const auto vars = cell.template load_vars<Vars>();
        const auto d1 = m_deriv.template diff1<Vars>(cell);
        const auto d2 = m_deriv.template diff2<WallRadiation::Diff2Vars>(cell);
        Vars<T> advec;
        advec.enum_mapping([](int, T &value) { value = 0.; }); // beta=0
        WallRadiation::Vars<T> rhs;
        m_matter.add_matter_rhs(rhs, vars, d1, d2, advec);
        // Dissipation acts only on evolving matter. The fixed metric is not
        // dissipated, even when its components vary in remapped coordinates.
        m_deriv.add_dissipation(rhs, cell, m_sigma);
        cell.store_vars(rhs);
    }
};

class AcousticDiagnostics
{
    WallRadiation m_matter;
    FourthOrderDerivatives m_deriv;
    AcousticGeometry m_geometry;
    double m_Pi0, m_kappa, m_time;
  public:
    AcousticDiagnostics(WallRadiation matter, double dx, AcousticGeometry geometry,
                        double Pi0, double kappa, double time)
        : m_matter(matter), m_deriv(dx), m_geometry(geometry),
          m_Pi0(Pi0), m_kappa(kappa), m_time(time) {}
    template <class T> void compute(Cell<T> cell) const
    {
        const auto vars = cell.template load_vars<FLRWDiagnostics::Vars>();
        const auto d1 = m_deriv.template diff1<FLRWDiagnostics::Vars>(cell);
        const auto d2 = m_deriv.template diff2<RadiationField::Diff2Vars>(cell);
        const auto h_UU = TensorAlgebra::compute_inverse_sym(vars.h);
        const auto chris = TensorAlgebra::compute_christoffel(d1.h, h_UU);
        const auto rad = m_matter.radiation().compute_emtensor(vars, d1, h_UU, chris.ULL);
        const T volume = sqrt(TensorAlgebra::compute_determinant_sym(vars.h)) *
                         pow(vars.chi, -1.5);
        const Coordinates<T> coords(cell.get_int_vect(), m_deriv.m_dx);
        const T f = m_geometry.jacobian(coords.x);
        const T fp = m_geometry.jacobian_derivative(coords.x);
        const T phase = m_geometry.k * m_geometry.X(coords.x);
        T s = 0.;
        FOR(i, j) { s += vars.chi * h_UU[i][j] * d1.psi_rad[i] * d1.psi_rad[j]; }
        const T X = 0.5 * (vars.Pi_rad * vars.Pi_rad - s);
        cell.store_vars(3. * m_kappa * X * X, c_rho_rest);
        cell.store_vars(volume * rad.rho, c_energy_proper);
        cell.store_vars(rad.Si[0], c_Srad1);
        cell.store_vars(rad.Sij[0][0], c_Srad11);
        cell.store_vars(rad.Sij[1][1], c_Srad22);
        cell.store_vars(rad.Sij[2][2], c_Srad33);
        // Physical-volume Fourier projections in X. Subtract the known
        // homogeneous part only for diagnostics; the evolved field is intact.
        cell.store_vars(volume * (vars.Pi_rad - m_Pi0) * cos(phase), c_Pi_cos);
        cell.store_vars(volume * (vars.Pi_rad - m_Pi0) * sin(phase), c_Pi_sin);
        cell.store_vars(volume * (vars.psi_rad - m_Pi0 * m_time) * cos(phase), c_psi_cos);
        cell.store_vars(volume * (vars.psi_rad - m_Pi0 * m_time) * sin(phase), c_psi_sin);
        // Avoid the generic SIMD abs wrapper's unqualified scalar overload,
        // which can truncate to int depending on the include environment.
        cell.store_vars(simd_max(d1.Pi_rad[0], -d1.Pi_rad[0]) / (f * m_Pi0), c_grad_Pi);
        // |d v_x(orthonormal)/dX| for this one-dimensional flat-space map.
        const T velocity_derivative = d2.psi_rad[0][0] - fp / f * d1.psi_rad[0] -
                                      d1.psi_rad[0] * d1.Pi_rad[0] / vars.Pi_rad;
        cell.store_vars(simd_max(velocity_derivative, -velocity_derivative) /
                        (f * f * simd_max(vars.Pi_rad, -vars.Pi_rad)), c_grad_velocity);
    }
};
#endif
