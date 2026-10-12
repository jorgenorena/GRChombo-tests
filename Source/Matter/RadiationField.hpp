/* GRChombo: irrotational radiation, P(X) = radiation_kappa * X^2.
 * X = -(grad psi_rad)^2/2, Pi_rad = n^mu partial_mu psi_rad.
 */
#ifndef RADIATIONFIELD_HPP_
#define RADIATIONFIELD_HPP_

#include "CCZ4Geometry.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp"
#include "VarsTools.hpp"

class RadiationField
{
    double m_kappa;
    bool m_enabled;

  public:
    RadiationField(double a_kappa, bool a_enabled = true)
        : m_kappa(a_kappa), m_enabled(a_enabled) {}

    template <class data_t> struct Vars
    {
        data_t psi_rad, Pi_rad;
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t f)
        {
            VarsTools::define_enum_mapping(f, c_psi_rad, psi_rad);
            VarsTools::define_enum_mapping(f, c_Pi_rad, Pi_rad);
        }
    };

    template <class data_t> struct Diff2Vars
    {
        data_t psi_rad;
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t f)
        {
            VarsTools::define_enum_mapping(f, c_psi_rad, psi_rad);
        }
    };

    template <class data_t, template <typename> class vars_t>
    emtensor_t<data_t> compute_emtensor(
        const vars_t<data_t> &vars, const vars_t<Tensor<1, data_t>> &d1,
        const Tensor<2, data_t> &h_UU,
        const Tensor<3, data_t> & /*conformal_chris*/) const
    {
        emtensor_t<data_t> out;
        out.rho = 0.;
        out.S = 0.;
        out.Si = 0.;
        out.Sij = 0.;
        if (!m_enabled)
            return out;

        data_t s = 0.;
        FOR(i, j)
        {
            s += vars.chi * h_UU[i][j] * d1.psi_rad[i] * d1.psi_rad[j];
        }
        const data_t X = 0.5 * (vars.Pi_rad * vars.Pi_rad - s);
        const data_t P = m_kappa * X * X;
        const data_t PX = 2. * m_kappa * X;
        out.rho = PX * vars.Pi_rad * vars.Pi_rad - P;
        FOR(i) { out.Si[i] = -PX * vars.Pi_rad * d1.psi_rad[i]; }
        FOR(i, j)
        {
            out.Sij[i][j] = PX * d1.psi_rad[i] * d1.psi_rad[j] +
                            P * vars.h[i][j] / vars.chi;
        }
        // Compute the trace independently; S-rho is a useful diagnostic.
        out.S = vars.chi * TensorAlgebra::compute_trace(out.Sij, h_UU);
        return out;
    }

    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t,
              template <typename> class rhs_vars_t>
    void add_matter_rhs(rhs_vars_t<data_t> &rhs, const vars_t<data_t> &vars,
                        const vars_t<Tensor<1, data_t>> &d1,
                        const diff2_vars_t<Tensor<2, data_t>> &d2,
                        const vars_t<data_t> &advec) const
    {
        // Exactly zero radiation is degenerate in the X^2 equation. The
        // explicitly disabled mode bypasses it, without a floor on X or Pi.
        rhs.psi_rad = 0.;
        rhs.Pi_rad = 0.;
        if (!m_enabled)
            return;

        const auto h_UU = TensorAlgebra::compute_inverse_sym(vars.h);
        const auto chris = TensorAlgebra::compute_christoffel(d1.h, h_UU);
        Tensor<1, data_t> w, w_U = 0., dchi_U = 0.;
        FOR(i) { w[i] = d1.psi_rad[i] / vars.Pi_rad; }
        FOR(i, j)
        {
            w_U[i] += vars.chi * h_UU[i][j] * w[j];
            dchi_U[i] += h_UU[i][j] * d1.chi[j];
        }
        data_t w2 = 0., laplacian = 0., ww_hessian = 0.;
        data_t w_dPi = 0., K_ww = 0., Phi_dalpha = 0.;
        FOR(i)
        {
            w2 += w[i] * w_U[i];
            w_dPi += w_U[i] * d1.Pi_rad[i];
            Phi_dalpha += vars.Pi_rad * w_U[i] * d1.lapse[i];
        }
        FOR(i, j)
        {
            data_t hessian = d2.psi_rad[i][j];
            FOR(k)
            {
                // Physical gamma_ij = h_ij/chi, not just the conformal
                // connection. The Kronecker terms include both chi gradients.
                const data_t physical_chris = chris.ULL[k][i][j] -
                    ((k == i ? d1.chi[j] : data_t(0.)) +
                     (k == j ? d1.chi[i] : data_t(0.)) -
                     vars.h[i][j] * dchi_U[k]) / (2. * vars.chi);
                hessian -= physical_chris * d1.psi_rad[k];
            }
            laplacian += vars.chi * h_UU[i][j] * hessian;
            ww_hessian += w_U[i] * w_U[j] * hessian;
            const data_t Kij =
                (vars.A[i][j] + vars.h[i][j] * vars.K / 3.) / vars.chi;
            K_ww += Kij * w_U[i] * w_U[j];
        }
        rhs.psi_rad = vars.lapse * vars.Pi_rad + advec.psi_rad;
        // Divide the full equation by Pi_rad^2 analytically. On the timelike
        // branch the denominator 3-w^2 is in (2,3], even as Pi redshifts.
        rhs.Pi_rad = advec.Pi_rad + Phi_dalpha + vars.lapse / (3. - w2) *
            ((1. - w2) * (laplacian + vars.K * vars.Pi_rad) -
             2. * ww_hessian + 4. * w_dPi + 2. * vars.Pi_rad * K_ww);
    }
};
#endif
