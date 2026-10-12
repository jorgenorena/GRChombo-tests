#ifndef RADIATIONGRADIENTTAGGING_HPP_
#define RADIATIONGRADIENTTAGGING_HPP_
#include "AcousticMatter.hpp"

// Dimensionless radiation-gradient indicator. It is invariant under field
// normalization changes and vanishes for homogeneous radiation. No density
// background is tagged merely because it is nonzero.
class RadiationGradientTagging
{
    FourthOrderDerivatives m_deriv;
    double m_Pi0, m_k;
  public:
    RadiationGradientTagging(double dx, double Pi0, double k)
        : m_deriv(dx), m_Pi0(Pi0), m_k(k) {}
    template <class T> void compute(Cell<T> cell) const
    {
        const auto vars = cell.template load_vars<FLRWDiagnostics::Vars>();
        const auto d1 = m_deriv.template diff1<RadiationField::Vars>(cell);
        const auto inv = TensorAlgebra::compute_inverse_sym(vars.h);
        T norm2 = 0.;
        FOR(i,j) { norm2 += vars.chi * inv[i][j] * d1.Pi_rad[i] * d1.Pi_rad[j]; }
        cell.store_vars(sqrt(simd_max(norm2, T(0.))) / (m_k * m_Pi0), c_radiation_tag);
    }
};
#endif
