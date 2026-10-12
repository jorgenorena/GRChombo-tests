#ifndef SCALARRADIATION_HPP_
#define SCALARRADIATION_HPP_

#include "RadiationField.hpp"
#include "ScalarField.hpp"

// The canonical wall and radiation exchange energy only through gravity.
template <class potential_t> class ScalarRadiation
{
    ScalarField<potential_t> m_wall;
    RadiationField m_radiation;

  public:
    ScalarRadiation(const potential_t &potential, double kappa,
                    bool radiation_enabled = true)
        : m_wall(potential), m_radiation(kappa, radiation_enabled) {}

    template <class data_t>
    struct Vars : ScalarField<potential_t>::template Vars<data_t>,
                  RadiationField::Vars<data_t>
    {
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t f)
        {
            ScalarField<potential_t>::template Vars<data_t>::enum_mapping(f);
            RadiationField::Vars<data_t>::enum_mapping(f);
        }
    };
    template <class data_t>
    struct Diff2Vars : ScalarField<potential_t>::template Diff2Vars<data_t>,
                      RadiationField::Diff2Vars<data_t>
    {
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t f)
        {
            ScalarField<potential_t>::template Diff2Vars<data_t>::enum_mapping(f);
            RadiationField::Diff2Vars<data_t>::enum_mapping(f);
        }
    };

    const ScalarField<potential_t> &wall() const { return m_wall; }
    const RadiationField &radiation() const { return m_radiation; }

    template <class data_t, template <typename> class vars_t>
    emtensor_t<data_t> compute_emtensor(
        const vars_t<data_t> &vars, const vars_t<Tensor<1, data_t>> &d1,
        const Tensor<2, data_t> &h_UU, const Tensor<3, data_t> &chris) const
    {
        auto out = m_wall.compute_emtensor(vars, d1, h_UU, chris);
        const auto rad = m_radiation.compute_emtensor(vars, d1, h_UU, chris);
        out.rho += rad.rho;
        out.S += rad.S;
        FOR(i) { out.Si[i] += rad.Si[i]; }
        FOR(i, j) { out.Sij[i][j] += rad.Sij[i][j]; }
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
        m_wall.add_matter_rhs(rhs, vars, d1, d2, advec);
        m_radiation.add_matter_rhs(rhs, vars, d1, d2, advec);
    }
};
#endif
