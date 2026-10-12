/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef POTENTIAL_HPP_
#define POTENTIAL_HPP_

#include "simd.hpp"
#include "DomainWallPotential.hpp"

class Potential
{
  public:
    struct params_t
    {
        // double scalar_mass;
        double v;
        double lambda;
        double q = 0.0;
    };

  private:
    params_t m_params;
    double m_phi_true;
    double m_phi_false;

  public:
    //! The constructor
    Potential(params_t a_params)
        : m_params(a_params),
          m_phi_true(DomainWallPotential::true_vacuum(a_params.v, a_params.lambda, a_params.q)),
          m_phi_false(DomainWallPotential::false_vacuum(a_params.v, m_phi_true)) {}

    //! Set the potential function for the scalar field here
    template <class data_t, template <typename> class vars_t>
    void compute_potential(data_t &V_of_phi, data_t &dVdphi,
                           const vars_t<data_t> &vars) const
    {
        V_of_phi = DomainWallPotential::value(vars.phi, m_params.v,
                                             m_params.lambda, m_params.q, m_phi_true);
        dVdphi = DomainWallPotential::derivative(vars.phi, m_params.v,
                                                m_params.lambda, m_params.q,
                                                m_phi_true, m_phi_false);
    }
};

#endif /* POTENTIAL_HPP_ */
