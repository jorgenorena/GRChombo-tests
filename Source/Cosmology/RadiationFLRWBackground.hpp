#ifndef RADIATIONFLRWBACKGROUND_HPP_
#define RADIATIONFLRWBACKGROUND_HPP_

#include "MayDay.H"
#include <cmath>

// Reference only: never used to overwrite the evolving matter or geometry.
// Coordinate time is proper time in exact FLRW with lapse=1.
class RadiationFLRWBackground
{
  public:
    struct params_t
    {
        double t0 = 0.;
        double rho0 = 0.;
        double kappa = 1.;
        double G = 1.;
    };
    struct state_t
    {
        double a, H, K, chi, rho, Pi, psi;
    };

    explicit RadiationFLRWBackground(params_t p) : m_p(p)
    {
        if (!(std::isfinite(p.t0) && std::isfinite(p.rho0) && p.rho0 >= 0. &&
              std::isfinite(p.kappa) && p.kappa > 0. &&
              std::isfinite(p.G) && p.G > 0.))
            MayDay::Error("Invalid radiation background normalization");
    }
    double H0() const { return std::sqrt(8. * M_PI * m_p.G * m_p.rho0 / 3.); }
    double Pi0() const { return std::pow(4. * m_p.rho0 / (3. * m_p.kappa), 0.25); }
    state_t at(double time) const
    {
        const double age = time - m_p.t0;
        const double a2 = 1. + 2. * H0() * age;
        if (!(std::isfinite(a2) && a2 > 0.))
            MayDay::Error("Radiation FLRW reference is outside its expanding branch");
        const double a = std::sqrt(a2);
        return {a, H0() / a2, -3. * H0() / a2, 1. / a2,
                m_p.rho0 / (a2 * a2), Pi0() / a,
                2. * Pi0() * age / (a + 1.)};
    }

  private:
    params_t m_p;
};
#endif
