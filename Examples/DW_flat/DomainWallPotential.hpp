#ifndef DOMAIN_WALL_POTENTIAL_HPP_
#define DOMAIN_WALL_POTENTIAL_HPP_

#include <cmath>
#include <stdexcept>

// Kept identical in the GRTresna and GRChombo DW_flat examples; the bias
// regression test checks both copies. q has the sign in -q * phi^3.
namespace DomainWallPotential
{
inline double true_vacuum(double v, double lambda, double q)
{
    if (!(std::isfinite(v) && v > 0.0 && std::isfinite(lambda) &&
          lambda > 0.0 && std::isfinite(q)))
        throw std::invalid_argument("DW potential requires finite v, lambda > 0 and finite q");
    const double shift = 1.5 * (q / lambda);
    return std::copysign(std::hypot(v, shift) + std::abs(shift), q);
}

inline double false_vacuum(double v, double phi_true)
{
    // Product of the two nonzero stationary points is -v^2. This form
    // avoids subtracting nearly equal numbers at large bias.
    return -v * (v / phi_true);
}

template <class T>
T value(const T &phi, double v, double lambda, double q, double phi_true)
{
    if (q == 0.0)
    {
        const T d = phi * phi - v * v;
        return 0.25 * lambda * d * d;
    }
    const double a = q / lambda;
    const T d = phi - phi_true;
    const T s = phi + phi_true - 2.0 * a;
    // V = U(phi)-U(phi_true), U=lambda/4*(phi^2-v^2)^2-q*phi^3.
    // Both terms in the bracket are nonnegative for the chosen global
    // minimum. No energy clipping or subtractive vacuum normalization.
    return 0.25 * lambda * d * d *
           (s * s + 2.0 * a * (phi_true - 2.0 * a));
}

template <class T>
T derivative(const T &phi, double v, double lambda, double q,
             double phi_true, double phi_false)
{
    if (q == 0.0)
        return lambda * phi * (phi * phi - v * v);
    return lambda * phi * (phi - phi_true) * (phi - phi_false);
}
} // namespace DomainWallPotential
#endif
