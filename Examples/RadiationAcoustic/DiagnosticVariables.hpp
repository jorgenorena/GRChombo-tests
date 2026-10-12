#ifndef DIAGNOSTICVARIABLES_HPP
#define DIAGNOSTICVARIABLES_HPP
#include <array>
#include <string>
// The first nine entries also support the shared Pass A diagnostic/check class.
enum
{
    c_Ham, c_Mom, c_rho_rad, c_rho_wall, c_S_minus_E_rad,
    c_sqrt_gamma, c_chi_proper, c_current_proper, c_timelike_margin,
    c_energy_proper, c_rho_rest, c_Srad1, c_Srad11, c_Srad22, c_Srad33,
    c_Pi_cos, c_Pi_sin, c_psi_cos, c_psi_sin, c_grad_Pi, c_grad_velocity,
    NUM_DIAGNOSTIC_VARS
};
namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {
    "Ham", "Mom", "rho_rad", "rho_wall", "S_minus_E_rad", "sqrt_gamma",
    "chi_proper", "current_proper", "timelike_margin", "energy_proper",
    "rho_rest", "Srad1", "Srad11", "Srad22", "Srad33", "Pi_cos", "Pi_sin",
    "psi_cos", "psi_sin", "grad_Pi", "grad_velocity"};
}
#endif
