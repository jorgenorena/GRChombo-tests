#ifndef DIAGNOSTICVARIABLES_HPP
#define DIAGNOSTICVARIABLES_HPP
#include <array>
#include <string>
enum
{
    c_Ham, c_Mom, c_rho_rad, c_rho_wall, c_S_minus_E_rad,
    c_sqrt_gamma, c_chi_proper, c_current_proper, c_timelike_margin,
    NUM_DIAGNOSTIC_VARS
};
namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {
    "Ham", "Mom", "rho_rad", "rho_wall", "S_minus_E_rad", "sqrt_gamma",
    "chi_proper", "current_proper", "timelike_margin"};
}
#endif
