#ifndef DIAGNOSTICVARIABLES_HPP
#define DIAGNOSTICVARIABLES_HPP
#include <array>
#include <string>
enum
{
    c_Ham, c_Mom, c_rho_rad, c_rho_wall, c_S_minus_E_rad,
    c_sqrt_gamma, c_chi_proper, c_current_proper, c_timelike_margin,
    c_Ham_abs, c_Mom1, c_Mom2, c_Mom3, c_Mom_abs, c_Ham_wall,
    c_rho_total, c_hybrid_defect,
    NUM_DIAGNOSTIC_VARS
};
namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {
    "Ham", "Mom", "rho_rad", "rho_wall", "S_minus_E_rad", "sqrt_gamma",
    "chi_proper", "current_proper", "timelike_margin",
    "Ham_abs", "Mom1", "Mom2", "Mom3", "Mom_abs", "Ham_wall", "rho_total", "hybrid_defect"};
}
#endif
