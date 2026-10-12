#ifndef USERVARIABLES_HPP
#define USERVARIABLES_HPP
#include "ArrayTools.hpp"
#include "CCZ4UserVariables.hpp"
#include "DiagnosticVariables.hpp"
enum
{
    c_phi = NUM_CCZ4_VARS,
    c_Pi,
    c_psi_rad,
    c_Pi_rad,
    NUM_VARS
};
static_assert(NUM_VARS == 29, "RadiationFLRW checkpoint layout changed");
namespace UserVariables
{
static const std::array<std::string, NUM_VARS - NUM_CCZ4_VARS>
    user_variable_names = {"phi", "Pi", "psi_rad", "Pi_rad"};
static const std::array<std::string, NUM_VARS> variable_names =
    ArrayTools::concatenate(ccz4_variable_names, user_variable_names);
}
#include "UserVariables.inc.hpp"
#endif
