#include "FArrayBox.H"
#include "MatterCCZ4RHS.hpp"
#include "RadiationField.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>

template <class T> using Vars = MatterCCZ4RHS<RadiationField>::Vars<T>;
template <class T> using Diff2 = MatterCCZ4RHS<RadiationField>::Diff2Vars<T>;
template <class V> V zero()
{
    V result;
    result.enum_mapping([](int, auto &value) { value = 0.; });
    return result;
}

int main()
{
    int checks = 0, failures = 0;
    auto check = [&](const char *label, double value, double expected, double tol = 3.e-13)
    {
        ++checks;
        if (!(std::isfinite(value) && std::abs(value - expected) <=
              tol * std::max(1., std::abs(expected))))
        {
            ++failures;
            std::cerr << label << ": got " << value << ", expected " << expected << std::endl;
        }
    };
    // Independent perfect-fluid projections, including nonzero momentum and
    // off-diagonal spatial stresses. E is not the rest-frame density.
    for (double speed : {0., 0.1, 0.5, 0.9})
        for (double Pi : {0.03, 1., 4.})
        {
            auto v = zero<Vars<double>>();
            auto d1 = zero<Vars<Tensor<1, double>>>();
            v.chi = 0.8;
            v.h[0][0] = 1.44; v.h[1][1] = 0.81; v.h[2][2] = 1.69;
            v.Pi_rad = Pi;
            const double kappa = 0.7;
            // Unit orthonormal direction (2/3,1/3,2/3), velocity=-Phi/Pi.
            const double direction[3] = {2./3., 1./3., 2./3.};
            double velocity[3];
            FOR(i)
            {
                velocity[i] = speed * direction[i] * std::sqrt(v.h[i][i] / v.chi);
                d1.psi_rad[i] = -Pi * velocity[i];
            }
            const auto inv = TensorAlgebra::compute_inverse_sym(v.h);
            const auto chris = TensorAlgebra::compute_christoffel(d1.h, inv);
            const auto em = RadiationField(kappa).compute_emtensor(v, d1, inv, chris.ULL);
            const double X = 0.5 * Pi * Pi * (1. - speed * speed);
            const double pressure = kappa * X * X;
            const double enthalpy_W2 = 4. * pressure / (1. - speed * speed);
            check("Eulerian energy", em.rho, enthalpy_W2 - pressure);
            check("radiation trace", em.S, em.rho);
            FOR(i)
            {
                check("fluid momentum", em.Si[i], enthalpy_W2 * velocity[i]);
                FOR(j)
                    check("fluid stress", em.Sij[i][j], enthalpy_W2 * velocity[i] * velocity[j] +
                          pressure * v.h[i][j] / v.chi);
            }
            for (double c : {0.3, 2., -1.5})
            {
                auto scaled_v = v;
                auto scaled_d1 = d1;
                scaled_v.Pi_rad *= c;
                FOR(i) { scaled_d1.psi_rad[i] *= c; }
                const auto scaled = RadiationField(kappa / std::pow(c, 4)).compute_emtensor(
                    scaled_v, scaled_d1, inv, chris.ULL);
                check("normalization energy", scaled.rho, em.rho);
                check("normalization trace", scaled.S, em.S);
                FOR(i)
                {
                    check("normalization momentum", scaled.Si[i], em.Si[i]);
                    FOR(j) { check("normalization stress", scaled.Sij[i][j], em.Sij[i][j]); }
                }
            }
        }

    auto v = zero<Vars<double>>();
    auto d1 = zero<Vars<Tensor<1, double>>>();
    auto d2 = zero<Diff2<Tensor<2, double>>>();
    const auto advec = zero<Vars<double>>();
    v.chi = v.lapse = v.Pi_rad = 1.;
    FOR(i) { v.h[i][i] = 1.; }
    d2.psi_rad[0][0] = -0.2;
    auto rhs = zero<Vars<double>>();
    RadiationField(1.).add_matter_rhs(rhs, v, d1, d2, advec);
    check("linear sound coefficient", rhs.Pi_rad, -0.2 / 3.);

    // Exact pointwise tensor transformation under X=F(x), with F'=f, F''=fp.
    // Nonzero lapse gradients and extrinsic curvature also transform here;
    // the PDE runs themselves use alpha=1 and Kij=0.
    v.Pi_rad = 1.2; v.lapse = 1.1; v.K = 0.12;
    v.A[0][0] = 0.02; v.A[1][1] = v.A[2][2] = -0.01;
    d1.psi_rad[0] = 0.23; d1.Pi_rad[0] = -0.17; d1.lapse[0] = 0.07;
    RadiationField(1.).add_matter_rhs(rhs, v, d1, d2, advec);
    double missing_chi_signal = 0.;
    for (double f : {0.7, 1.3})
    {
        const double fp = -0.4;
        const double chi = std::pow(f, -2./3.);
        auto mapped_v = v;
        auto mapped_d1 = d1;
        auto mapped_d2 = d2;
        auto mapped_rhs = zero<Vars<double>>();
        mapped_v.chi = chi;
        FOR(i)
        {
            const double factor = i == 0 ? f * f : 1.;
            mapped_v.h[i][i] = chi * factor;
            mapped_v.A[i][i] = chi * factor * v.A[i][i];
            mapped_d1.h[i][i][0] = mapped_v.h[i][i] * (i == 0 ? 4./3. : -2./3.) * fp/f;
        }
        mapped_d1.chi[0] = -2./3. * chi * fp/f;
        mapped_d1.psi_rad[0] *= f;
        mapped_d1.Pi_rad[0] *= f;
        mapped_d1.lapse[0] *= f;
        mapped_d2.psi_rad[0][0] = f*f*d2.psi_rad[0][0] + fp*d1.psi_rad[0];
        RadiationField(1.).add_matter_rhs(mapped_rhs, mapped_v, mapped_d1, mapped_d2, advec);
        check("coordinate invariant field RHS", mapped_rhs.psi_rad, rhs.psi_rad);
        check("physical connection and Kij RHS", mapped_rhs.Pi_rad, rhs.Pi_rad);
        // Sensitivity control: retaining only the conformal connection must
        // not accidentally pass the coordinate-invariance test.
        mapped_d1.chi[0] = 0.;
        RadiationField(1.).add_matter_rhs(mapped_rhs, mapped_v, mapped_d1, mapped_d2, advec);
        missing_chi_signal = std::max(missing_chi_signal, std::abs(mapped_rhs.Pi_rad - rhs.Pi_rad));
    }
    check("connection test detects missing chi gradient", missing_chi_signal > 1.e-3 ? 1. : 0., 1.);
    for (double c : {0.3, 2., -1.5})
    {
        auto sv = v;
        auto sd1 = d1;
        auto sd2 = d2;
        sv.Pi_rad *= c;
        FOR(i)
        {
            sd1.psi_rad[i] *= c;
            sd1.Pi_rad[i] *= c;
            FOR(j) { sd2.psi_rad[i][j] *= c; }
        }
        auto srhs = zero<Vars<double>>();
        RadiationField(1./std::pow(c,4)).add_matter_rhs(srhs, sv, sd1, sd2, advec);
        check("RHS normalization field", srhs.psi_rad, c * rhs.psi_rad);
        check("RHS normalization momentum", srhs.Pi_rad, c * rhs.Pi_rad);
    }
    std::cout << "RadiationMatterTest: " << checks-failures << "/" << checks
              << " passed; missing-chi sensitivity=" << missing_chi_signal << std::endl;
    return failures ? 1 : 0;
}
