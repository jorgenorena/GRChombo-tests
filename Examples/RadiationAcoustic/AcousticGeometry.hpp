#ifndef ACOUSTICGEOMETRY_HPP_
#define ACOUSTICGEOMETRY_HPP_
#include "simd.hpp"
#include <cmath>

// A periodic relabelling of flat space: X=x+b*sin(k*x)/k, Y=y, Z=z.
// f=dX/dx is positive. gamma=diag(f^2,1,1), chi=f^(-2/3).
struct AcousticGeometry
{
    double k, amplitude;
    template <class T> T X(const T &x) const
    {
        return x + amplitude * sin(k * x) / k;
    }
    template <class T> T jacobian(const T &x) const
    {
        return 1. + amplitude * cos(k * x);
    }
    template <class T> T jacobian_derivative(const T &x) const
    {
        return -amplitude * k * sin(k * x);
    }
};
#endif
