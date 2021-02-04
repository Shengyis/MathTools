#pragma once
#include "wrappers.h"

namespace mathtool
{
    // calculate anti_derivative of f(x), x in [a, b]
    // return F(x) = \int_a^x f(y) dy.
    // input function is f(x(:)) as a vector
    // input x is nodes in [a, b] as a vector, where x(0) = a, x(end) = b.
    // basic data type is double
    inline void anti_derivative_inplace(const vec_d &in, vec_d &out, const vec_d &x)
    {
        out(0) = 0.0;
        for (int k = 1; k < x.size(); ++k)
            out(k) = out(k-1) + (in(k) + in(k-1)) / 2.0 * (x(k) - x(k-1)); 
    }

    inline vec_d anti_derivative(const vec_d &in, const vec_d &x)
    {
        vec_d out = vec_d::Zero(x.size());
        anti_derivative_inplace(in, out, x);
        return out;
    }
};