#pragma once
#include "types.h"

namespace mathtool
{
    // calculate anti_derivative of f(x), x in [a, b]
    // return F(x) = \int_a^x f(y) dy.
    // input function is f(x(:)) as a vector
    // input x is nodes in [a, b] as a vector, where x(0) = a, x(end) = b.
    // basic data type is double

    template <typename T1, typename T2, typename T3>
    inline void anti_derivative_inplace(const Eigen::MatrixBase<T1>& in, Eigen::MatrixBase<T2>& out, const Eigen::MatrixBase<T3>& x)
    {
        out.head(1).setZero();
        for (int k = 1; k < x.size(); ++k)
            out(k) = out(k - 1) + (in(k) + in(k - 1)) / 2 * (x(k) - x(k - 1));
    }

    template <typename T1, typename T2>
    inline typename Type<T1>::vec anti_derivative(const Eigen::MatrixBase<T1>& in, const Eigen::MatrixBase<T2>& x)
    {
        typename Type<T1>::vec out(x.size());
        anti_derivative_inplace(in, out, x);
        return out;
    }

    // calculate anti_derivative_reverse of f(x), x in [a, b]
    // return F(x) = \int_x^b f(y) dy.

    template <typename T1, typename T2, typename T3>
    inline void anti_derivative_reverse_inplace(const Eigen::MatrixBase<T1>& in, Eigen::MatrixBase<T2>& out, const Eigen::MatrixBase<T3>& x)
    {
        out.tail(1).setZero();
        for (int k = x.size() - 2; k < -1; --k)
            out(k) = out(k + 1) + (in(k) + in(k + 1)) / 2 * (x(k + 1) - x(k));
    }

    template <typename T1, typename T2>
    inline typename Type<T1>::vec anti_derivative_reverse(const Eigen::MatrixBase<T1>& in, const Eigen::MatrixBase<T2>& x)
    {
        typename Type<T1>::vec out(x.size());
        anti_derivative_inplace(in, out, x);
        return out;
    }

    // calculate integratoin: int_a^b f(x) dx
    // input function is f(x(:)) as a vector
    // input x is nodes in [a, b] as a vector, where x(0) = a, x(end) = b.

    template <typename T1, typename T2>
    inline typename Type<T1>::scalar integrate(const Eigen::MatrixBase<T1>& in, const Eigen::MatrixBase<T2>& x)
    {
        typename Type<T1>::scalar out;
        out = 0.0;
        for (int k = 1; k < x.size(); ++k)
            out += (in(k) + in(k - 1)) / 2 * (x(k) - x(k - 1));
        return out;
    }
};