#pragma once
#include <tuple>
#include "wrappers.h"

namespace mathtool
{
    template <typename T, typename scalar = typename Type<T>::scalar, typename... Args>
    std::tuple<scalar, T> LMA(
        typename FuncType<T, Args...>::crT2T f,
        typename FuncType<T, Args...>::dcrT2T df,
        T x0,
        scalar a,
        const scalar &tor,
        const Args... args)
    {
        typedef typename std::result_of<typename FuncType<T, Args...>::dcrT2T(T, Args...)>::type dT;
        scalar err, res0, res1, mu;
        T p, x1, f0, f1, b;
        dT D = df(x0, args...);
        dT A = transpose(D) * D;
        dT I = identity<dT>(size(A) / 2, size(A) / 2);
        mu = 1.0;
        A += mu * I;
        f0 = f(x0, args...);
        b = D * f0;
        res0 = norm(f0);
        do
        {
            p = linearSolve(A, b);
            err = norm(p);
            x1 = x0 - p;
            f1 = f(x1, args...);
            res1 = norm(f1);
            if (res1 < res0)
            {
                mu *= a;
                x0 = x1;
                f0 = f1;
                res0 = res1;
                D = df(x0, args...);
                A = transpose(D) * D + mu * I;
                b = D * f0;
            }
            else
            {
                mu *= 2.0 - a;
                A += (1.0 - a) / (2.0 - a) * I;
            }
        } while (err > tor);
        return std::make_tuple(res0, x0);
    }
}; // namespace mathtool