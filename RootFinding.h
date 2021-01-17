#pragma once
#include <tuple>
#include "wrappers.h"

namespace mathtool
{
    template <typename T, typename scalar = typename Type<T>::scalar, typename... Args>
    std::tuple<scalar, T> Newton(
        typename FuncType<T, Args...>::crT2T f,
        typename FuncType<T, Args...>::dcrT2T df,
        T x0,
        const scalar &tor,
        const Args... args)
    {
        scalar err;
        T p;
        do
        {
            p = linearSolve(df(x0, args...), f(x0, args...));
            err = norm(p);
            x0 -= p;
        } while (err > tor);
        err = norm(f(x0, args...));
        return std::make_tuple(err, x0);
    }
} // namespace mathtool