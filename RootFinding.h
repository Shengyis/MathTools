#pragma once
#include <tuple>
#include "wrappers.h"

namespace mathtool
{
    template <typename F, typename DF, typename Tx, typename Te, typename... Args>
    std::tuple<Te, Tx> Newton(F f, DF df, Tx x0, const Te &tor, const Args... args)
    {
        Te err;
        Tx p;
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