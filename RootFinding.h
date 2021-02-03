#pragma once
#include <tuple>
#include "wrappers.h"

namespace mathtool
{
    template <typename T, typename dT = typename Type<T>::derivative, typename scalar = typename Type<T>::scalar, typename... Args>
    inline scalar Newton_inplace_debug(
        typename FuncType<T, Args...>::rT2rT_inplace f,
        typename FuncType<T, Args...>::drT2rT_inplace df,
        T &x,
        T &y,
        dT &dy,
        const scalar &tor,
        Args&&... args)
    {
        scalar err;
        T p;
        do
        {
            f(x, y, std::forward<Args>(args)...);
            df(x, dy, std::forward<Args>(args)...);
            p = linearSolve(dy, y);
            x -= p;
            err = norm(y);
            std::cout << err << std::endl;
            if (err > 1e2)
            {
                std::cout << "Netwon Method fails." << std::endl; 
                break;
            }
        } while (err > tor);
        f(x, y, std::forward<Args>(args)...);
        err = norm(y);
        std::cout << "done!" << std::endl;
        return err;
    }

    template <typename T, typename dT = typename Type<T>::derivative, typename scalar = typename Type<T>::scalar, typename... Args>
    inline scalar Newton_inplace(
        typename FuncType<T, Args...>::rT2rT_inplace f,
        typename FuncType<T, Args...>::drT2rT_inplace df,
        T &x,
        T &y,
        dT &dy,
        const scalar &tor,
        Args&&... args)
    {
        scalar err;
        T p;
        do
        {
            f(x, y, std::forward<Args>(args)...);
            df(x, dy, std::forward<Args>(args)...);
            p = linearSolve(dy, y);
            x -= p;
            err = norm(y);
            if (err > 1e2)
            {
                std::cout << "Netwon Method fails." << std::endl; 
                break;
            }
        } while (err > tor);
        f(x, y, std::forward<Args>(args)...);
        err = norm(y);
        return err;
    }

    template <typename T, typename scalar = typename Type<T>::scalar, typename... Args>
    inline std::tuple<scalar, T> Newton(
        typename FuncType<T, Args...>::crT2T f,
        typename FuncType<T, Args...>::dcrT2T df,
        T x0,
        const scalar &tor,
        Args &&... args)
    {
        scalar err;
        T p;
        do
        {
            p = linearSolve(df(x0, std::forward<Args>(args)...), f(x0, std::forward<Args>(args)...));
            err = norm(p);
            x0 -= p;
            if (err > 1e2)
            {
                std::cout << "Netwon Method fails." << std::endl; 
                break;
            }
        } while (err > tor);
        err = norm(f(x0, std::forward<Args>(args)...));
        return std::make_tuple(err, x0);
    }
} // namespace mathtool