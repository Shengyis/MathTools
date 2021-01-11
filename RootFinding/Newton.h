#pragma once
#include <tuple>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/QR>

namespace mathtool
{
    template <typename Tx, typename Ty, typename... Args>
    struct Type
    {
        typedef Ty (*pF)(Tx, Args...);
    };

    template <typename Te, typename Tx, typename... Args>
    std::tuple<Te, Eigen::Vector<Tx, -1>> Newton(
        typename Type<Eigen::Vector<Tx, -1>, Eigen::Vector<Tx, -1>, Args...>::pF f,
        typename Type<Eigen::Vector<Tx, -1>, Eigen::Matrix<Tx, -1, -1>, Args...>::pF df,
        Eigen::Vector<Tx, -1> x0,
        const Te &tor,
        const Args... args)
    {
        Te err;
        Eigen::Vector<Tx, -1> p;
        do
        {
            p = df(x0, args...).colPivHouseholderQr().solve(f(x0, args...));
            err = p.norm();
            x0 -= p;
        } while (err > tor);
        err = f(x0, args...).norm();
        return std::make_tuple(err, x0);
    }

    template <typename Te, typename Tx, typename... Args>
    std::tuple<Te, Tx> Newton(
        typename Type<Tx, Tx, Args...>::pF f,
        typename Type<Tx, Tx, Args...>::pF df,
        Tx x0,
        const Te &tor,
        const Args... args)
    {
        Te err;
        Tx p;
        do
        {
            p = f(x0, args...) / df(x0, args...);
            err = std::fabs(p);
            x0 -= p;
        } while (err > tor);
        err = std::fabs(f(x0, args...));
        return std::make_tuple(err, x0);
    }
} // namespace mathtool