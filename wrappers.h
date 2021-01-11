#pragma once
#ifdef USE_EIGEN
#include <eigen3/Eigen/Dense>
#endif

namespace mathtool
{
    template <typename T>
    auto norm(const T &x)
    {
        return std::fabs(x);
    }

    template <typename T>
    auto norm(const Eigen::Vector<T, -1> &x)
    {
        return x.norm();
    }

    template <typename T>
    T linearSolve(const T &A, const T &b)
    {
        return b / A;
    }

    template <typename T>
    Eigen::Vector<T, -1> linearSolve(const Eigen::Matrix<T, -1, -1> &A, const Eigen::Vector<T, -1> &b)
    {
        return A.colPivHouseholderQr().solve(b);
    }
}; // namespace mathtool