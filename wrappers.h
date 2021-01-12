// function wrappers.
// default to use eigen for vectors/matrix calculation. If other matrix class is used, for example Armadillo, please implement all wrapper functions.
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
    T linearSolve(const T &A, const T &b)
    {
        return b / A;
    }

    template <typename T>
    auto norm(const Eigen::MatrixBase<T> &x)
    {
        return x.norm();
    }

    template <typename T>
    auto transpose(const Eigen::MatrixBase<T> &x)
    {
        return x.transpose();
    }

    template <typename T>
    auto operator%(const Eigen::MatrixBase<T> &x, const Eigen::MatrixBase<T> &y)
    {
        return x.cwiseProduct(y);
    }

    template <typename T>
    auto operator/(const Eigen::MatrixBase<T> &x, const Eigen::MatrixBase<T> &y)
    {
        return x.cwiseQuotient(y);
    }

    template <typename T>
    Eigen::Vector<T, -1> linearSolve(const Eigen::Matrix<T, -1, -1> &A, const Eigen::Vector<T, -1> &b)
    {
        return A.colPivHouseholderQr().solve(b);
    }
}; // namespace mathtool
