// function wrappers.
// default to use eigen for vectors/matrix calculation. If other matrix class is used, for example Armadillo, please implement all wrapper functions.
#pragma once
#include <cmath>
#include "type.h"

namespace mathtool
{
    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value, typename Type<T>::scalar>::type
    norm(const T &x) { return std::fabs(x); }

    template <typename T>
    typename std::enable_if<std::is_integral<T>::value || is_complex<T>::value, typename Type<T>::scalar>::type
    norm(const T &x) { return std::abs(x); }

    template <typename TA, typename Tb>
    typename std::enable_if<is_numeric<TA>::value, Tb>::type
    linearSolve(const TA &A, const Tb &b) { return b / A; }

    template <typename T>
    typename std::enable_if<is_numeric<T>::value, int>::type
    size(const T &x) { return 1; }

    template <typename T>
    typename std::enable_if<is_numeric<T>::value, T>::type
    identity(const int &r, const int &i) { return T(1); }
}; // namespace mathtool

#ifdef USE_OTHER_MATRIX_LIB
#else
namespace mathtool
{
    template <typename T>
    typename std::enable_if<is_matrix_related<T>::value, typename Type<T>::scalar>::type
    norm(const T &x) { return x.norm(); }

    template <typename TA, typename Tb>
    Tb
    linearSolve(const Eigen::MatrixBase<TA> &A, const Eigen::MatrixBase<Tb> &b) { return A.colPivHouseholderQr().solve(b); }

    template <typename T>
    auto
    transpose(const Eigen::DenseBase<T> &x) { return x.transpose(); }

    template <typename T>
    auto
    size(const Eigen::EigenBase<T> &x) { return x.size(); }

    template <typename T>
    typename std::enable_if<is_matrix_related<T>::value, T>::type
    identity(const int &r, const int &c) { return T::Identity(r, c); }

    template <typename T>
    typename Type<T>::scalar
    det(const Eigen::MatrixBase<T> &x) { return x.determinant(); }
}; // namespace mathtool
#endif