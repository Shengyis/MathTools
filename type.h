#pragma once
#include <type_traits>
#include <complex>

using namespace std::complex_literals;
typedef std::complex<double> cd;

// declaration part
namespace mathtool
{
    template <typename T>
    struct is_complex;

    template <typename T>
    struct is_matrix_related;

    template <typename T>
    struct is_numeric;

    template <typename T>
    struct is_calculable;

    template <typename T, bool>
    struct Scalar;

    template <typename T, bool>
    struct RealScalar;

    template <typename T, bool>
    struct Derivative;

    template <typename T>
    struct ToMat;

    template <typename T>
    struct ToVec;

    template <typename T>
    struct Type;

    template <typename T, typename... Args>
    struct FuncType;
}; // namespace mathtool

// implement
namespace mathtool
{
    template <typename T>
    struct is_complex : std::false_type
    {
    };

    template <typename T>
    struct is_complex<std::complex<T>> : std::true_type
    {
    };

    template <typename T>
    struct is_numeric : std::integral_constant<bool, is_complex<T>::value || std::is_arithmetic<T>::value>
    {
    };

    template <typename T>
    struct is_calculable : std::integral_constant<bool, is_numeric<T>::value || is_matrix_related<T>::value>
    {
    };

    template <typename T, bool = is_numeric<T>::value>
    struct Scalar
    {
        typedef T type;
    };

    template <typename T, bool = is_numeric<T>::value>
    struct RealScalar
    {
        typedef T type;
    };

    template <typename T>
    struct RealScalar<std::complex<T>, 1>
    {
        typedef T type;
    };

    template <typename T, bool = is_numeric<T>::value>
    struct Derivative
    {
        typedef T type;
    };

    template <typename T>
    struct Type
    {
        static_assert(is_calculable<T>::value, "T is not a calculable type.");
        typedef typename Scalar<T>::type scalar;
        typedef typename RealScalar<T>::type real_scalar;
        typedef typename Derivative<T>::type derivative;
    };

}; // namespace mathtool

#ifdef USE_OTHER_MATRIX_LIB
#else
#define EIGEN_MATRIXBASE_PLUGIN <MathTools/MatrixBaseAddons.h>
#include <eigen3/Eigen/Dense>
namespace mathtool
{
    template <typename T>
    struct is_matrix_related : std::integral_constant<bool, std::is_base_of<Eigen::EigenBase<T>, T>::value>
    {
    };

    template <typename T>
    struct Scalar<T, 0>
    {
        typedef typename T::value_type type;
    };

    template <typename T>
    struct RealScalar<T, 0>
    {
        typedef typename T::RealScalar type;
    };

    template <typename T>
    struct Derivative<T, 0>
    {
        typedef typename Eigen::Matrix<typename T::value_type, -1, -1> type;
    };

    template <typename T>
    struct ToMat
    {
        typedef typename Eigen::Matrix<typename Type<T>::scalar, -1, -1> type;
    };

    template <typename T>
    struct ToVec
    {
        typedef typename Eigen::Vector<typename Type<T>::scalar, -1> type;
    };
}; // namespace mathtool
#endif

namespace mathtool
{
    template <typename T, typename... Args>
    struct FuncType
    {
        typedef void (*rT2rT_inplace)(T &, T &, Args&&...);
        typedef T (*crT2T)(const T &, Args&&...);
        typedef void (*drT2rT_inplace)(T &, typename Type<T>::derivative &, Args&&...);
        typedef typename Type<T>::derivative (*dcrT2T)(const T &, Args&&...);
    };
}; // namespace mathtool