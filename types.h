#pragma once
#ifndef EIGEN_MATRIXBASE_PLUGIN
#define EIGEN_MATRIXBASE_PLUGIN <MathTools/MatrixBaseAddons.h>
#endif

#ifndef PREC
#define PREC 128
#endif

#include <type_traits>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/numeric/interval.hpp>
#include <eigen3/Eigen/Eigen>

using namespace std::complex_literals;

namespace mathtool
{
    typedef typename boost::multiprecision::number<boost::multiprecision::cpp_dec_float<PREC>> mp_f;
    static mp_f mp_pi = boost::math::constants::pi<mp_f>();
    static double pi = boost::math::constants::pi<double>();

    template <typename T>
    struct is_matrix
        : std::integral_constant<bool, std::is_base_of<Eigen::EigenBase<T>, T>::value>
    {
    };

    template <typename T>
    struct is_mp_f
        : std::integral_constant<bool, std::is_same<mp_f, T>::value>
    {
    };

    template <typename T>
    struct is_mp_i
        : std::false_type
    {
    };

    template <typename T>
    struct is_mp_ra
        : std::false_type
    {
    };

    template <typename T>
    struct is_mp
        : std::integral_constant<bool, is_mp_f<T>::value || is_mp_i<T>::value || is_mp_ra<T>::value>
    {
    };

    template <typename T>
    struct is_real
        : std::integral_constant<bool, std::is_arithmetic<T>::value || is_mp<T>::value>
    {
    };

    template <typename T>
    struct is_complex : std::false_type
    {
    };

    template <typename T>
    struct is_complex<std::complex<T>>
        : std::integral_constant<bool, is_real<T>::value>
    {
    };

    template <typename T>
    struct is_scalar
        : std::integral_constant<bool, is_real<T>::value || is_complex<T>::value>
    {
    };

    template <typename T, bool = is_scalar<T>::value>
    struct Scalar
    {
        typedef T type;
    };

    template <typename T>
    struct Scalar<T, 0>
    {
        static_assert(is_matrix<T>::value, "T is neither a scalar type nor a matrix type.");
        typedef typename T::Scalar type;
    };

    template <typename T, bool = is_scalar<T>::value>
    struct RealScalar
    {
        typedef T type;
    };

    template <typename T>
    struct RealScalar<T, 0>
    {
        static_assert(is_matrix<T>::value, "T is neither a scalar type nor a matrix type.");
        typedef typename T::RealScalar type;
    };

    template <typename T>
    struct RealScalar<std::complex<T>, 1>
    {
        typedef T type;
    };

    template <typename T>
    struct Type;
    template <typename T>
    struct Mat;
    template <typename T>
    struct Vec;
    template <typename T>
    struct Mat_c;
    template <typename T>
    struct Vec_c;

    template <typename T>
    struct Type
    {
        static_assert(is_scalar<T>::value || is_matrix<T>::value, "T is neither a scalar type nor a matrix type.");
        typedef typename Scalar<T>::type scalar;
        typedef typename RealScalar<T>::type real_scalar;
        typedef typename std::complex<typename RealScalar<T>::type> complex_scalar;
        typedef typename Mat<T>::type mat;
        typedef typename Vec<T>::type vec;
        typedef typename Mat_c<T>::type mat_c;
        typedef typename Vec_c<T>::type vec_c;
    };

    template <typename T>
    struct Mat
    {
        typedef typename Eigen::Matrix<typename Type<T>::scalar, -1, -1> type;
    };

    template <typename T>
    struct Vec
    {
        typedef typename Eigen::Matrix<typename Type<T>::scalar, -1, 1> type;
    };

    template <typename T>
    struct Mat_c
    {
        typedef typename Eigen::Matrix<typename Type<T>::complex_scalar, -1, -1> type;
    };

    template <typename T>
    struct Vec_c
    {
        typedef typename Eigen::Matrix<typename Type<T>::complex_scalar, -1, 1> type;
    };

    typedef std::complex<double> cd;
    typedef Eigen::MatrixXd mat;
    typedef Eigen::MatrixXi mat_i;
    typedef Eigen::MatrixXcd mat_cd;
    typedef Eigen::VectorXd vec;
    typedef Eigen::VectorXi vec_i;
    typedef Eigen::VectorXcd vec_cd;

    typedef std::complex<mp_f> mp_c;
    typedef Eigen::Matrix<mp_f, -1, -1> mp_mat;
    typedef Eigen::Matrix<mp_f, -1, 1> mp_vec;
    typedef Eigen::Matrix<mp_c, -1, -1> mp_mat_c;
    typedef Eigen::Matrix<mp_c, -1, 1> mp_vec_c;
}; // namespace mathtool

namespace boost
{
    namespace numeric
    {
        namespace interval_lib
        {
            namespace constants
            {
                template <>
                inline mathtool::mp_f pi_lower<mathtool::mp_f>() { return mathtool::mp_pi; }
                template <>
                inline mathtool::mp_f pi_upper<mathtool::mp_f>() { return mathtool::mp_pi; }
                template <>
                inline mathtool::mp_f pi_twice_lower<mathtool::mp_f>() { return mathtool::mp_pi * 2; }
                template <>
                inline mathtool::mp_f pi_twice_upper<mathtool::mp_f>() { return mathtool::mp_pi * 2; }
                template <>
                inline mathtool::mp_f pi_half_lower<mathtool::mp_f>() { return mathtool::mp_pi / 2; }
                template <>
                inline mathtool::mp_f pi_half_upper<mathtool::mp_f>() { return mathtool::mp_pi / 2; }
            };
        };
    };
};
