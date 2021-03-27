#pragma once
#include "types.h"
#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
#include <boost/math/interpolators/barycentric_rational.hpp>

namespace mathtool
{
    enum interp_method
    {
        cubic,
        rational
    };

    template <interp_method method, typename real>
    struct InterpMethod
    {
    };

    template <typename real>
    struct InterpMethod<cubic, real>
    {
        typedef boost::math::interpolators::cardinal_cubic_b_spline<real> plan_type;
        static auto get_plan(const typename Type<real>::vec &x, const typename Type<real>::vec &y,
                             real left_endpoint_derivative = std::numeric_limits<real>::quiet_NaN(),
                             real right_endpoint_derivative = std::numeric_limits<real>::quiet_NaN())
        {
            return boost::math::interpolators::cardinal_cubic_b_spline<real>(y.data(), y.size(), x(0), x(1) - x(0),
                                                                             left_endpoint_derivative, right_endpoint_derivative);
        }
    };

    template <typename real>
    struct InterpMethod<rational, real>
    {
        typedef boost::math::barycentric_rational<real> plan_type;
        static auto get_plan(const typename Type<real>::vec &x, const typename Type<real>::vec &y,
                             real left_endpoint_derivative = std::numeric_limits<real>::quiet_NaN(),
                             real right_endpoint_derivative = std::numeric_limits<real>::quiet_NaN())
        {
            return boost::math::barycentric_rational<real>(x.data(), y.data(), x.size());
        }
    };

    template <interp_method method, typename T>
    inline typename InterpMethod<method, typename Type<T>::scalar>::plan_type interp1d_plan(
        const T &x, const T &y,
        typename Type<T>::scalar left_endpoint_derivative = std::numeric_limits<typename Type<T>::scalar>::quiet_NaN(),
        typename Type<T>::scalar right_endpoint_derivative = std::numeric_limits<typename Type<T>::scalar>::quiet_NaN())
    {
        return InterpMethod<method, typename Type<T>::scalar>::get_plan(x, y, left_endpoint_derivative, right_endpoint_derivative);
    }

    template <typename Tx, typename Tout, typename PlanType>
    inline void interp1d_inplace(const Eigen::MatrixBase<Tx> &x, Eigen::MatrixBase<Tout> &out, PlanType &plan)
    {
        #pragma omp parallel for
        for (int k = 0; k < out.size(); ++k)
            out(k) = plan(x(k));
    }

    template <typename T, typename PlanType>
    inline typename std::enable_if<is_real<T>::value, void>::type interp1d_inplace(const T &x, T &out, PlanType &plan)
    {
        out = plan(x);
    }

    template <typename T, typename PlanType>
    inline typename Type<T>::vec interp1d(const Eigen::MatrixBase<T> &x, PlanType &plan)
    {
        typename Type<T>::vec out(x.size());
        interp1d_inplace(x, out, plan);
        return out;
    }

    template <typename T, typename PlanType>
    inline typename std::enable_if<is_real<T>::value, T>::type interp1d(const T &x, PlanType &plan)
    {
        T out;
        interp1d_inplace(x, out, plan);
        return out;
    }
};