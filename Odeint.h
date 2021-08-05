#pragma once
#include "types.h"
#include <tuple>

namespace mathtool
{
    typedef Eigen::Matrix<vec, -1, 1> vec_vec;

    class rk45 {
    public:
        static const vec A, C, CH, CT;
        static const mat B;
        static vec k0, k1, k2, k3, k4, k5;

        template <class ODEFun>
        static void one_step(double& dt, double& t, double& TE, vec& in, vec& out, ODEFun& f)
        {
            k0 = dt * f(t + A(0) * dt, in);
            k1 = dt * f(t + A(1) * dt, in + B(1, 0) * k0);
            k2 = dt * f(t + A(2) * dt, in + B(2, 0) * k0 + B(2, 1) * k1);
            k3 = dt * f(t + A(3) * dt, in + B(3, 0) * k0 + B(3, 1) * k1 + B(3, 2) * k2);
            k4 = dt * f(t + A(4) * dt, in + B(4, 0) * k0 + B(4, 1) * k1 + B(4, 2) * k2 + B(4, 3) * k3);
            k5 = dt * f(t + A(5) * dt, in + B(5, 0) * k0 + B(5, 1) * k1 + B(5, 2) * k2 + B(5, 3) * k3 + B(5, 4) * k4);
            out = in + CH(0) * k0 + CH(1) * k1 + CH(2) * k2 + CH(3) * k3 + CH(4) * k4 + CH(5) * k5;
            TE = (CT(0) * k0 + CT(1) * k1 + CT(2) * k2 + CT(3) * k3 + CT(4) * k4 + CT(5) * k5).norm();
        }

        template <class ODEFun>
        static auto solve(double dt, double t0, double t_end, vec& y0, ODEFun& f, double tor = 1e-8)
        {
            vec_vec Y(10);
            vec t(10);
            int size = 10;
            int Nt = 1;
            Y(0) = y0;
            t(0) = t0;
            double TE = 0;
            vec yt(y0.size());
            while (t(Nt - 1) < t_end)
            {
                one_step(dt, t(Nt - 1), TE, Y(Nt - 1), yt, f);
                if (TE > tor)
                {
                    dt = 0.9 * dt * pow(tor / TE, 0.2);
                }
                else
                {
                    ++Nt;
                    if (size < Nt)
                    {
                        size *= 2;
                        t.conservativeResize(size);
                        Y.conservativeResize(size);
                    }
                    t(Nt - 1) = t(Nt - 2) + dt;
                    Y(Nt - 1) = yt;
                    dt = 0.9 * dt * pow(tor / TE, 0.2);
                }
            }
            t.conservativeResize(Nt);
            Y.conservativeResize(Nt);
            auto sol = std::make_tuple(t, Y);
            return sol;
        }
    };
    const vec rk45::A = (vec(6) << 0, 2.0 / 9.0, 1.0 / 3.0, 3.0 / 4.0, 1.0, 5.0 / 6.0).finished();
    const mat rk45::B = (mat(6, 5)
        << 0, 0, 0, 0, 0,
        2.0 / 9.0, 0, 0, 0, 0,
        1.0 / 12.0, 1.0 / 4.0, 0, 0, 0,
        69.0 / 128, -243.0 / 128.0, 135.0 / 64, 0, 0,
        -17.0 / 12.0, 27.0 / 4.0, -27.0 / 5.0, 16.0 / 15.0, 0,
        65.0 / 432.0, -5.0 / 16.0, 13.0 / 16.0, 4.0 / 27.0, 5.0 / 144.0).finished();
    const vec rk45::C = (vec(5) << 1.0 / 9.0, 0, 9.0 / 20.0, 16.0 / 45.0, 1.0 / 12.0).finished();
    const vec rk45::CH = (vec(6) << 47.0 / 450.0, 0, 12.0 / 25.0, 32.0 / 225.0, 1.0 / 30.0, 6.0 / 25.0).finished();
    const vec rk45::CT = (vec(6) << -1.0 / 150.0, 0, 3.0 / 100, -16.0 / 75.0, -1.0 / 20.0, 6.0 / 25.0).finished();
    vec rk45::k0, rk45::k1, rk45::k2, rk45::k3, rk45::k4, rk45::k5;
};