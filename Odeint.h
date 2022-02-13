#pragma once
#include "types.h"
#include <tuple>

namespace mathtool
{
    typedef Eigen::Matrix<vec, -1, 1> vec_vec;

    template <class T>
    class rk45 
    {
    public:
        static const vec A, C, CH, CT;
        static const mat B;
        static T k0, k1, k2, k3, k4, k5;

        template <class ODEFun>
        static void one_step(double& dt, double& t, double& TE, T& in, T& out, ODEFun& f)
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
        static auto solve(double dt, double t0, double t_end, T& y0, ODEFun& f, double tor = 1e-8)
        {
            Eigen::Matrix<T, -1, 1> Y(10);
            vec t(10);
            int size = 10;
            int Nt = 1;
            Y(0) = y0;
            t(0) = t0;
            double TE = 0;
            T yt = y0;
            while (t(Nt - 1) < t_end)
            {
                if (dt < 1e-12)
                {
                    std::cout << "min time step 1e-12 is reached, system may blow up" << std::endl;
                    break;
                }
                one_step(dt, t(Nt - 1), TE, Y(Nt - 1), yt, f);
                if (TE > tor)
                {
                    dt = 0.9 * dt * pow(tor / TE, 0.2);
                    std::cout << "t = " << t(Nt - 1) << ". current try fails, new dt = " << dt << std::endl;
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
                    std::cout << "t = " << t(Nt - 1) << ". current try accepts, next dt = " << dt << std::endl;
                }
            }
            t.conservativeResize(Nt);
            Y.conservativeResize(Nt);
            auto sol = std::make_tuple(t, Y);
            return sol;
        }

        template <class ODEFun>
        static auto solve(double dt, double t0, double t_end, vec& t, T& y0, ODEFun& f, double tor = 1e-8)
        {
            int Nt = t.size();
            Eigen::Matrix<T, -1, 1> Y(Nt);
            int size = 10;
            Y(0) = y0;
            t(0) = t0;
            double TE = 0;
            double time = t0;
            T yc = y0;
            T yt = y0;
            int ind = 1;
            while (time < t_end)
            {
                if (dt < 1e-12)
                {
                    std::cout << "min time step 1e-12 is reached, system may blow up" << std::endl;
                    break;
                }
                one_step(dt, time, TE, yc, yt, f);
                if (TE > tor)
                {
                    dt = 0.9 * dt * pow(tor / TE, 0.2);
                    std::cout << "t = " << time << ". current try fails, new dt = " << dt << std::endl;
                }
                else
                {
                    std::cout << "t = " << time;
                    time += dt;
                    dt = 0.9 * dt * pow(tor / TE, 0.2);
                    if (t(ind) < time)
                    {
                        time -= dt;
                        dt = t(ind) - time;
                        one_step(dt, time, TE, yc, yt, f);
                        time = t(ind);
                        Y(ind) = yt;
                        //t(ind) = time;
                        //Y(ind) = yc;
                        ++ind;
                    }
                    yc = yt;
                    std::cout << ". current try accepts, next dt = " << dt << std::endl;
                }
            }
            auto sol = std::make_tuple(t, Y);
            return sol;
        }
    };
    template <class T>
    const vec rk45<T>::A = (vec(6) << 0, 2.0 / 9.0, 1.0 / 3.0, 3.0 / 4.0, 1.0, 5.0 / 6.0).finished();
    template <class T>
    const mat rk45<T>::B = (mat(6, 5)
        << 0, 0, 0, 0, 0,
        2.0 / 9.0, 0, 0, 0, 0,
        1.0 / 12.0, 1.0 / 4.0, 0, 0, 0,
        69.0 / 128, -243.0 / 128.0, 135.0 / 64, 0, 0,
        -17.0 / 12.0, 27.0 / 4.0, -27.0 / 5.0, 16.0 / 15.0, 0,
        65.0 / 432.0, -5.0 / 16.0, 13.0 / 16.0, 4.0 / 27.0, 5.0 / 144.0).finished();
    template <class T>
    const vec rk45<T>::C = (vec(5) << 1.0 / 9.0, 0, 9.0 / 20.0, 16.0 / 45.0, 1.0 / 12.0).finished();
    template <class T>
    const vec rk45<T>::CH = (vec(6) << 47.0 / 450.0, 0, 12.0 / 25.0, 32.0 / 225.0, 1.0 / 30.0, 6.0 / 25.0).finished();
    template <class T>
    const vec rk45<T>::CT = (vec(6) << -1.0 / 150.0, 0, 3.0 / 100, -16.0 / 75.0, -1.0 / 20.0, 6.0 / 25.0).finished();
    template <class T> T rk45<T>::k0;
    template <class T> T rk45<T>::k1;
    template <class T> T rk45<T>::k2;
    template <class T> T rk45<T>::k3;
    template <class T> T rk45<T>::k4;
    template <class T> T rk45<T>::k5;
};