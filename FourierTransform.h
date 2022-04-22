#pragma once
#include "types.h"
#include "fft_helper.h"

namespace mathtool
{
    template <typename T>
    inline auto fft(const Eigen::MatrixBase<T>& in)
    {
        typedef FFT_PLAN<typename Type<T>::real_scalar> FFTP;
        return FFTP::fft.fwd(in);
    }

    template <typename T>
    inline auto ifft(const Eigen::MatrixBase<T>& in)
    {
        typedef FFT_PLAN<typename Type<T>::real_scalar> FFTP;
        return FFTP::fft.inv(in);
    }

    template <typename T>
    inline auto fft2(const Eigen::MatrixBase<T>& in)
    {
        int N1 = in.rows();
        int N2 = in.cols();
        mat_cd tmp(N1, N2), out(N2, N1);

        int k = 0;
        tmp.col(k) = fft(in.col(k).eval());
#pragma omp parallel for
        for (k = 1; k < N2; ++k)
        {
            tmp.col(k) = fft(in.col(k).eval());
        }
        tmp.transposeInPlace();
        k = 0;
        out.col(k) = fft(tmp.col(k).eval());
#pragma omp parallel for
        for (k = 1; k < N1; ++k)
        {
            out.col(k) = fft(tmp.col(k).eval());
        }
        out.transposeInPlace();
        return out;
    }

    template <typename T>
    inline auto ifft2(const Eigen::MatrixBase<T>& in)
    {
        int N1 = in.rows();
        int N2 = in.cols();
        mat_cd tmp(N1, N2);
        mat out(N2, N1);

        int k = 0;
        tmp.col(k) = ifft(in.col(k).eval());
#pragma omp parallel for
        for (k = 1; k < N2; ++k)
        {
            tmp.col(k) = ifft(in.col(k).eval());
        }
        tmp.transposeInPlace();
        k = 0;
        out.col(k) = ifft(tmp.col(k).eval());
#pragma omp parallel for
        for (k = 1; k < N1; ++k)
        {
            out.col(k) = ifft(tmp.col(k).eval());
        }
        out.transposeInPlace();
        return out;
    }

    /* template <typename T1, typename T2>
    inline void fft2_inplace(const Eigen::MatrixBase<T1>& in, Eigen::MatrixBase<T2>& out)
    {
        static_assert(std::is_same<typename Type<T1>::real_scalar, typename Type<T2>::real_scalar>::value == 1, "in and out data have different base data type.");
        typedef FFT_PLAN<typename Type<T1>::real_scalar> FFTP;
        typedef typename Type<T1>::vec vec_input;
        typedef typename Type<T1>::vec_c vec_F_side;

        int N1 = in.rows();
        int N2 = in.cols();
        typename Type<T1>::mat_c tmp(N1, N2);

        int k = 0;
        tmp.col(k) = FFTP::fft.fwd(vec_input(in.col(k)));
#pragma omp parallel for
        for (k = 1; k < N2; ++k)
        {
            tmp.col(k) = FFTP::fft.fwd(vec_input(in.col(k)));
        }
        out.transposeInPlace();
        tmp.transposeInPlace();

        k = 0;
        out.col(k) = FFTP::fft.fwd(vec_F_side(tmp.col(k)));
#pragma omp parallel for
        for (k = 1; k < N1; ++k)
        {
            out.col(k) = FFTP::fft.fwd(vec_F_side(tmp.col(k)));
        }
        out.transposeInPlace();
    }

    template <typename T>
    inline auto fft2(const Eigen::MatrixBase<T>& in)
    {
        typename Type<T>::mat_c out(in.rows(), in.cols());
        fft2_inplace(in, out);
        return out;
    }

    template <typename T1, typename T2>
    inline void ifft2_inplace(const Eigen::MatrixBase<T1>& in, Eigen::MatrixBase<T2>& out)
    {
        static_assert(std::is_same<typename Type<T1>::real_scalar, typename Type<T2>::real_scalar>::value == 1, "in and out data have different base data type.");
        typedef FFT_PLAN<typename Type<T1>::real_scalar> FFTP;
        typedef typename Type<T1>::vec vec_F_side;

        int N1 = in.rows();
        int N2 = in.cols();
        typename Type<T1>::mat_c tmp(N1, N2);

        int k = 0;
        tmp.col(k) = FFTP::fft.inv(vec_F_side(in.col(k)));
#pragma omp parallel for
        for (k = 1; k < N2; ++k)
        {
            tmp.col(k) = FFTP::fft.inv(vec_F_side(in.col(k)));
        }
        out.transposeInPlace();
        tmp.transposeInPlace();

        k = 0;
        out.col(k) = FFTP::fft.inv(vec_F_side(tmp.col(k)));
#pragma omp parallel for
        for (k = 1; k < N1; ++k)
        {
            out.col(k) = FFTP::fft.inv(vec_F_side(tmp.col(k)));
        }
        out.transposeInPlace();
    }

    template <typename T>
    inline auto ifft2_c2r(const Eigen::MatrixBase<T>& in)
    {
        mat out(in.rows(), in.cols());
        ifft2_inplace(in, out);
        return out;
    }

    template <typename T>
    inline auto ifft2_c2c(const Eigen::MatrixBase<T>& in)
    {
        mat_cd out(in.rows(), in.cols());
        ifft2_inplace(in, out);
        return out;
    } */
};