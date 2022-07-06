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
        out.col(k) = fft(tmp.col(k));
#pragma omp parallel for
        for (k = 1; k < N1; ++k)
        {
            out.col(k) = fft(tmp.col(k));
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
        out.col(k) = ifft(tmp.col(k));
#pragma omp parallel for
        for (k = 1; k < N1; ++k)
        {
            out.col(k) = ifft(tmp.col(k));
        }
        out.transposeInPlace();
        return out;
    }
};