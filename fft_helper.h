#pragma once
#include "types.h"
#include <eigen3/unsupported/Eigen/FFT>

namespace mathtool 
{
    template <typename T>
    struct FFT_PLAN
    {
        static Eigen::FFT<typename Type<T>::real_scalar> fft;
    };

    template <typename T>
    Eigen::FFT<typename Type<T>::real_scalar> FFT_PLAN<T>::fft;
};