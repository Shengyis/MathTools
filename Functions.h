#pragma once
#include <tuple>
#include "types.h"

template <typename Tx, typename Ty>
auto meshgrid(const Eigen::MatrixBase<Tx>& x, const Eigen::MatrixBase<Ty>& y)
{
    typedef typename mathtool::Type<Tx>::mat Mx;
    typedef typename mathtool::Type<Ty>::mat My;
    int cols = x.size();
    int rows = y.size();
    Mx X(rows, cols);
    My Y(rows, cols);
    for (int k = 0; k < rows; ++k)
        X.row(k) = x.transpose();
    for (int k = 0; k < cols; ++k)
        Y.col(k) = y;
    return std::make_tuple(X, Y);
}

template <typename T>
auto meshgrid(const Eigen::MatrixBase<T>& x)
{
    return meshgrid(x, x);
}