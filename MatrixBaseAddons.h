using Base::Constant;

Derived &
operator+=(const Scalar &scalar)
{
    *this += Constant(rows(), cols(), scalar);
    return derived();
}

Derived &
operator-=(const Scalar &scalar)
{
    *this -= Constant(rows(), cols(), scalar);
    return derived();
}

template <typename OtherDerived>
Derived &
operator%=(const MatrixBase<OtherDerived> &mat)
{
    *this = *this % mat;
    return derived();
}

template <typename OtherDerived>
Derived &
operator/=(const MatrixBase<OtherDerived> &mat)
{
    static_assert(rows() == mat.rows() && cols() == mat.cols(), "both sides should have the same size.");
    *this = *this / mat;
    return derived();
}

friend const CwiseUnaryOp<internal::scalar_cos_op<Scalar>, const Derived>
cos(const MatrixBase<Derived> &mat)
{
    return  CwiseUnaryOp<internal::scalar_cos_op<Scalar>, const Derived>(mat.derived());
}

friend const CwiseUnaryOp<internal::scalar_sin_op<Scalar>, const Derived>
sin(const MatrixBase<Derived> &mat)
{
    return  CwiseUnaryOp<internal::scalar_sin_op<Scalar>, const Derived>(mat.derived());
}

friend const CwiseUnaryOp<internal::scalar_tan_op<Scalar>, const Derived>
tan(const MatrixBase<Derived> &mat)
{
    return  CwiseUnaryOp<internal::scalar_tan_op<Scalar>, const Derived>(mat.derived());
}

friend const CwiseBinaryOp<internal::scalar_sum_op<Scalar>, const Derived, const ConstantReturnType>
operator+(const MatrixBase<Derived> &mat, const Scalar &scalar)
{
    return CwiseBinaryOp<internal::scalar_sum_op<Scalar>, const Derived, const ConstantReturnType>(mat.derived(), Constant(mat.rows(), mat.cols(), scalar));
}

friend const CwiseBinaryOp<internal::scalar_sum_op<Scalar>, const ConstantReturnType, const Derived>
operator+(const Scalar &scalar, const MatrixBase<Derived> &mat)
{
    return CwiseBinaryOp<internal::scalar_sum_op<Scalar>, const ConstantReturnType, const Derived>(Constant(mat.rows(), mat.cols(), scalar), mat.derived());
}

friend const CwiseBinaryOp<internal::scalar_difference_op<Scalar>, const Derived, const ConstantReturnType>
operator-(const MatrixBase<Derived> &mat, const Scalar &scalar)
{
    return CwiseBinaryOp<internal::scalar_difference_op<Scalar>, const Derived, const ConstantReturnType>(mat.derived(), Constant(mat.rows(), mat.cols(), scalar));
}

friend const CwiseBinaryOp<internal::scalar_difference_op<Scalar>, const ConstantReturnType, const Derived>
operator-(const Scalar &scalar, const MatrixBase<Derived> &mat)
{
    return CwiseBinaryOp<internal::scalar_difference_op<Scalar>, const ConstantReturnType, const Derived>(Constant(mat.rows(), mat.cols(), scalar), mat.derived());
}

template <typename OtherDerived>
friend const CwiseBinaryOp<internal::scalar_product_op<Scalar>, const Derived, const OtherDerived>
operator%(const MatrixBase<Derived> &lhs, const MatrixBase<OtherDerived> &rhs)
{
    return CwiseBinaryOp<internal::scalar_product_op<Scalar>, const Derived, const OtherDerived>(lhs.derived(), rhs.derived());
}

friend const CwiseBinaryOp<internal::scalar_quotient_op<Scalar>, const ConstantReturnType, const Derived>
operator/(const Scalar &scalar, const MatrixBase<Derived> &mat)
{
    return CwiseBinaryOp<internal::scalar_quotient_op<Scalar>, const ConstantReturnType, const Derived>(Constant(mat.rows(), mat.cols(), scalar), mat.derived());
}

template <typename OtherDerived>
friend const CwiseBinaryOp<internal::scalar_quotient_op<Scalar>, const Derived, const OtherDerived>
operator/(const MatrixBase<Derived> &lhs, const MatrixBase<OtherDerived> &rhs)
{
    return CwiseBinaryOp<internal::scalar_quotient_op<Scalar>, const Derived, const OtherDerived>(lhs.derived(), rhs.derived());
}