#ifndef INCLUDE_FELT_PUBLIC_UTIL_HPP_
#define INCLUDE_FELT_PUBLIC_UTIL_HPP_

#include <eigen3/Eigen/Dense>
#include <Felt/Impl/Common.hpp>

namespace Felt
{
/**
 * String format a vector (useful for logging).
 *
 * @param vec_ vector to stringify.
 * @return formatted vector.
 */
template<class VecType>
std::string format(const VecType& vec_)
{
	using namespace Eigen;
	IOFormat fmt(
		StreamPrecision, DontAlignCols, ", ", ", ",
		"", "", "(", ")"
	);
	std::stringstream str;
	str << vec_.format(fmt);
	return str.str();
}

/**
 * Get the sign of a value (+/-1).
 *
 * @param val_ value to get signum for.
 * @return -1 for negative, +1 for positive.
 */
template <typename T> INT sgn(T val_)
{
	return (T(0) < val_) - (val_ < T(0));
}

/**
 * Round float accuracy position to integer accuracy.
 *
 * @tparam D the dimension of the vector.
 * @param pos_ float vector to round
 * @return rounded integer vector (away from zero).
 */
template <INT D>
VecDi<D> round(const VecDf<D>& pos_)
{
	VecDi<D> pos_rounded;
	for (UINT dim = 0; dim < pos_.size(); dim++)
		pos_rounded(dim) = (INT)(pos_(dim) + sgn(pos_(dim)) * 0.5f);
	return pos_rounded;
}

/**
 * Call std::floor on each element of float vector to give integer vector.
 *
 * @tparam D the dimension of the vector.
 * @param pos_ float vector
 * @return floored integer vector (away from zero).
 */
template <INT D>
VecDi<D> floor(const VecDf<D>& pos_)
{
	VecDi<D> pos_rounded;
	for (UINT dim = 0; dim < pos_.size(); dim++)
		pos_rounded(dim) = std::floor(pos_(dim));
	return pos_rounded;
}

/**
 * Call std::ceil on each element of float vector to give integer vector.
 *
 * @tparam D the dimension of the vector.
 * @param pos_ float vector
 * @return ceiling of integer vector (away from zero).
 */
template <INT D>
VecDi<D> ceil(const VecDf<D>& pos_)
{
	VecDi<D> pos_rounded;
	for (UINT dim = 0; dim < pos_.size(); dim++)
		pos_rounded(dim) = std::ceil(pos_(dim));
	return pos_rounded;
}

/**
 * Call std::floor on each element of float vector to give float vector.
 *
 * @tparam D the dimension of the vector.
 * @param pos_ float vector
 * @return floored float vector (away from zero).
 */
template <INT D>
VecDf<D> floorf(const VecDf<D>& pos_)
{
	VecDf<D> pos_rounded;
	for (UINT dim = 0; dim < pos_.size(); dim++)
		pos_rounded(dim) = std::floor(pos_(dim));
	return pos_rounded;
}
} // Felt.

#endif /* INCLUDE_FELT_PUBLIC_UTIL_HPP_ */
