#ifndef INCLUDE_FELT_UTIL_HPP_
#define INCLUDE_FELT_FELT_HPP_

#include <eigen3/Eigen/Dense>

namespace felt
{

/// Shorthand cast to CRTP derived type.
#define self static_cast<DerivedType*>(this)
/// Shorthand cast to const CRTP derived type.
#define cself static_cast<const DerivedType*>(this)

/**
 * Use 32 bit float by default.
 */
using FLOAT = float;

/**
 * Use 32 bit int by default.
 */
using INT = int;

/**
 * Use 32 bit unsigned int by default.
 */
using UINT = unsigned;

/**
 * Shorthand for D-dimensional vector with elements of T type.
 */
template <typename T, UINT D>
using VecDT = Eigen::Matrix<T, D, 1>;
/**
 * Shorthand for D-dimensional float vector.
 */
template <UINT D>
using VecDf = VecDT<FLOAT, D>;
/**
 * Shorthand for D-dimensional integer vector.
 */
template <UINT D>
using VecDi = VecDT<INT, D>;
/**
 * Shorthand for D-dimensional unsigned integer vector.
 */
template <UINT D>
using VecDu = VecDT<UINT, D>;

/**
 * Shorthand for 2D float vector.
 */
using Vec2f = VecDf<2>;
/**
 * Shorthand for 2D unsigned integer vector.
 */
using Vec2u = VecDu<2>;
/**
 * Shorthand for 2D integer vector.
 */
using Vec2i = VecDi<2>;
/**
 * Shorthand for 3D float vector.
 */
using Vec3f = VecDf<3>;
/**
 * Shorthand for 3D unsigned integer vector.
 */
using Vec3u = VecDu<3>;
/**
 * Shorthand for 3D integer vector.
 */
using Vec3i = VecDi<3>;

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

} // End namespace felt.
#endif /* INCLUDE_FELT_UTIL_HPP_ */
