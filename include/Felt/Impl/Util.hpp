#ifndef INCLUDE_FELT_PUBLIC_UTIL_HPP_
#define INCLUDE_FELT_PUBLIC_UTIL_HPP_

#include <eigen3/Eigen/Dense>
#include <Felt/Impl/Common.hpp>

namespace Felt
{

/**
 * Get index in data array of position vector.
 *
 * The grid is packed in a 1D array, so this method is required to
 * get the index in that array of the D-dimensional position.
 *
 * @snippet test_Grid.cpp Position index
 *
 * @param pos_ position in grid.
 * @param size_ size of grid.
 * @param offset_ spatial offset of grid.
 * @return index in data array of pos in grid of given size and offset.
 */
template <Dim D>
PosIdx index (const VecDi<D>& pos_, const VecDi<D>& size_, const VecDi<D>& offset_)
{
	PosIdx pos_idx = 0;
	for (Dim i = 0; i < D; i++)
	{
		PosIdx pos_idx_axis = PosIdx(pos_(i) - offset_(i));

		for (Dim j = i+1; j < D; j++)
			pos_idx_axis *= PosIdx(size_(j));

		pos_idx += pos_idx_axis;
	}
	return pos_idx;
}

/**
 * Get position of index.
 *
 * Given an index and the dimensions and offset of a grid, calculate
 * the position vector that the index pertains to in a representative
 * 1D array.
 *
 * @snippet test_Grid.cpp Position index
 *
 * @param idx_ index in to query.
 * @param size_ size of grid.
 * @param offset_ spatial offset of grid.
 * @return position that the given index would represent in a grid of given size and offset.
 */
template <Dim D>
VecDi<D> index (
	PosIdx idx_, const VecDi<D>& size_, const VecDi<D>& offset_ = VecDi<D>::Zero()
) {
/*
Eg. 2D: row major order (3x4=12): (x,y)[idx] =>
(0,0)[0], (0,1)[1], (0,2)[2],  (0,3)[3]
(1,0)[4], (1,1)[5], (1,2)[6],  (1,3)[7]
(2,0)[8], (2,1)[9], (2,2)[10], (2,3)[11]

E.g. 3D:
z = idx % Dz
y = (idx/Dz) % Dy
x = (idx/Dz)/Dy % Dx
*/
	using AxisCoord = typename VecDi<D>::Scalar;
	VecDi<D> pos;

	// Note: since `Dim` is unsigned, we cannot allow `axis` to decrement below zero.
	for (Dim axis = D-1; axis != 0; axis--)
	{
		pos(axis) = AxisCoord(idx_ % PosIdx(size_(axis))) + offset_(axis);
		idx_ /= PosIdx(size_(axis));
	}
	pos(0) = AxisCoord(idx_ % PosIdx(size_(0))) + offset_(0);

	return pos;
}


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
		pos_rounded(dim) = INT(std::floor(pos_(dim)));
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
		pos_rounded(dim) = INT(std::ceil(pos_(dim)));
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
