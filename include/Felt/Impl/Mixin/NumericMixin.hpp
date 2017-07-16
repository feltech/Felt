#ifndef INCLUDE_FELT_IMPL_MIXIN_NUMERICMIXIN_HPP_
#define INCLUDE_FELT_IMPL_MIXIN_NUMERICMIXIN_HPP_

#include <Felt/Impl/Common.hpp>
namespace Felt
{
namespace Impl
{
namespace Mixin
{
namespace Numeric
{

template <class Derived>
class Snapshot
{
private:
	/// Type of data to store in grid nodes.
	using LeafType = typename Traits<Derived>::LeafType;
protected:
	/// Map of of POD to Eigen::Array for manipulation using Eigen BLAS methods.
	using VArrayData = Eigen::Map< Eigen::Array<LeafType, 1, Eigen::Dynamic> >;

	/**
	 * Map the raw data to an Eigen::Map, which can be used for BLAS arithmetic.
	 *
	 * @return Eigen compatible vector of data array.
	 */
	VArrayData vdata()
	{
		return VArrayData(pself->data().data(), pself->data().size());
	}
};


template <class Derived>
class Spatial
{
private:
	using TraitsType = Impl::Traits<Derived>;
	using LeafType = typename TraitsType::LeafType;
	static constexpr UINT Dims = TraitsType::Dims;
	using VecDf = Felt::VecDf<Dims>;
	using VecDi = Felt::VecDi<Dims>;
	using VecDT = Felt::VecDT<LeafType, Dims>;

	FLOAT m_dx;
	const VecDi m_pos_min;
	const VecDi m_pos_max;

protected:
	Spatial() : m_dx{1.0f}, m_pos_min{pself->offset()}, m_pos_max{pself->offset() + pself->size()}
	{

	}

	/**
	 * Forward difference gradient.
	 *
	 * @param pos_ position in grid to query.
	 * @return vector with forward difference gradient.
	 */
	template <typename PosType>
	VecDT gradF (const Felt::VecDT<PosType, Dims>& pos_) const
	{
		// Value at this point.
		const LeafType centre = pself->get(pos_);
		// Reference to GridBase dimensions.
		const VecDi& size = pself->size();
		// Vector to store gradient calculation.
		VecDT vec_grad;
		// Position for look-ahead.
		Felt::VecDT<PosType, Dims> pos_neigh(pos_);

		for (INT axis = 0; axis < size.size(); axis++)
		{
			pos_neigh(axis) += 1;
			vec_grad(axis) = pself->get(pos_neigh) - centre;
			pos_neigh(axis) -= 1;
		}

		return vec_grad / m_dx;
	}

	/**
	 * Backward difference gradient, \f$ \nabla \phi \f$ .
	 *
	 * @param pos_ position in grid to query.
	 * @return vector with backward difference gradient.
	 * @return
	 */
	template <typename PosType>
	VecDT gradB (const Felt::VecDT<PosType, Dims>& pos_) const
	{
		// pself->getue at pself point.
		const LeafType centre = pself->get(pos_);
		// Reference to GridBase dimensions.
		const VecDi& size = pself->size();
		// Vector to store gradient calculation.
		VecDT vec_grad;
		// Position for look-behind.
		Felt::VecDT<PosType, Dims> vec_dir(pos_);

		for (INT axis = 0; axis < size.size(); axis++) {
			vec_dir(axis) -= 1;
			vec_grad(axis) = centre - pself->get(vec_dir);
			vec_dir(axis) += 1;
		}

		return vec_grad / m_dx;
	}

	/**
	 * Central difference gradient, \f$ \nabla \phi \f$ .
	 *
	 * @param pos_ position in grid to query.
	 * @return vector with central difference gradient.
	 */
	template <typename PosType>
	VecDT gradC (const Felt::VecDT<PosType, Dims>& pos_) const
	{
		// Reference to GridBase dimensions.
		const VecDi& size = pself->size();
		// Vector to store gradient calculation.
		VecDT vec_grad;
		// Position for look-around.
		Felt::VecDT<PosType, Dims> vec_dir(pos_);

		for (INT axis = 0; axis < size.size(); axis++) {
			vec_dir(axis) -= 1;
			const LeafType back = pself->get(vec_dir);
			vec_dir(axis) += 2;
			const LeafType forward = pself->get(vec_dir);
			vec_dir(axis) -= 1;

			vec_grad(axis) = (forward - back) /  2;
		}

		return vec_grad / m_dx;
	}

	/**
	 * Get grid's delta x, \f$ \Delta x \f$ .
	 *
	 * @snippet test_Grid.cpp Delta x setter
	 *
	 * @return representative spatial size of a leaf node.
	 */
	const FLOAT dx () const
	{
		return m_dx;
	}

	/**
	 * Set grid's delta x, \f$ \Delta x \f$.
	 *
	 * @snippet test_Grid.cpp Delta x setter
	 *
	 * @param dx_ the new representative spatial size of a leaf node.
	 */
	void dx (const FLOAT dx_)
	{
		m_dx = dx_;
	}

	/**
	 * Get interpolated grid value.
	 *
	 * Passing a floating point position vector will initiate a linear interpolation of the grid at
	 * that real-valued location.
	 *
	 * @param pos_ position in grid to query.
	 * @return linearly interpolated value at given position.
	 */
	LeafType get (const VecDf& pos_) const
	{
		return interp(pos_);
	}

	/**
	 * Linear interpolation.
	 *
	 * @snippet test_Grid.cpp Interpolation

	 * @param pos_ position in grid to query.
	 * @return interpolated value at given position.
	 */
	LeafType interp (const VecDf& pos_) const
	{
		const VecDi& size = pself->size();
		const VecDi& pos_floor = Felt::floor(pos_);

		// Store all 2^d corners.
		std::vector< LeafType > val_corners(1 << size.size());

		// Get all corners of containing cell.
		for (UINT i = 0; i < val_corners.size(); i++)
		{
			// 0 = 00 => (x,y)
			// 1 = 01 => (x+1,y)
			// 2 = 10 => (x,y+1)
			// 3 = 11 => (x+1,y+1)

			VecDi pos_corner(pos_floor);
			for (INT axis = 0; axis < pos_corner.size(); axis++)
			{
				const INT dir = (i >> axis) & 1;
				pos_corner(axis) += dir;
			}

			val_corners[i] = pself->get(pos_corner);
		}

		// Translate position vector into 'hypercube space',
		// so 0 <= v(x) <= 1.
		VecDf pos_centred = pos_ - Felt::floorf(pos_);

		// Repeatedly reduce along axes,
		// i.e. hypercube -> cube -> square -> line -> point
		while (val_corners.size() > 1)
			interp(val_corners, pos_centred);

		return val_corners[0];
	}

	/**
	 * Interpolate down one dimension.
	 *
	 * The values of val_corners_ are interpolated to one dimension smaller than they are
	 * currently (cube->square, square->line, line->point).
	 *
	 * @param val_corners_ list of corner values.
	 * @param pos_ real-valued position to interpolate to.
	 * @return list of interpolated values
	 */
	static void interp (std::vector<LeafType>& val_corners_, const VecDf& pos_)
	{
		const UINT num_corners = val_corners_.size();

		// Number of values returned.
		// This is a power of 2 less than input dimensions
		// (cube becomes square, square becomes line, line becomes point).
		const UINT num_out = num_corners >> 1;

		// The axis along which to interpolate.
		// This is computed from the dimensions of the original input and
		// the dimensions of the intended output.
		const UINT axis_idx = pos_.size() - log2(num_corners);

		// The weighting to be used in interpolating each pair of points.
		// This is the position along the axis of interpolation.
		const FLOAT axis_pos = pos_(axis_idx);

		for (UINT i = 0; i < num_out; i++)
		{
			const LeafType low = val_corners_[(i << 1)];
			const LeafType high = val_corners_[(i << 1) + 1];
			const LeafType val = axis_pos*high + (1.0f-axis_pos)*low;
			val_corners_[i] = val;
		}
		val_corners_.resize(num_out);
	}
};

} // Numeric.
} // Mixin.
} // Impl.
} // Felt.


#endif /* INCLUDE_FELT_IMPL_MIXIN_NUMERICMIXIN_HPP_ */
