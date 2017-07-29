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
		return VArrayData(pself->data().data(), Eigen::Index(pself->data().size()));
	}
};


template <class Derived>
class Spatial
{
private:
	using TraitsType = Impl::Traits<Derived>;
	using LeafType = typename TraitsType::LeafType;
	static constexpr Dim t_dims = TraitsType::t_dims;
	using VecDf = Felt::VecDf<t_dims>;
	using VecDi = Felt::VecDi<t_dims>;
	using VecDT = Felt::VecDT<LeafType, t_dims>;

	FLOAT m_dx;
	const VecDi m_pos_min;
	const VecDi m_pos_max;

protected:
	Spatial() : m_dx{1.0f}, m_pos_min{pself->offset()}, m_pos_max{pself->offset() + pself->size()}
	{

	}

	/**
	 * Mean curvature,
	 * \f$ \frac{1}{2} \nabla \bullet \frac{\nabla\phi}{\left|\nabla\phi\right|} \f$ .
	 *
	 * Based on difference of normals method.
	 *
	 * @param pos_ position in grid to query.
	 * @return curvature value
	 */
	template <typename PosType>
	LeafType curv (const Felt::VecDT<PosType, t_dims>& pos_) const
	{
		using VecDp = Felt::VecDT<PosType, t_dims>;

		const LeafType val_centre = pself->get(pos_);
		const VecDi& size = pself->size();
		VecDp dir(pos_);

		// Forward directed principal normal.
		VecDT n_forward;

		for (INT axis = 0; axis < size.size(); axis++)
		{
			dir(axis) += 1;

			const LeafType val_axis = pself->get(dir) - val_centre;
			LeafType val_neighs_sq = 0;

			// Loop other dimensions to get central difference across them.
			for (INT axis_neigh = 0; axis_neigh < size.size(); axis_neigh++)
			{
				// Only getting differences across other axes.
				if (axis_neigh != axis)
				{
					// Central difference across pself forward point.
					VecDp dir_neigh(dir);
					dir_neigh(axis_neigh) -= 1;
					const LeafType val_low = pself->get(dir_neigh);
					dir_neigh(axis_neigh) += 2;
					const LeafType val_high = pself->get(dir_neigh);

					const LeafType val_neigh = (val_high - val_low) / 2;
					val_neighs_sq += val_neigh*val_neigh;
				}
			}

			n_forward(axis) =		val_axis /
						sqrtf(val_axis*val_axis + val_neighs_sq);

			dir(axis) -= 1;
		}


		// Backward directed principal normal.
		VecDT n_backward;

		for (INT axis = 0; axis < size.size(); axis++)
		{
			dir(axis) -= 1;

			const LeafType val_axis = val_centre - pself->get(dir);
			LeafType val_neighs_sq = 0;

			// Loop other dimensions to get central difference across them.
			for (
				INT axis_neigh = 0; axis_neigh < size.size(); axis_neigh++
			) {
				// Only getting differences across other axes.
				if (axis_neigh != axis)
				{
					// Central difference across pself backward point.
					VecDp dir_neigh(dir);
					dir_neigh(axis_neigh) -= 1;
					const LeafType val_low = pself->get(dir_neigh);
					dir_neigh(axis_neigh) += 2;
					const LeafType val_high = pself->get(dir_neigh);

					const LeafType val_neigh = (val_high - val_low) / 2;
					val_neighs_sq += val_neigh*val_neigh;
				}
			}

			n_backward(axis) =		val_axis /
						sqrtf(val_axis*val_axis + val_neighs_sq);

			dir(axis) += 1;
		}

		const VecDT dn_by_dx = (n_forward - n_backward);

		LeafType curvature = dn_by_dx.sum() / 2;

		return curvature;
	}

	/**
	 * Calculate 2nd order divergence \f$ \nabla \bullet \nabla \phi \f$.
	 *
	 * @param pos_ position in grid to query.
	 * @return divergence value
	 */
	template <typename PosType>
	LeafType divergence (const Felt::VecDT<PosType, t_dims>& pos_) const
	{
		const VecDT vec_grad_f = gradF(pos_);
		const VecDT vec_grad_b = gradB(pos_);
		const VecDT vec_grad_diff = vec_grad_f - vec_grad_b;

		// Component-wise sum.
		const LeafType val = vec_grad_diff.sum();

		return val / (m_dx*m_dx);
	}

	/**
	 * Safe gradient, \f$ \nabla \phi \f$ .
	 *
	 * Will calculate central, forward or backward difference along each
	 * axis, depending what grid values are available.
	 * That is, for grid points at the edge of the grid it will
	 * return forward/backward differences.
	 *
	 * @param pos_ position in grid to query.
	 * @return vector with 2nd order if possible, 1st order if not, gradient.
	 */
	template <typename PosType>
	VecDT grad (const Felt::VecDT<PosType, t_dims>& pos_) const
	{
		using VecDR = Felt::VecDT<PosType, t_dims>;
		// Reference to GridBase dimensions.
		const VecDi& size = pself->size();
		// Vector to store gradient calculation.
		VecDT vec_grad;
		// Position for look-around.
		VecDR pos_test(pos_);

		// Central value.
		LeafType centre = pself->get(pos_);

		for (Dim axis = 0; axis < size.size(); axis++)
		{
			LeafType back = centre;
			LeafType forward = centre;
			Dim order = 0;
			// Check if backward value is within GridBase.
			pos_test(axis) -= 1;
			if (pself->inside(pos_test))
			{
				back = pself->get(pos_test);
				order++;
			}
			// Check if forward value is within GridBase.
			pos_test(axis) += 2;
			if (pself->inside(pos_test))
			{
				forward = pself->get(pos_test);
				order++;
			}
			pos_test(axis) -= 1;
			// Calculate central/forward/backward difference along pself
			// axis.
			if (order != 0)
				vec_grad(axis) = (forward - back) / LeafType(order);
			else
				vec_grad(axis) = 0;
		}

		return vec_grad / m_dx;
	}

	/**
	 * Entropy satisfying gradient, \f$ \nabla \phi \f$ .
	 *
	 * Use first order upwind scheme to select from forward or backward difference gradient along
	 * each cardinal direction.
	 *
	 * @param pos_ position in grid to query.
	 * @return vector of entropy satisfying gradient.
	 */
	template <typename PosType>
	VecDT gradE (const Felt::VecDT<PosType, t_dims>& pos_) const
	{
		using VecDp = Felt::VecDT<PosType, t_dims>;
		// Value at pself point.
		const LeafType centre = pself->get(pos_);
		// Reference to grid dimensions.
		const VecDi& size = pself->size();
		// Vector to store gradient calculation.
		VecDT vec_grad;
		// Position for look-around.
		VecDp pos_test(pos_);

		for (INT axis = 0; axis < size.size(); axis++)
		{
			pos_test(axis) -= 1;
			LeafType back = pself->get(pos_test);
			pos_test(axis) += 2;
			LeafType forward = pself->get(pos_test);
			pos_test(axis) -= 1;

			back = std::max(centre - back, 0.0f);
			forward = std::min(forward - centre, 0.0f);

			vec_grad(axis) = forward + back;
		}

		return vec_grad / m_dx;
	}

	/**
	 * Forward difference gradient.
	 *
	 * @param pos_ position in grid to query.
	 * @return vector with forward difference gradient.
	 */
	template <typename PosType>
	VecDT gradF (const Felt::VecDT<PosType, t_dims>& pos_) const
	{
		// Value at this point.
		const LeafType centre = pself->get(pos_);
		// Reference to GridBase dimensions.
		const VecDi& size = pself->size();
		// Vector to store gradient calculation.
		VecDT vec_grad;
		// Position for look-ahead.
		Felt::VecDT<PosType, t_dims> pos_neigh(pos_);

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
	VecDT gradB (const Felt::VecDT<PosType, t_dims>& pos_) const
	{
		// pself->getue at pself point.
		const LeafType centre = pself->get(pos_);
		// Reference to GridBase dimensions.
		const VecDi& size = pself->size();
		// Vector to store gradient calculation.
		VecDT vec_grad;
		// Position for look-behind.
		Felt::VecDT<PosType, t_dims> vec_dir(pos_);

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
	VecDT gradC (const Felt::VecDT<PosType, t_dims>& pos_) const
	{
		// Reference to GridBase dimensions.
		const VecDi& size = pself->size();
		// Vector to store gradient calculation.
		VecDT vec_grad;
		// Position for look-around.
		Felt::VecDT<PosType, t_dims> vec_dir(pos_);

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
	FLOAT dx () const
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
		const VecDi& pos_floor = pos_.array().floor().matrix().template cast<INT>();

		// Store all 2^d corners.
		DataArray< LeafType > val_corners(PosIdx(1 << size.size()));

		// Get all corners of containing cell.
		for (ListIdx i = 0; i < val_corners.size(); i++)
		{
			// 0 = 00 => (x,y)
			// 1 = 01 => (x+1,y)
			// 2 = 10 => (x,y+1)
			// 3 = 11 => (x+1,y+1)

			VecDi pos_corner(pos_floor);
			for (Dim axis = 0; axis < pos_corner.size(); axis++)
			{
				const INT dir = (i >> axis) & 1;
				pos_corner(axis) += dir;
			}

			val_corners[i] = pself->get(pos_corner);
		}

		// Translate position vector into 'hypercube space',
		// so 0 <= v(x) <= 1.
		VecDf pos_centred = pos_ - pos_.array().floor().matrix();

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
	static void interp (DataArray<LeafType>& val_corners_, const VecDf& pos_)
	{
		const ListIdx num_corners = val_corners_.size();

		// Number of values returned.
		// This is a power of 2 less than input dimensions
		// (cube becomes square, square becomes line, line becomes point).
		const ListIdx num_out = num_corners >> 1;

		static_assert(
			std::is_same<ListIdx, unsigned long>::value,
			"num_corners must be unsigned to use __builtin_ctzl" //__builtin_ffsl
		);
		// The axis along which to interpolate.
		// This is computed from the dimensions of the original input and
		// the dimensions of the intended output.
		const Dim axis_idx = pos_.size() - Dim(__builtin_ctzl(num_corners));

		// The weighting to be used in interpolating each pair of points.
		// This is the position along the axis of interpolation.
		const FLOAT axis_pos = pos_(axis_idx);

		for (ListIdx i = 0; i < num_out; i++)
		{
			const LeafType low = val_corners_[(i << 1)];
			const LeafType high = val_corners_[(i << 1) + 1];
			const LeafType val = axis_pos*high + (1.0f-axis_pos)*low;
			val_corners_[i] = val;
		}
		val_corners_.resize(num_out);
	}
	
	/**
	 * Call a lambda passing neighbours of a position in the cardinal directions.
	 *
	 * @param pos_ position to search around
	 * @param fn_ lambda function
	 */
	template <typename Fn>
	static void neighs (VecDi pos_, Fn fn_)
	{
		for (Dim axis = 0; axis < t_dims; axis++)
		{
			pos_(axis) -= 1;
			fn_(pos_);
			pos_(axis) += 2;
			fn_(pos_);
			pos_(axis) -= 1;
		}
	}
};

} // Numeric.
} // Mixin.
} // Impl.
} // Felt.


#endif /* INCLUDE_FELT_IMPL_MIXIN_NUMERICMIXIN_HPP_ */
