#ifndef Grid_hpp
#define Grid_hpp

#include <inttypes.h>
#include <cmath>
#include <vector>
#include <functional>
#include <eigen3/Eigen/Dense>
#include <boost/iterator/iterator_facade.hpp>

namespace felt
{
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

/** @defgroup Grids
 *
 *  Base arbitrarily dimensioned grid classes storing arbitrary data types.
 *  @{
 */

/**
 * Traits for classes CRTP derived from GridBase
 */
template <class Derived> struct GridBaseTraits {};


/**
 * Abstract base class for n-dimensional grid.
 *
 * Subclasses can override the get() method to return a datatype other
 * than that stored in the grid, so that grid values can be mutated
 * before e.g. using in gradient calculations.
 *
 * @tparam Derived the CRTP derived class
 */
template <class Derived>
class GridBase
{
public:
	/// This GridBase type.
	using ThisType = GridBase<Derived>;
	/// CRTP derived class.
	using DerivedType = typename GridBaseTraits<Derived>::ThisType;
	/// Dimension of the grid.
	static const UINT Dims = GridBaseTraits<Derived>::Dims;
	/// Type of data to store in grid nodes.
	using LeafType = typename GridBaseTraits<Derived>::LeafType;
	/// Type of data to return when grid nodes are queried (usually same as LeafType).
	using RetType = typename GridBaseTraits<Derived>::RetType;

	/**
	 * D-dimensional unsigned integer vector.
	 */
	using VecDu = felt::VecDu<Dims>;
	/**
	 * D-dimensional integer vector.
	 */
	using VecDi = felt::VecDi<Dims>;
	/**
	 * D-dimensional float vector.
	 */
	using VecDf = felt::VecDf<Dims>;
	/**
	 * D-dimensional vector of type LeafType.
	 */
	using VecDT = felt::VecDT<LeafType, Dims>;
	/**
	 * Dynamic 1D vector (a resizeable array of data) for storage of grid data.
	 */
	using ArrayData = Eigen::Array<LeafType, 1, Eigen::Dynamic>;
	/**
	 * Resizeable array of VecDi (grid locations).
	 *
	 * Uses Eigen::aligned_allocator for optimal address alignment.
	 */
	using PosArray = std::vector<VecDi, Eigen::aligned_allocator<VecDi> >;


	/**
	 * Iterator for contiguous cycling over entire grid.
	 */
	class iterator : public boost::iterator_facade<
		GridBase::iterator, const VecDi&, boost::forward_traversal_tag
	> {
	private:
		friend class boost::iterator_core_access;
		/// Index in grid data array currently pointed to.
		UINT			m_idx;
		/// Position in grid currently pointed to.
		VecDi			m_pos;
		/// Size of the grid.
		const VecDu&	m_dims;
		/// Offset of grid from bottom-front-left to (0,0,0).
		const VecDi&	m_offset;
	public:
		/**
		 * Default constructor for empty iterator.
		 */
		iterator() : m_idx(0), m_dims(VecDu::Constant(0)), m_offset(VecDi::Constant(0))
		{}

		/**
		 * Construct an iterator from given grid beginning at startIdx in the data array.
		 * @param grid_ grid to iterator over.
		 * @param start_idx_ index in underlying data to start at.
		 */
		iterator(const GridBase& grid_, const UINT& start_idx_ = 0)
		: m_idx(start_idx_), m_dims(grid_.dims()), m_offset(grid_.offset()),
		  m_pos(ThisType::index(start_idx_, grid_.dims(), grid_.offset()))
		{}

	private:

		/**
		 * Increment the iterator to next position in the grid data.
		 */
		void increment()
		{
			m_idx++;
			m_pos = ThisType::index(m_idx, m_dims, m_offset);
		}

		/**
		 * Check for equality with another grid data iterator.
		 *
		 * @param other_ iterator to compare with.
		 * @return true if equal, false if not equal.
		 */
		bool equal(iterator const& other_) const
		{
			return this->m_idx == other_.m_idx;
		}

		/**
		 * Get the position currently pointed to by the iterator.
		 *
		 * @return position currently represented by iterator.
		 */
		const VecDi& dereference() const
		{
			return m_pos;
		}
	};


protected:
	/**
	 * The translational offset of the grid's zero coordinate.
	 */
	VecDi m_offset;
	/**
	 * The dimensions (size) of the grid.
	 */
	VecDu m_dims;

	/// Minimum position stored in grid (equal to m_offset).
	VecDi m_pos_min;
	/// One more than maximum position stored in grid (equal to m_offset + m_dims).
	VecDi m_pos_max;

	/**
	 * The physical size of a grid node \f$\Delta x\f$ (used for spatial derivatives).
	 */
	FLOAT m_dx;
	/**
	 * The actual grid data store.
	 */
	ArrayData m_data;

public:
	/**
	 * Trivial destructor.
	 */
	~GridBase ()
	{}

	/**
	 * Initialise a zero-size grid.
	 *
	 * If custom initialisation is required, then use this constructor
	 * and call init() in the derived class.
	 */
	GridBase () :
	m_offset(VecDi::Zero()),
	m_dims(VecDu::Zero()),
	m_dx(1)
	{
	}

	/**
	 * Initialise a grid with given dimension, offset and delta x.
	 *
	 * @param size_ size of grid.
	 * @param offset_ spatial offset of grid.
	 * @param delta_ spatial size of a leaf node, for spatial derivative calculation.
	 */
	GridBase (
		const VecDu& size_, const VecDi& offset_ = VecDi::Zero(),
		const FLOAT& delta_ = 1
	) :
	m_dx(1)
	{
		this->init(size_, offset_, delta_);
	}

	/**
	 * Initialise the grid dimensions and offset.
	 *
	 * @param size_ size of grid.
	 * @param offset_ spatial offset of grid.
	 * @param delta_ representative spatial size of a leaf node, for spatial derivative calculation.
	 */
	void init (
		const VecDu& size_, const VecDi& offset_ = VecDi::Zero(),
		const FLOAT& delta_ = 1
	) {
		DerivedType* self = static_cast<DerivedType*>(this);

		self->dx(delta_);
		self->dims(size_);
		self->offset(offset_);
	}

	/**
	 * Set grid offset.
	 *
	 * The offset is used to 'centre' the grid, so that e.g. negative
	 * grid positions can be used. It is equal to the spatial position
	 * of the zero coordinate.
	 *
	 * @param offset_ spatial offset of grid.
	 */
	void offset (const VecDi& offset_)
	{
		m_offset = offset_;
		m_pos_min = offset_;
		m_pos_max = m_offset + m_dims.template cast<INT>();
	}

	/**
	 * Get the grid offset parameter.
	 *
	 * @return spatial offset of grid.
	 */
	const VecDi& offset () const
	{
		return m_offset;
	}

	/**
	 * Get grid's delta x, \f$ \Delta x \f$ .
	 *
	 * @return representative spatial size of a leaf node.
	 */
	inline const FLOAT& dx () const
	{
		return m_dx;
	}

	/**
	 * Set grid's delta x, \f$ \Delta x \f$ .
	 *
	 * @param dx_ the new representative spatial size of a leaf node.
	 */
	void dx (const FLOAT& dx_)
	{
		m_dx = dx_;
	}

	/**
	 * Shorthand to access grid values.
	 *
	 * @param pos_ position in grid to query.
	 * @return value represented at given position.
	 */
	RetType& operator() (const VecDi& pos_)
	{
		return this->get(pos_);
	}

	/**
	 * Shorthand to access grid values (const version).
	 *
	 * @param pos_ position in grid to query.
	 * @return value represented at given position.
	 */
	const RetType& operator() (const VecDi& pos_) const
	{
		return this->get(pos_);
	}

	/**
	 * Get interpolated grid value.
	 *
	 * Passing a floating point position vector will initiate a linear
	 * interpolation of the grid at that real-valued location.
	 *
	 * @param pos_ position in grid to query.
	 * @return linearly interpolated value at given position.
	 */
	const LeafType operator() (const VecDf& pos_) const
	{
		return this->interp(pos_);
	}

	/**
	 * Get interpolated grid value.
	 *
	 * Convenience operator to return linearly interpolated value given a real-valued location.
	 *
	 * @param pos_ position in grid to query.
	 * @return linearly interpolated value at given position.
	 */
	LeafType operator() (const VecDf& pos_)
	{
		return this->interp(pos_);
	}

	/**
	 * Abstract getter to be overriden by subclasses.
	 *
	 * This allows the subclass to mutate the value stored in the grid before it is used.
	 *
	 * @param pos_ position in grid to query.
	 * @return value represented at given position.
	 */
	RetType& get (const VecDi& pos_)
	{
		return static_cast<DerivedType*>(this)->get(pos_);
	}

	/**
	 * Abstract getter to be overriden by subclasses (const version).
	 *
	 * @param pos_ position in grid to query.
	 * @return value represented at given position.
	 */
	const RetType& get (const VecDi& pos_) const
	{
		return static_cast<const DerivedType*>(this)->get(pos_);
	}

	/**
	 * Get interpolated grid value.
	 *
	 * Passing a floating point position vector will initiate a linear
	 * interpolation of the grid at that real-valued location.
	 *
	 * @param pos_ position in grid to query.
	 * @return linearly interpolated value at given position.
	 */
	const LeafType get (const VecDf& pos_) const
	{
		return this->interp(pos_);
	}

	/**
	 * Get interpolated grid value.
	 *
	 * Passing a floating point position vector will initiate a linear
	 * interpolation of the grid at that real-valued location.
	 *
	 * @param pos_ position in grid to query.
	 * @return linearly interpolated value at given position.
	 */
	LeafType get (const VecDf& pos_)
	{
		return this->interp(pos_);
	}

	/**
	 * Get index in data array of position vector.
	 *
	 * The grid is packed in a 1D array, so this method is required to
	 * get the index in that array of the D-dimensional position.
	 *
	 * @param pos_ position in grid to query.
	 * @return index in internal data array of this grid position.
	 */
	UINT index (const VecDi& pos_) const
	{
		return ThisType::index(pos_, this->dims(), this->offset());
	}

	/**
	 * Get index in data array of position vector.
	 *
	 * The grid is packed in a 1D array, so this method is required to
	 * get the index in that array of the D-dimensional position.
	 *
	 * @param pos_ position in grid.
	 * @param size_ size of grid.
	 * @param offset_ spatial offset of grid.
	 * @return index in data array of pos in grid of given size and offset.
	 */
	static UINT index (
		const VecDi& pos_, const VecDu& size_, const VecDi& offset_ = VecDi::Constant(0)
	) {
		UINT idx = 0;
		for (INT i = 0; i < size_.size(); i++)
		{
			INT u_pos = pos_(i) - offset_(i);
			for (INT j = i+1; j < size_.size(); j++)
			{
				u_pos *= size_(j);
			}
			idx += u_pos;
		}
		return idx;
	}

	/**
	 * Get position of index.
	 *
	 * Given an index in the 1D grid data array, calculate the position vector that it pertains to.
	 *
	 * @param idx_ index in internal data array to query.
	 * @return the position in the grid represented in the data array at given index.
	 */
	VecDi index (const UINT& idx_) const
	{
		return ThisType::index(idx_, this->dims(), this->offset());
	}

	/**
	 * Get position of index.
	 *
	 * Given an index and the dimensions and offset of a grid, calculate
	 * the position vector that the index pertains to in a representative
	 * 1D array.
	 *
	 * @param idx_ index in to query.
	 * @param size_ size of grid.
	 * @param offset_ spatial offset of grid.
	 * @return position that the given index would represent in a grid of given size and offset.
	 */
	static VecDi index (UINT idx_, const VecDu& size_, const VecDi& offset_ = VecDi::Zero())
	{
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
		VecDi pos;

		for (INT axis = size_.size()-1; axis >= 0; axis--)
		{
			pos(axis) = idx_ % size_(axis) + offset_(axis);
			idx_ /= size_(axis);
		}

		return pos;
	}

	/**
	 * Retrieve a reference to the raw grid data array.
	 *
	 * @return the array of data representing nodes in the grid.
	 */
	ArrayData& data ()
	{
		return m_data;
	}

	/**
	 * Retrieve a reference to the raw grid data array (const version).
	 *
	 * @return the array of data representing nodes in the grid.
	 */
	const ArrayData& data () const
	{
		return m_data;
	}

	/**
	 * Get an iterator wrapping a position at 0 in data index.
	 *
	 * @return an iterator referencing the first element in the internal data array.
	 */
	const iterator begin () const
	{
		return iterator(*this, 0);
	}

	/**
	 * Get an iterator wrapping a position at one greater than the size of the data.
	 *
	 * @return an iterator pointing off the end of the internal data array.
	 */
	const iterator end () const
	{
		return iterator(*this, this->data().size());
	}

	/**
	 * Set the dimensions of the grid and resize it.
	 *
	 * Values will be default initialised (or not at all).
	 *
	 * @param size_ new size of the grid.
	 * @return
	 */
	void dims (const VecDu& size_)
	{
		m_dims = size_;

		INT arr_size = m_dims(0);
		for (INT i = 1; i < m_dims.size(); i++)
		{
			arr_size *= m_dims(i);
		}
		m_data.resize(arr_size);
	}

	/**
	 * Get grid size.
	 *
	 * @return the size of the grid in grid nodes.
	 */
	const VecDu& dims () const
	{
		return m_dims;
	}

	/**
	 * Fill grid with a single value.
	 *
	 * @param val_ value to fill the grid with.
	 */
	void fill (const LeafType& val_)
	{
		this->data().fill(val_);
	}

	/**
	 * Test if a position is inside the grid bounds.
	 *
	 * @tparam PosType the type of position vector (i.e. float vs. int).
	 * @param pos_ position in grid to query.
	 * @return true if position lies inside the grid, false otherwise.
	 */
	template <typename PosType>
	bool inside (const felt::VecDT<PosType, Dims>& pos_) const
	{
		return inside(
			pos_, this->offset(),
			this->offset() + this->dims().template cast<INT>()
		);
	}

	/**
	 * Test if a position is inside given bounds.
	 *
	 * @tparam PosType the type of position vector (i.e. float vs. int).
	 * @param pos_ position in grid to query.
	 * @param pos_min_ minimum allowed position.
	 * @param pos_max_ one more than the maximum allowed position.
	 * @return true if position lies inside the grid, false otherwise.
	 */
	template <typename PosType>
	static bool inside (
		const felt::VecDT<PosType, Dims>& pos_, const VecDi& pos_min_, const VecDi& pos_max_
	) {

		for (INT i = 0; i < pos_.size(); i++)
		{
			if (pos_(i) >= (PosType)pos_max_(i))
				return false;
			if (pos_(i) < (PosType)pos_min_(i))
				return false;
		}
		return true;
	}

	/**
	 * Get the neighbouring positions in the cardinal directions.
	 *
	 * Neighbour positions will be added to aout_.
	 *
	 * If bcheck is true then duplicates will not be allowed in vout,
	 * using a linear search to ensure so.
	 *
	 * @param pos_ position in grid to query.
	 * @param aout_ list to add found neighbouring positions to.
	 * @param bcheck_ set to true to ensure no duplicates in aout_.
	 */
	void neighs (const VecDi& pos_, PosArray& aout_, bool bcheck = false) const
	{
		// Reference to GridBase dimensions.
		const VecDu& dims = this->dims();
		// Most likely all neighbours are valid.
		aout_.reserve(2*Dims);
		// Position for look-around.
		VecDi vec_dir(pos_);
		for (INT axis = 0; axis < dims.size(); axis++) {
			// Check if backward value is within GridBase.
			vec_dir(axis) -= 1;
			if (this->inside(vec_dir) && (!bcheck ||
				std::find(aout_.begin(), aout_.end(), vec_dir) == aout_.end())
			)
			{
				aout_.push_back(vec_dir);
			}
			// Check if forward value is within GridBase.
			vec_dir(axis) += 2;
			if (this->inside(vec_dir) && (!bcheck ||
				std::find(aout_.begin(), aout_.end(), vec_dir) == aout_.end())
			)
			{
				aout_.push_back(vec_dir);
			}
			vec_dir(axis) -= 1;
		}
	}

	/**
	 * Call a lambda passing neighbours of a position in the cardinal directions.
	 *
	 * @param pos_ position to search around
	 * @param fn_ lambda function
	 */
	void neighs (const VecDi& pos_, const std::function<void(const VecDi&)>& fn_) const
	{
		// Position for look-around.
		VecDi vec_dir(pos_);
		for (INT axis = 0; axis < this->dims().size(); axis++)
		{
			// Check if backward value is within GridBase.
			vec_dir(axis) -= 1;
			if (this->inside(vec_dir))
				fn_(vec_dir);

			// Check if forward value is within GridBase.
			vec_dir(axis) += 2;

			if (this->inside(vec_dir))
				fn_(vec_dir);

			vec_dir(axis) -= 1;
		}
	}

	/**
	 * Forward difference gradient.
	 *
	 * @param pos_ position in grid to query.
	 * @return vector with forward difference gradient.
	 */
	template <typename PosType>
	VecDT gradF (const felt::VecDT<PosType, Dims>& pos_) const
	{
		// Value at this point.
		const LeafType val_centre = (*this)(pos_);
		// Reference to GridBase dimensions.
		const VecDu& dims = this->dims();
		// Vector to store gradient calculation.
		VecDT vec_grad(dims.size());
		// Position for look-ahead.
		Eigen::Matrix<PosType, Dims, 1> vec_dir(pos_);

		for (INT axis = 0; axis < dims.size(); axis++)
		{
			vec_dir(axis) += 1;
			vec_grad(axis) = (*this)(vec_dir) - val_centre;
			vec_dir(axis) -= 1;
		}

		return vec_grad / this->dx();
	}

	/**
	 * Backward difference gradient, \f$ \nabla \phi \f$ .
	 *
	 * @param pos_ position in grid to query.
	 * @return vector with backward difference gradient.
	 * @return
	 */
	template <typename PosType>
	VecDT gradB (const felt::VecDT<PosType, Dims>& pos_) const
	{
		// Value at this point.
		const LeafType val_centre = (*this)(pos_);
		// Reference to GridBase dimensions.
		const VecDu& dims = this->dims();
		// Vector to store gradient calculation.
		VecDT vec_grad(dims.size());
		// Position for look-behind.
		Eigen::Matrix<PosType, Dims, 1> vec_dir(pos_);

		for (INT axis = 0; axis < dims.size(); axis++) {
			vec_dir(axis) -= 1;
			vec_grad(axis) = val_centre - (*this)(vec_dir);
			vec_dir(axis) += 1;
		}

		return vec_grad / this->dx();
	}

	/**
	 * Central difference gradient, \f$ \nabla \phi \f$ .
	 *
	 * @param pos_ position in grid to query.
	 * @return vector with central difference gradient.
	 */
	template <typename PosType>
	VecDT gradC (const felt::VecDT<PosType, Dims>& pos_) const
	{
		// Reference to GridBase dimensions.
		const VecDu& dims = this->dims();
		// Vector to store gradient calculation.
		VecDT vec_grad(dims.size());
		// Position for look-around.
		felt::VecDT<PosType, Dims> vec_dir(pos_);

		for (INT axis = 0; axis < dims.size(); axis++) {
			vec_dir(axis) -= 1;
			const LeafType back = (*this)(vec_dir);
			vec_dir(axis) += 2;
			const LeafType forward = (*this)(vec_dir);
			vec_dir(axis) -= 1;

			vec_grad(axis) = (forward - back) /  2;
		}

		return vec_grad / this->dx();
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
	VecDT grad (const felt::VecDT<PosType, Dims>& pos_) const
	{
		using VecDR = felt::VecDT<PosType, Dims>;
		// Reference to GridBase dimensions.
		const VecDu& dims = this->dims();
		// Vector to store gradient calculation.
		VecDT vec_grad(dims.size());
		// Position for look-around.
		VecDR pos_test(pos_);

		// Central value.
		LeafType centre = this->get(pos_);

		for (INT axis = 0; axis < dims.size(); axis++)
		{
			LeafType back = centre;
			LeafType forward = centre;
			UINT order = 0;
			// Check if backward value is within GridBase.
			pos_test(axis) -= 1;
			if (this->inside(pos_test))
			{
				back = this->get(pos_test);
				order++;
			}
			// Check if forward value is within GridBase.
			pos_test(axis) += 2;
			if (this->inside(pos_test))
			{
				forward = this->get(pos_test);
				order++;
			}
			pos_test(axis) -= 1;
			// Calculate central/forward/backward difference along this
			// axis.
			if (order != 0)
				vec_grad(axis) = (forward - back) / order;
			else
				vec_grad(axis) = 0;
		}

		return vec_grad / this->dx();
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
	VecDT gradE (const Eigen::Matrix<PosType, Dims, 1>& pos) const
	{
		typedef Eigen::Matrix<PosType, Dims, 1> VecDp;
		// Value at this point.
		const LeafType centre = (*this)(pos);
		// Reference to GridBase dimensions.
		const VecDu& dims = this->dims();
		// Vector to store gradient calculation.
		VecDT vec_grad(dims.size());
		// Position for look-around.
		VecDp pos_test(pos);

		for (INT axis = 0; axis < dims.size(); axis++)
		{
			pos_test(axis) -= 1;
			LeafType back = this->get(pos_test);
			pos_test(axis) += 2;
			LeafType forward = this->get(pos_test);
			pos_test(axis) -= 1;

			back = std::max((centre - back), 0.0f);
			forward = std::min(forward - centre, 0.0f);

			vec_grad(axis) = forward + back;
		}

		return vec_grad / this->dx();
	}

	/**
	 * Linear interpolation.
	 *
	 * @param pos_ position in grid to query.
	 * @return interpolated value at given position.
	 */
	LeafType interp (const VecDf& pos_) const
	{
		const VecDu& dims = this->dims();

		// Store all 2^d corners.
		std::vector< LeafType > val_corners(1 << dims.size());

		// Get all corners of containing cell.
		for (UINT i = 0; i < val_corners.size(); i++)
		{
			// 0 = 00 => (x,y)
			// 1 = 01 => (x+1,y)
			// 2 = 10 => (x,y+1)
			// 3 = 11 => (x+1,y+1)

			VecDi pos_corner(dims.size());
			for (INT dim = 0; dim < pos_corner.size(); dim++)
			{
				INT pos = (INT)std::floor(pos_(dim));
				const INT dir = (i >> dim) & 1;
				if (dir)
				{
					pos += dir;
					if (m_pos_min(dim) > pos || pos >= m_pos_max(dim))
						pos -= dir;
				}
				pos_corner(dim) = pos;
			}

			val_corners[i] = (*this)(pos_corner);
		}

		// Translate position vector into 'hypercube space',
		// so 0 <= v(x) <= 1.
		VecDf dir = pos_ - floorf(pos_);

		// Repeatedly reduce along axes,
		// i.e. hypercube -> cube -> square -> line -> point
		while (val_corners.size() > 1)
		{
			val_corners = this->interp(val_corners, dir);
		}

		return val_corners[0];
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
	LeafType curv (const felt::VecDT<PosType, Dims>& pos_) const
	{
		using VecDp = felt::VecDT<PosType, Dims>;

		const LeafType val_centre = (*this)(pos_);
		const VecDu& dims = this->dims();
		VecDp dir(pos_);

		// Forward directed principal normal.
		VecDT n_forward;

		for (INT axis = 0; axis < dims.size(); axis++)
		{
			dir(axis) += 1;

			const LeafType val_axis = (*this)(dir) - val_centre;
			LeafType val_neighs_sq = 0;

			// Loop other dimensions to get central difference across them.
			for (INT axis_neigh = 0; axis_neigh < dims.size(); axis_neigh++)
			{
				// Only getting differences across other axes.
				if (axis_neigh != axis)
				{
					// Central difference across this forward point.
					VecDp dir_neigh(dir);
					dir_neigh(axis_neigh) -= 1;
					const LeafType val_low = this->get(dir_neigh);
					dir_neigh(axis_neigh) += 2;
					const LeafType val_high = this->get(dir_neigh);

					const LeafType val_neigh = (val_high - val_low) / 2;
					val_neighs_sq += val_neigh*val_neigh;
				}
			}

			n_forward(axis) =		val_axis /
						sqrt(val_axis*val_axis + val_neighs_sq);

			dir(axis) -= 1;
		}


		// Backward directed principal normal.
		VecDT n_backward;

		for (INT axis = 0; axis < dims.size(); axis++)
		{
			dir(axis) -= 1;

			const LeafType val_axis = val_centre - (*this)(dir);
			LeafType val_neighs_sq = 0;

			// Loop other dimensions to get central difference across them.
			for (
				INT axis_neigh = 0; axis_neigh < dims.size(); axis_neigh++
			) {
				// Only getting differences across other axes.
				if (axis_neigh != axis)
				{
					// Central difference across this backward point.
					VecDp dir_neigh(dir);
					dir_neigh(axis_neigh) -= 1;
					const LeafType val_low = this->get(dir_neigh);
					dir_neigh(axis_neigh) += 2;
					const LeafType val_high = this->get(dir_neigh);

					const LeafType val_neigh = (val_high - val_low) / 2;
					val_neighs_sq += val_neigh*val_neigh;
				}
			}

			n_backward(axis) =		val_axis /
						sqrt(val_axis*val_axis + val_neighs_sq);

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
	LeafType divergence (const felt::VecDT<PosType, Dims>& pos_) const
	{
		const VecDT vec_grad_f = this->gradF(pos_);
		const VecDT vec_grad_b = this->gradB(pos_);
		const VecDT vec_grad_diff = vec_grad_b - vec_grad_f;

		// Component-wise sum.
		const LeafType val = vec_grad_diff.sum();

		return val / (this->dx()*this->dx());
	}

	/**
	 * Get the value stored in the grid, circumventing subclass's mutation.
	 *
	 * @param pos_ position in grid to query.
	 * @return internally stored value at given grid position
	 */
	LeafType& get_internal (const VecDi& pos_) {
		#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
		assert_pos_bounds(pos_, "get_internal: ");
		#endif
		const UINT& idx = this->index(pos_);
		return this->data()(idx);
	}

	/**
	 * Get the value stored in the grid, circumventing subclass's mutation
	 * (const version).
	 *
	 * @param pos_ position in grid to query.
	 * @return internally stored value at given grid position
	 */
	const LeafType& get_internal (const VecDi& pos) const
	{
		#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
		assert_pos_bounds(pos, "get_internal: ");
		#endif
		const UINT& idx = this->index(pos);
		return this->data()(idx);
	}

	/**
	 * Check if given position is within the grid and raise a domain_error if not.
	 *
	 * @param pos_ position in grid to query.
	 * @param title_ message to include in generated exception.
	 */
	void assert_pos_bounds (const VecDi& pos_, std::string title_) const
	{
		const DerivedType* self = static_cast<const DerivedType*>(this);

		if (!self->inside(pos_))
		{
			const VecDi& pos_min = self->offset();
			const VecDi& pos_max = (
				self->dims().template cast<INT>() + pos_min - VecDi::Constant(1)
			);
			std::stringstream err;
			err << title_ << format(pos_.transpose())
				<< " is outside grid "
				<< format(pos_min) << "-" << format(pos_max)
				<< std::endl;
			std::string err_str = err.str();
			throw std::domain_error(err_str);
		}
	}

#ifndef _TESTING
protected:
#else
public:
#endif

	/**
	 * Interpolate down one dimension.
	 *
	 * The values of val_corners_in_ are interpolated to one dimension smaller than they are
	 * currently (cube->square, square->line, line->point).
	 *
	 * @param val_corners_in_ list of corner values.
	 * @param pos_ real-valued position to interpolate to.
	 * @return
	 */
	std::vector<LeafType> interp (
		const std::vector<LeafType>& val_corners_in_, const VecDf& pos_
	) const
	{
		const size_t num_corners = val_corners_in_.size();

		// Number of values returned.
		// This is a power of 2 less than input dimensions
		// (cube becomes square, square becomes line, line becomes point).
		const size_t num_out = num_corners >> 1;

		// The axis along which to interpolate.
		// This is computed from the dimensions of the original input and
		// the dimensions of the intended output.
		const size_t axis_idx = pos_.size() - log2(num_corners);

		// The weighting to be used in interpolating each pair of points.
		// This is the position along the axis of interpolation.
		const FLOAT axis_pos = pos_(axis_idx);

		std::vector<LeafType> val_corners_out(num_out);

		for (size_t i = 0; i < num_out; i++)
		{
			const LeafType low = val_corners_in_[(i << 1)];
			const LeafType high = val_corners_in_[(i << 1) + 1];
			const LeafType val = axis_pos*high + (1.0f-axis_pos)*low;
			val_corners_out[i] = val;
		}

		return val_corners_out;
	}

};


/**
 * A standard D-dimensional grid for storing values of type T.
 *
 * @tparam T the type of data to store in the grid.
 * @tparam D the dimensions of the grid.
 */
template <typename T, UINT D>
class Grid : public GridBase<Grid<T, D> >
{
public:
	using ThisType = Grid<T, D>;
	using Base = GridBase<ThisType>;
	using typename Base::VecDu;
	using typename Base::VecDi;
	using typename Base::VecDf;
	using Base::GridBase;

	~Grid ()
	{}

	/**
	 * Override GridBase::get to simply return the value stored in the grid.
	 *
	 * @param pos_ position in grid to query.
	 * @return internally stored value at given grid position
	 */
	T& get (const VecDi& pos_)
	{
		return this->get_internal(pos_);
	}

	/**
	 * Override GridBase::get to simply return the value stored in the grid
	 * (const version).
	 *
	 * @param pos_ position in grid to query.
	 * @return internally stored value at given grid position
	 */
	const T& get (const VecDi& pos_) const
	{
		return this->get_internal(pos_);
	}
};

/**
 * Traits for GridBase to understand Grid.
 *
 * @tparam T the type of data to store in the grid.
 * @tparam D the number of dimensions of the grid.
 */
template <typename T, UINT D>
struct GridBaseTraits<Grid<T, D> >
{
	/// The class inheriting from the base.
	using ThisType = Grid<T, D>;
	/// Dimensions of the grid, from template parameter.
	static const UINT Dims = D;
	/// The data type to store at leaf grid nodes, from template parameter.
	using LeafType = T;
	/// The data type to return when leaf grid nodes are queried, from template parameter.
	using RetType = T;
};

/** @} */ // End group Grid.

}// End namespace felt.
#endif


