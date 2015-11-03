#ifndef MAPPEDGRID_HPP_
#define MAPPEDGRID_HPP_

#include <array>
#include <sstream>
#include <mutex>
#include "Grid.hpp"

namespace felt
{

/** @defgroup LookupGrids
 *
 *  Bi-directional lookup between tracking list and grid locations - grids index lists and lists
 *  index grids.
 *  @{
 */

/**
 * Base traits class for classes CRTP derived from LookupGridBase.
 */
template <class Derived> struct LookupGridBaseTraits;


/**
 * Base class for a lookup grid.
 * 
 * Array elements store grid positions and grid nodes store array indices.  Each grid node holds an
 * n-tuple of indices, and there are n arrays associated with this lookup grid.  An alternative 
 * would be to use multiple grids each tracking a single list, but a multi-list structure allows
 * related lists to be spatially tracked in contiguous memory.
 *
 * @tparam Derived the CRTP derived class
 */
template <class Derived>
class LookupGridBase : public GridBase<LookupGridBase<Derived> >
{
public:
	using ThisType = LookupGridBase<Derived>;
	/// Type of data stored in grid nodes. For lookup grids this is an N-tuple of array indices.
	using LeafType = typename LookupGridBaseTraits<Derived>::LeafType;
	/// Type of data to return when grid nodes are queried.  Either array index or indices.
	using RetType = typename LookupGridBaseTraits<Derived>::RetType;
	/// GridBase base class.
	using Base = GridBase<ThisType>;
	using typename Base::VecDu;
	using typename Base::VecDi;
	using typename Base::PosArray;

	/// A tuple of array indices indicating all NULL indices (nothing pointed to).
	static const LeafType	NULL_IDX_TUPLE;
	/// An array index indicating a NULL index (nothing pointed to).
	static const UINT		NULL_IDX;
	/// Number of lists tracked in this grid.
	static const UINT 		NUM_LISTS = LookupGridBaseTraits<Derived>::NumLists;

protected:
	/// N-tuple of lists of grid positions - the tracking lists.
	std::array<PosArray, NUM_LISTS>	m_a_pos;
	/// Mutex for use by other classes where multiple threads hold a reference to this grid.
	std::mutex						m_mutex;
public:

	/**
	 * Default (empty) destructor.
	 */
	~LookupGridBase() {}

	/**
	 * Explicitly defined default constructor.
	 */
	LookupGridBase () = default;
	
	/**
	 * Construct a lookup grid of given size and spatial offset.
	 *
	 * @param size_ size of grid.
	 * @param offset_ spatial offset of grid.
	 */
	LookupGridBase (const VecDu& size_, const VecDi& offset_ = VecDi::Zero()) : Base()
	{
		this->init(size_, offset_);
	}

	/**
	 * Copy constructor.
	 *
	 * Required since destructor defined and mutex is non-copyable.
	 *
	 * @param other_ grid to copy.
	 */
	LookupGridBase(const ThisType& other_)
	{
		m_a_pos = other_.m_a_pos;
		this->m_data = other_.m_data;
		this->m_offset = other_.m_offset;
		this->m_dims = other_.m_dims;
	}

	/**
	 * Move constructor,
	 *
	 * Required since destructor defined and mutex is non-movable.
	 *
	 * @param other_ grid to move.
	 */
	LookupGridBase(ThisType&& other_)
	{
		m_a_pos = std::move(other_.m_a_pos);
		this->m_data = std::move(other_.m_data);
		this->m_offset = std::move(other_.m_offset);
		this->m_dims = std::move(other_.m_dims);
	}

	/**
	 * Copy assigment operator.
	 *
	 * Required since destructor defined and mutex is non-movable.
	 *
	 * @param other_ grid to copy.
	 */
	ThisType& operator=(const ThisType& other_)
	{
		m_a_pos = other_.m_a_pos;
		this->m_data = other_.m_data;
		this->m_offset = other_.m_offset;
		this->m_dims = other_.m_dims;
		return *this;
	}

	/**
	 * Move assigment operator.
	 *
	 * Required since destructor defined and mutex is non-movable.
	 *
	 * @param other_ grid to move.
	 */
	ThisType& operator=(ThisType&& other_)
	{
		m_a_pos = std::move(other_.m_a_pos);
		this->m_data = std::move(other_.m_data);
		this->m_offset = std::move(other_.m_offset);
		this->m_dims = std::move(other_.m_dims);
		return *this;
	}

	/**
	 * Return tuple of indices stored at given position in grid.
	 *
	 * Overloads get method of GridBase to simply return value stored.
	 *
	 * @param pos_ position in grid to query.
	 * @return tuple of indices at this grid position.
	 */
	RetType& get (const VecDi& pos_)
	{
		return this->get_internal(pos_);
	}

	/**
	 * Return tuple of indices stored at given position in grid.
	 *
	 * Overloads get method of GridBase to simply return value stored.
	 *
	 * @param pos_ position in grid to query.
	 * @return tuple of indices at this grid position.
	 */
	const RetType& get (const VecDi& pos_) const
	{
		return this->get_internal(pos_);
	}

	/**
	 * Get mutex associated with this lookup grid - only used externally.
	 *
	 * Adding and removing elements from the tracking list/grid is not thread safe, so if multiple
	 * threads use this lookup grid, this mutex must be used.
	 *
	 * @return mutex for use by external objects.
	 */
	std::mutex& mutex ()
	{
		return m_mutex;
	}

	/**
	 * Reshape grid and fill with NULL indices.
	 *
	 * @param size_ new size of the grid.
	 */
	void dims (const VecDu& size_)
	{
		Base::dims(size_);
		this->fill(NULL_IDX_TUPLE);
	}

	using Base::dims;

	/**
	 * Get tracking list by id.
	 *
	 * @param arr_idx_ id of tracking list to get.
	 * @return tracking list at given index.
	 */
	inline PosArray& list (const UINT& arr_idx_ = 0)
	{
		return m_a_pos[arr_idx_];
	}

	/**
	 * Get tracking list by id.
	 *
	 * @param arr_idx_ id of tracking list to get.
	 * @return tracking list at given index.
	 */
	inline const PosArray& list (const UINT& arr_idx_ = 0) const
	{
		return m_a_pos[arr_idx_];
	}

	/**
	 * Return true if position currently tracked for given list id, false
	 * otherwise.
	 *
	 * @param pos_ position in grid to query.
	 * @param arr_idx_ id of tracking list to test.
	 * @return true if grid position tracked, false otherwise.
	 */
	const bool is_active (const VecDi& pos, const UINT& arr_idx_ = 0) const
	{
		return this->get(pos)[arr_idx_] != NULL_IDX;
	}

	/**
	 * Add position to tracking list and store index in tracking list in
	 * grid.
	 *
	 * Does nothing if grid node already set with an index in tracking list.
	 *
	 * @param pos_ position in grid to track.
	 * @param arr_idx_ id of tracking list to add to.
	 * @return true if grid node set and position added to list,
	 * false if grid node was already set so position already in a list.
	 */
	bool add (const VecDi& pos, const UINT& arr_idx_ = 0)
	{
		return add(pos, arr_idx_, arr_idx_);
	}

	/**
	 * Remove an element from a tracking list by index and set it's
	 * corresponding grid node to NULL index.
	 *
	 * @param idx_ index in tracking list.
	 * @param arr_idx_ tracking list id.
	 */
	void remove (const UINT& idx, const UINT& arr_idx_ = 0)
	{
		const VecDi& pos = this->list(arr_idx_)[idx];
		remove(idx, pos, arr_idx_, arr_idx_);
	}

	/**
	 * Look up tracking list index in grid, remove from list and set grid
	 * node to NULL index.
	 *
	 * @param pos_ position in lookup grid.
	 * @param arr_idx_ tracking list id.
	 */
	void remove (const VecDi& pos, const UINT& arr_idx_ = 0)
	{
		const UINT& idx = this->get_internal(pos)[arr_idx_];
		if (idx == NULL_IDX)
			return;
		remove(idx, pos, arr_idx_, arr_idx_);
	}

	/**
	 * Clear tracking list and reset every grid point to NULL index.
	 *
	 * @param arr_idx_ tracking list id to clear.
	 */
	void reset (const UINT& arr_idx_ = 0)
	{
		reset(arr_idx_, arr_idx_);
	}

protected:

	/**
	 * Add pos to tracking list and set pos in grid to index in tracking
	 * list.
	 *
	 * If a grid node has a non-NULL index then does nothing.
	 *
	 * The arr_idx and lookup_idx can be different values - used in
	 * subclass overrides.  See LookupGrid.
	 *
	 * @param pos_ position in lookup grid.
	 * @param arr_idx_ tracking list id.
	 * @param lookup_idx_ lookup grid id.
	 * @return true if grid node was set and position added to list,
	 * false if grid node was already set so position already in a list.
	 */
	bool add(const VecDi& pos_, const UINT& arr_idx_, const UINT& lookup_idx_)
	{
		#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
		this->assert_pos_bounds(pos_, "add: ");
		#endif

		UINT& idx = this->get_internal(pos_)[lookup_idx_];
		// Do not allow duplicates.
		if (idx != NULL_IDX)
			return false;
		// idx is by reference, so this sets the value in grid.
		idx = this->list(arr_idx_).size();
		this->list(arr_idx_).push_back(pos_);
		return true;
	}

	/**
	 * Set all lookup grid nodes to NULL index and clear all lists.
	 *
	 * @param arr_idx_ tracking list id.
	 */
	void reset_all ()
	{
		for (UINT idx = 0; idx < m_a_pos.size(); idx++)
			reset(idx, idx);
	}

	/**
	 * For given tracking list, set all lookup grid nodes to NULL index and
	 * clear the list.
	 *
	 * @param arr_idx_ tracking list id.
	 * @param lookup_idx_ lookup grid node tuple index.
	 */
	void reset(const UINT& arr_idx_, const UINT& lookup_idx_)
	{
		for (VecDi pos : m_a_pos[arr_idx_])
			this->get_internal(pos)[lookup_idx_] = NULL_IDX;
		this->list(arr_idx_).clear();
	}

	/**
	 * Remove pos at index idx in the array and set lookup at pos to NULL index.
	 *
	 * NOTE: idx passed by value since it changes indirectly via grid lookup.
	 *
	 * @param idx_ index of element in tracking list to remove.
	 * @param pos_ position in grid matching index in tracking list to remove.
	 * @param arr_idx_ tracking list id to remove element from.
	 * @param lookup_idx_ index in tuple stored at grid node that references element to remove.
	 */
	void remove(
		const UINT idx_, const VecDi& pos_, const UINT& arr_idx_, const UINT& lookup_idx_
	) {
		#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
		this->assert_pos_bounds(pos_, "remove: ");
		#endif

		// Set index lookup to null value.
		this->get_internal(pos_)[lookup_idx_] = NULL_IDX;

		// If this is not the last remaining position in the array, then
		// we must move the last position to this position and update the
		// lookup grid.
		const UINT& size = this->list(arr_idx_).size();
		if (idx_ < size - 1)
		{
			// Duplicate last element into this index.
			const VecDi& pos_last = this->list(arr_idx_)[size - 1];
			this->list(arr_idx_)[idx_] = pos_last;
			// Set the lookup grid to reference the new index in the array.
			this->get_internal(pos_last)[lookup_idx_] = idx_;
		}
		// Remove the last element in the array (which is at this point
		// either the last remaining element or a duplicate).
		this->list(arr_idx_).pop_back();
	}
};


template <class Derived> const typename LookupGridBase<Derived>::LeafType
LookupGridBase<Derived>::NULL_IDX_TUPLE = LookupGridBase<Derived>::LeafType::Constant(
	std::numeric_limits<UINT>::max()
);


template <class Derived> const UINT
LookupGridBase<Derived>::NULL_IDX = std::numeric_limits<UINT>::max();


/**
 * Traits for GridBase to understand LookupGridBase.
 *
 * Just forward the traits defined for LookupGridBase subclasses.
 */
template <class Derived>
struct GridBaseTraits<LookupGridBase<Derived> >
{
	/// Dimension of grid.
	static const UINT Dims = LookupGridBaseTraits<Derived>::Dims;
	using LeafType = typename LookupGridBaseTraits<Derived>::LeafType;
	using RetType = typename LookupGridBaseTraits<Derived>::RetType;
	using ThisType = typename LookupGridBaseTraits<Derived>::ThisType;
};


/**
 * Standard lookup grid.
 *
 * Holds a set of tracking lists storing grid positions, and a corresponding grid storing tuples of
 * list indices.
 *
 * @tparam D the dimension of the grid.
 * @tparam N the number of tracking lists to use.
 */
template <UINT D, UINT N=1>
class LookupGrid : public LookupGridBase<LookupGrid<D, N> >
{
protected:
	using ThisType = LookupGrid<D, N>;
	using Base = LookupGridBase<ThisType>;
	using typename Base::VecDu;
	using typename Base::VecDi;

public:
	using Base::LookupGridBase;
	using Base::reset;
};


/**
 * Traits of standard LookupGrid for CRTP inheritance from LookupGridBase.
 *
 * @tparam D the dimension of the grid.
 * @tparam N the number of tracking lists to use.
 */
template <UINT D, UINT N>
struct LookupGridBaseTraits<LookupGrid<D, N> >
{
	using ThisType = LookupGrid<D, N>;
	/// Leaf grid nodes store N-tuple of list indices.
	using LeafType = VecDu<N>;
	/// Grid returns the same type as is stored, N-tuple of list indices.
	using RetType = LeafType;
	/// Dimension of the grid taken from template parameter.
	static const UINT Dims = D;
	/// Number of tracking lists taken from template parameter.
	static const UINT NumLists = N;
};


/**
 * Base traits class for classes CRTP derived from SharedLookupGridBase.
 */
template <class Derived> struct SharedLookupGridBaseTraits {};


/**
 * Similar to LookupGrid but grid nodes store only a single list index.
 *
 * Useful in cases where grid nodes cannot be in more than one list.
 *
 * @tparam D the dimension of the grid.
 * @tparam N the number of tracking lists to use.
 */
template <class Derived>
class SharedLookupGridBase : public LookupGridBase<SharedLookupGridBase<Derived> >
{
public:
	using ThisType = SharedLookupGridBase<Derived>;
	using Base = LookupGridBase<ThisType>;
	using typename Base::RetType;
	using typename Base::VecDu;
	using typename Base::VecDi;
protected:
	using Base::m_a_pos;
public:
	using Base::LookupGridBase;

	/**
	 * Return index in associated list of grid node.
	 *
	 * Overloads get method of GridBase to return the first (and only)
	 * element in the tuple of indices.
	 *
	 * @param pos_ position in grid to query.
	 * @return list index stored at given position
	 */
	RetType& get (const VecDi& pos_)
	{
		return this->get_internal(pos_)[0];
	}

	/**
	 * Return index in associated list of grid node.
	 *
	 * Overloads get method of GridBase to return the first (and only)
	 * element in the tuple of indices.
	 *
	 * @param pos_ position in grid to query.
	 * @return list index stored at given position
	 */
	const RetType& get (const VecDi& pos_) const
	{
		return this->get_internal(pos_)[0];
	}

	/**
	 * Return true if position currently tracked for given list id, false
	 * otherwise.
	 *
	 * @param pos_ position in grid to query.
	 * @return true if position is tracked, false otherwise.
	 */
	const bool is_active (const VecDi& pos_) const
	{
		return this->get(pos_) != Base::NULL_IDX;
	}

	/**
	 * Add position to tracking list and store index in tracking list in grid.
	 *
	 * Overloads LookupGridBase to place lookup index in first (and only) element of tuple of
	 * indices at the grid position.
	 *
	 * @param pos_ position in grid to track.
	 * @param arr_idx_ tracking list id to add position to.
	 * @return true if grid node set and position added to list, false if grid node was already set
	 * so position already in a list.
	 */
	bool add (const VecDi& pos_, const UINT& arr_idx_ = 0)
	{
		return Base::add(pos_, arr_idx_, 0u);
	}

	/**
	 * For given tracking list, set all lookup grid nodes to NULL index and
	 * clear the list.
	 *
	 * Overloads LookupGridBase to reset first (and only) element of tuple
	 * of indices for each grid position referenced in given tracking list.
	 *
	 * @param arr_idx_ tracking list id.
	 */
	void reset (const UINT& arr_idx_ = 0)
	{
		Base::reset(arr_idx_, 0u);
	}

	/**
	 * Set all lookup grid nodes to NULL index and clear all lists.
	 */
	void reset_all ()
	{
		for (UINT idx = 0; idx < m_a_pos.size(); idx++)
			Base::reset(idx, 0u);
	}

	/**
	 * Remove an element from a tracking list by index and set it's
	 * corresponding grid node to NULL index.
	 *
	 * Overloads LookupGridBase to set the first (and only) element of
	 * tuple of indices in the grid to NULL index.
	 *
	 * @param idx_ index in tracking list.
	 * @param arr_idx_ tracking list id.
	 */
	void remove (const UINT& idx_, const UINT& arr_idx_ = 0)
	{
		const VecDi& pos = this->list(arr_idx_)[idx_];
		Base::remove(idx_, pos, arr_idx_, 0);
	}

	/**
	 * Look up tracking list index in grid, remove from list and set grid
	 * node to NULL index.
	 *
	 * Overloads LookupGridBase to set the first (and only) element of
	 * tuple of indices in the grid to NULL index.
	 *
	 * @param pos_ position in lookup grid.
	 * @param arr_idx_ tracking list id.
	 */
	void remove (const VecDi& pos_, const UINT& arr_idx_ = 0)
	{
		const UINT& idx = this->get_internal(pos_)[0];
		Base::remove(idx, pos_, arr_idx_, 0);
	}
};


/**
 * Traits for LookupGridBase to understand SharedLookupGridBase.
 *
 * Just forward the traits defined for SharedLookupGridBase subclasses.
 */
template <class Derived>
struct LookupGridBaseTraits<SharedLookupGridBase<Derived> >
{
	using ThisType = typename SharedLookupGridBaseTraits<Derived>::ThisType;
	using LeafType = typename SharedLookupGridBaseTraits<Derived>::LeafType;
	using RetType = typename SharedLookupGridBaseTraits<Derived>::RetType;
	static const UINT Dims = SharedLookupGridBaseTraits<Derived>::Dims;
	static const UINT NumLists = SharedLookupGridBaseTraits<Derived>::NumLists;
};


/**
 * Concrete definition of standard SharedLookupGrid from SharedLookupGridBase.
 *
 * A simple stub to expose SharedLookupGridBase, which is also CRTP derived from elsewhere
 * - @see SharedTrackedPartitionedGrid.
 *
 * @tparam D the dimension of the grid.
 * @tparam N the number of tracking lists to use.
 */
template <UINT D, UINT N=1>
class SharedLookupGrid : public SharedLookupGridBase<SharedLookupGrid<D, N> >
{
public:
	using ThisType = SharedLookupGrid<D, N>;
	using Base = SharedLookupGridBase<ThisType>;
	using Base::SharedLookupGridBase;
};


/**
 * Traits of SharedLookupGrid for CRTP inheritance from SharedLookupGridBase.
 *
 * @tparam D dimension of the grid.
 * @tparam N number of tracking lists to use.
 */
template <UINT D, UINT N>
struct SharedLookupGridBaseTraits<SharedLookupGrid<D, N> >
{
	using ThisType = SharedLookupGrid<D, N>;
	/// Shared lookup grids only store a single index at each grid position.
	using LeafType = VecDu<1>;
	/// Grid queries return the index stored at that position.
	using RetType = UINT;
	/// Dimension of the grid taken from template parameter.
	static const UINT Dims = D;
	/// Number of tracking lists taken from template parameter.
	static const UINT NumLists = N;

};

/** @} */ // End group LookupGrid.

/** @defgroup TrackedGrids
 *
 *  Grids storing arbitrary data, using paired lookup grids to track active points.
 *
 *  @{
 */


/**
 * Base traits class for classes CRTP derived from TrackedGridBase.
 */
template <class Derived> struct TrackedGridBaseTraits {};



/**
 * Base class for a tracking grid
 *
 * Grid nodes store arbitrary values and active nodes are tracked by a LookupGrid.
 *
 * @See TrackedGrid and SharedTrackedGrid.
 */
template <class Derived>
class TrackedGridBase : public GridBase<TrackedGridBase<Derived> >
{
public:
	/// GridBase base class.
	using Base = GridBase<TrackedGridBase<Derived> >;
	/// Lookup grid type to use for tracking active grid positions.
	using Lookup = typename TrackedGridBaseTraits<Derived>::LookupType;
	/// Type of data to store in the main grid.
	using LeafType = typename TrackedGridBaseTraits<Derived>::LeafType;
	/// Tracking list of grid positions.
	using PosArray = typename Lookup::PosArray;
	using typename Base::VecDu;
	using typename Base::VecDi;
protected:
	/// Mutex for use by other classes where multiple threads hold a reference to this grid.
	std::mutex	m_mutex;
	/// Internal lookup grid to track active grid positions.
	Lookup		m_grid_lookup;
public:
	using Base::offset;
	using Base::dims;

	/**
	 * Explicitly defined default constructor.
	 */
	TrackedGridBase() = default;

	/**
	 * Construct a grid with given size and spatial offset.
	 *
	 * @param size_ size of grid.
	 * @param offset_ spatial offset of grid.
	 */
	TrackedGridBase(const VecDu& size_, const VecDi& offset_) : Base(), m_grid_lookup()
	{
		this->init(size_, offset_);
	}

	/**
	 * Get mutex associated with this lookup grid - only used externally.
	 *
	 * Adding and removing elements from the tracking list/grid is not thread safe, so if multiple
	 * threads use this lookup grid, this mutex must be used.
	 *
	 * @return mutex for use by external objects.
	 */
	std::mutex& mutex ()
	{
		return m_mutex;
	}

	/**
	 * Reshape both grid of values and lookup grid.
	 *
	 * Lookup grid nodes will be reset to NULL.
	 *
	 * @param size_ new size of the grid.
	 */
	void dims (const VecDu& size_)
	{
		Base::dims(size_);
		m_grid_lookup.dims(size_);
	}

	/**
	 * Set offset of grid of values and lookup grid.
	 *
	 * @param offset_ spatial offset of grid.
	 */
	void offset (const VecDi& offset_)
	{
		Base::offset(offset_);
		m_grid_lookup.offset(offset_);
	}

	/**
	 * Override GridBase::get to simply return the value stored in the grid.
	 *
	 * @param pos_ position in grid to query.
	 * @return value stored in main grid at given position.
	 */
	LeafType& get (const VecDi& pos_)
	{
		return this->get_internal(pos_);
	}

	/**
	 * Override GridBase::get to simply return the value stored in the grid
	 * (const version).
	 *
	 * @param pos_ position in grid to query.
	 * @return value stored in main grid at given position.
	 */
	const LeafType& get (const VecDi& pos_) const
	{
		return this->get_internal(pos_);
	}

	/**
	 * Get lookup grid.
	 *
	 * @return the internal lookup grid tracking active grid positions.
	 */
	Lookup& lookup()
	{
		return m_grid_lookup;
	}

	/**
	 * Get lookup grid.
	 *
	 * @return the internal lookup grid tracking active grid positions.
	 */
	const Lookup& lookup() const
	{
		return m_grid_lookup;
	}

	/**
	 * Get list of active grid points from lookup grid.
	 *
	 * @param arr_idx_ tracking list id.
	 * @return the tracking list of active grid positions from the internal lookup grid.
	 */
	PosArray& list(const UINT& arr_idx_ = 0)
	{
		return m_grid_lookup.list(arr_idx_);
	}

	/**
	 * Get list of active grid points from lookup grid.
	 *
	 * @param arr_idx_ tracking list id.
	 * @return the tracking list of active grid positions from the internal lookup grid.
	 */
	const PosArray& list(const UINT& arr_idx_ = 0) const
	{
		return m_grid_lookup.list(arr_idx_);
	}

	/**
	 * Set value in grid at given position and add position to lookup grid.
	 *
	 * Will set value regardless whether lookup grid already set for given
	 * position + tracking list.
	 *
	 * @param pos_ position in grid.
	 * @param val_ value to set.
	 * @param arr_idx_ tracking list id.
	 * @return true if grid node set in lookup grid and position added to
	 * tracking list, false if grid node was already set so position already
	 * in a list.
	 */
	bool add(const VecDi& pos_, const LeafType& val_, const UINT& arr_idx_ = 0)
	{
		this->get(pos_) = val_;
		return add(pos_, arr_idx_);
	}

	/**
	 * Add a position to the lookup grid.
	 *
	 * @param pos_ position in the grid to add.
	 * @param arr_idx_ tracking list id.
	 * @return true if grid node set in lookup grid and position added to
	 * tracking list, false if grid node was already set so position already
	 * in a list.
	 */
	bool add(const VecDi& pos_, const UINT& arr_idx_ = 0)
	{
		return m_grid_lookup.add(pos_, arr_idx_);
	}

	/**
	 * Set every active grid node (those referenced by lookup grid)
	 * to given value and reset the lookup grid.
	 *
	 * Lookup grid will then be full of NULL indices and it's tracking
	 * list(s) will be empty.
	 *
	 * @param val_ value to set in main grid.
	 * @param arr_idx_ tracking list id to cycle over and clear.
	 */
	void reset(const LeafType& val_, const UINT& arr_idx_ = 0)
	{
		for (VecDi pos : m_grid_lookup.list(arr_idx_))
			this->get(pos) = val_;
		reset(arr_idx_);
	}

	/**
	 * Reset a tracking list on the lookup grid.
	 *
	 * @param arr_idx_ tracking list to clear.
	 */
	void reset(const UINT& arr_idx_ = 0)
	{
		m_grid_lookup.reset(arr_idx_);
	}

	/**
	 * Remove an element from a tracking list by index and set it's
	 * corresponding grid node in the tracking grid to NULL index.
	 *
	 * @param idx_ index in tracking list.
	 * @param arr_idx_ tracking list id.
	 */
	void remove(const UINT& idx_, const UINT& arr_idx_ = 0)
	{
		m_grid_lookup.remove(idx_, arr_idx_);
	}

	/**
	 * Look up tracking list index in grid, remove from list and set lookup
	 * grid node to NULL index.
	 *
	 * @param pos_ position in lookup grid.
	 * @param arr_idx_ tracking list id.
	 */
	void remove (const VecDi& pos_, const UINT& arr_idx_ = 0)
	{
		m_grid_lookup.remove(pos_, arr_idx_);
	}

	/**
	 * Return true if position currently tracked for given list id, false
	 * otherwise.
	 *
	 * @param pos_ position in grid to query.
	 * @param arr_idx_ tracking list id to query.
	 * @return true if position currently tracked, false otherwise.
	 */
	const bool is_active (const VecDi& pos_, const UINT& arr_idx_ = 0) const
	{
		return m_grid_lookup.is_active(pos_, arr_idx_);
	}
};


/**
 * Traits for GridBase to understand TrackedGridBase.
 *
 * Just forward the traits defined for TrackedGridBase subclasses.
 */
template <class Derived>
struct GridBaseTraits<TrackedGridBase<Derived> >
{
	using ThisType = typename TrackedGridBaseTraits<Derived>::ThisType;
	using LeafType = typename TrackedGridBaseTraits<Derived>::LeafType;
	using RetType = LeafType;
	static const UINT Dims = TrackedGridBaseTraits<Derived>::Dims;
};


/**
 * Standard tracked grid.
 *
 * A grid of arbitrary data, with active positions tracked by an internal LookupGrid.
 *
 * Supports multiple tracking lists that can overlap (the same grid position can be found in
 * multiple lists).
 *
 * This is just a stub exposing TrackedGridBase, relying on the associated traits to differentiate
 * behaviour.
 *
 * @tparam T the type of data to store in the grid.
 * @tparam D the dimension of the grid.
 * @tparam N the number of tracking lists to use.
 */
template <typename T, UINT D, UINT N=1>
class TrackedGrid : public TrackedGridBase<TrackedGrid<T, D, N> >
{
public:
	using TrackedGridBase<TrackedGrid<T, D, N> >::TrackedGridBase;
};


/**
 * Traits of TrackedGrid for CRTP inheritance from TrackedGridBase.
 *
 * @tparam T the type of data to store in the grid.
 * @tparam D the dimension of the grid.
 * @tparam N the number of tracking lists to use.
 */
template <typename T, UINT D, UINT N>
struct TrackedGridBaseTraits<TrackedGrid<T, D, N> >
{
	using ThisType = TrackedGrid<T, D, N>;
	/// Type of data to store in main grid taken from template parameter.
	using LeafType = T;
	/// Dimension of the grid taken from template parameter.
	static const UINT Dims = D;
	/// Type of lookup grid to use.  This is what differentiates this from SharedTrackedGrid.
	using LookupType = LookupGrid<D, N>;
};


/**
 * A tracked grid that assumes non-overlapping tracking lists.
 *
 * A grid of arbitrary data, with active positions tracked by an internal SharedLookupGrid.
 *
 * Each node of the associated lookup grid stores only a single list index. A significant memory
 * saving when a grid node can only be in one of the tracking lists.
 *
 * This is just a stub exposing TrackedGridBase, relying on the associated traits to differentiate
 * behaviour.
 *
 * @tparam T the type of data to store in the grid.
 * @tparam D the dimension of the grid.
 * @tparam N the number of tracking lists to use.
 */
template <typename T, UINT D, UINT N=1>
class SharedTrackedGrid : public TrackedGridBase<SharedTrackedGrid<T, D, N> >
{
public:
	using TrackedGridBase<SharedTrackedGrid<T, D, N> >::TrackedGridBase;
};


/**
 * Traits of SharedTrackedGrid for CRTP inheritance from TrackedGridBase.
 *
 * @tparam T the type of data to store in the grid.
 * @tparam D the dimension of the grid.
 * @tparam N the number of tracking lists to use.
 */
template <typename T, UINT D, UINT N>
struct TrackedGridBaseTraits<SharedTrackedGrid<T, D, N> >
{
	using ThisType = SharedTrackedGrid<T, D, N>;
	/// Type of data to store in main grid taken from template parameter.
	using LeafType = T;
	/// Dimension of the grid taken from template parameter.
	static const UINT Dims = D;
	/// Type of lookup grid to use.  This is what differentiates this from TrackedGrid.
	using LookupType = SharedLookupGrid<D, N>;
};

/** @} */ // End group LookupGrid.

}
#endif /* MAPPEDGRID_HPP_ */
