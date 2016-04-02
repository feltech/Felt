/*
 * TrackedGrid.hpp
 *
 *  Created on: 10 Jan 2016
 *      Author: dave
 */

#ifndef INCLUDE_FELT_TRACKEDGRID_HPP_
#define INCLUDE_FELT_TRACKEDGRID_HPP_

#include <array>
#include <sstream>
#include <mutex>
#include "SharedLookupGrid.hpp"

namespace felt
{

/**
 * Base class for a tracking grid
 *
 * Grid nodes store arbitrary values and active nodes are tracked by a LookupGrid.
 *
 * @See TrackedGrid and SharedTrackedGrid.
 */
template <class Derived, bool IsLazy=false>
class TrackedGridBase : public GridBase<TrackedGridBase<Derived>, IsLazy>
{
public:
	/// GridBase base class.
	using Base = GridBase<TrackedGridBase<Derived>, IsLazy>;
	/// Lookup grid type to use for tracking active grid positions.
	using Lookup = typename GridTraits<Derived>::LookupType;
	/// Type of data to store in the main grid.
	using LeafType = typename GridTraits<Derived>::LeafType;
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
	using Base::size;
	using Base::GridBase;

	/**
	 * Explicitly defined default constructor.
	 */
	TrackedGridBase() : Base() {}

	TrackedGridBase (
		const VecDu& size_, const VecDi& offset_, const LeafType& background_
	):
		Base(size_, offset_, background_),
		m_grid_lookup()
	{
		this->init(size_, offset_, background_);
	}

	TrackedGridBase(const TrackedGridBase&& other) : Base(other)
	{
		m_grid_lookup = std::move(other.m_grid_lookup);
	}

	TrackedGridBase(const TrackedGridBase& other) : Base(other)
	{
		m_grid_lookup = other.m_grid_lookup;
	}

	void operator=(const TrackedGridBase& other)
	{
		Base::operator=(other);
		m_grid_lookup = other.m_grid_lookup;
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
	void size (const VecDu& size_)
	{
		Base::size(size_);
		m_grid_lookup.size(size_);
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
	PosArray& list(const UINT arr_idx_ = 0)
	{
		return m_grid_lookup.list(arr_idx_);
	}

	/**
	 * Get list of active grid points from lookup grid.
	 *
	 * @param arr_idx_ tracking list id.
	 * @return the tracking list of active grid positions from the internal lookup grid.
	 */
	const PosArray& list(const UINT arr_idx_ = 0) const
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
	bool add(const VecDi& pos_, const LeafType& val_, const UINT arr_idx_ = 0)
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
	bool add(const VecDi& pos_, const UINT arr_idx_ = 0)
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
	void reset(const LeafType& val_, const UINT arr_idx_ = 0)
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
	void reset(const UINT arr_idx_ = 0)
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
	void remove(const UINT idx_, const UINT arr_idx_ = 0)
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
	void remove (const VecDi& pos_, const UINT arr_idx_ = 0)
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
	const bool is_active (const VecDi& pos_, const UINT arr_idx_ = 0) const
	{
		return m_grid_lookup.is_active(pos_, arr_idx_);
	}
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
 * Traits of TrackedGridBase.
 *
 * Just forward the traits defined for TrackedGridBase subclasses.
 */
template <class Derived>
struct GridTraits<TrackedGridBase<Derived> > : GridTraits<Derived>
{};


/**
 * Traits of TrackedGrid.
 *
 * @tparam T the type of data to store in the grid.
 * @tparam D the dimension of the grid.
 * @tparam N the number of tracking lists to use.
 */
template <typename T, UINT D, UINT N>
struct GridTraits<TrackedGrid<T, D, N> > : DefaultGridTraits<T, D>
{
	using ThisType = TrackedGrid<T, D, N>;
	/// Type of lookup grid to use.  This is what differentiates this from SharedTrackedGrid.
	using LookupType = LookupGrid<D, N>;
};

}

#endif /* INCLUDE_FELT_TRACKEDGRID_HPP_ */
