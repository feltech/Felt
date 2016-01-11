/*
 * SharedLookupGrid.hpp
 *
 *  Created on: 10 Jan 2016
 *      Author: dave
 */

#ifndef INCLUDE_FELT_SHAREDLOOKUPGRID_HPP_
#define INCLUDE_FELT_SHAREDLOOKUPGRID_HPP_

#include <array>
#include <sstream>
#include <mutex>
#include "LookupGrid.hpp"

namespace felt
{

/** @addtogroup LookupGrids
 *  @{
 */


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
	friend class LookupGridBase<SharedLookupGridBase<Derived> >;
	using ThisType = SharedLookupGridBase<Derived>;
	using Base = LookupGridBase<ThisType>;
	using Base::NULL_IDX;
	using typename Base::LeafType;
	using typename Base::VecDu;
	using typename Base::VecDi;
protected:
	using Base::m_a_pos;
public:
	using Base::LookupGridBase;

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
		const UINT& idx = idx_from_pos(pos_);
		Base::remove(idx, pos_, arr_idx_, 0);
	}

protected:

	/**
	 * Reset the entire grid to null indices.
	 */
	void clear ()
	{
		this->fill(NULL_IDX);
	}

	/**
	 * Get index in a tracking list from position.
	 *
	 * @param pos_ position in grid to find index data.
	 * @param arr_idx_ tracking list id.
	 * @return
	 */
	UINT& idx_from_pos(const VecDi& pos_, const UINT& arr_idx_ = 0)
	{
		return this->get(pos_);
	}

	/**
	 * @copydoc idx_from_pos(const VecDi&,const UINT&)
	 */
	const UINT& idx_from_pos(const VecDi& pos_, const UINT& arr_idx_) const
	{
		return this->get(pos_);
	}
};


/**
 * Traits for LookupGridBase to understand SharedLookupGridBase.
 *
 * Just forward the traits defined for SharedLookupGridBase subclasses.
 */
template <class Derived>
struct GridTraits<SharedLookupGridBase<Derived> > : GridTraits<Derived>
{};


/**
 * Standard traits for all classes CRTP derived from LookupGridBase
 *
 * @tparam D the dimension of the grid.
 * @tparam N the number of tracking lists to use.
 */
template <UINT D, UINT N>
struct DefaultSharedLookupGridTraits : DefaultGridTraits<UINT, D >
{
	/// Null index grid value in data array.
	static const UINT NULL_IDX_DATA;
	/// Number of tracking lists taken from template parameter.
	static const UINT NumLists = N;
};

template <UINT D, UINT N>
const UINT DefaultSharedLookupGridTraits<D, N>::NULL_IDX_DATA = std::numeric_limits<UINT>::max();


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
struct GridTraits<SharedLookupGrid<D, N> > : DefaultSharedLookupGridTraits<D, N>
{
	/// The derived type.
	using ThisType = SharedLookupGrid<D, N>;
};

/** @} */ // End group LookupGrid.

} // End namespace felt


#endif /* INCLUDE_FELT_SHAREDLOOKUPGRID_HPP_ */
