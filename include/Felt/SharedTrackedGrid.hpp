#ifndef INCLUDE_FELT_SHAREDTRACKEDGRID_HPP_
#define INCLUDE_FELT_SHAREDTRACKEDGRID_HPP_

#include "TrackedGrid.hpp"

namespace felt
{

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
	using ThisType = SharedTrackedGrid<T, D, N>;
	using Base = TrackedGridBase<ThisType>;
	using TrackedGridBase<ThisType>::TrackedGridBase;
};

/**
 * A lazy tracked grid that assumes non-overlapping tracking lists.
 *
 * @copydetails SharedTrackedGrid
 *
 * Lazy variant, that can be activated and deactivated (data array destroyed and created).
 *
 * @tparam T the type of data to store in the grid.
 * @tparam D the dimension of the grid.
 * @tparam N the number of tracking lists to use.
 */
template <typename T, UINT D, UINT N=1>
class LazySharedTrackedGrid : public TrackedGridBase<LazySharedTrackedGrid<T, D, N>, true>
{
public:
	using ThisType = LazySharedTrackedGrid<T, D, N>;
	using Base = TrackedGridBase<ThisType, true>;
	using Base::TrackedGridBase;
	using Base::Base::is_active;

	/**
	 * Create the internal data array and fill with background value.
	 *
	 * Also activates lookup grid.
	 *
	 * @snippet test_MappedGrid.cpp LazySharedTrackedGrid activate
	 */
	void activate()
	{
		Base::activate();
		this->m_grid_lookup.activate();
	}

	/**
	 * Destroy the internal data array.
	 *
	 * Also destroys lookup grid.
	 *
	 * @snippet test_MappedGrid.cpp LazySharedTrackedGrid deactivate
	 */
	void deactivate()
	{
		Base::deactivate();
		this->m_grid_lookup.deactivate();
	}
};


/**
 * Traits of SharedTrackedGrid.
 *
 * @tparam T the type of data to store in the grid.
 * @tparam D the dimension of the grid.
 * @tparam N the number of tracking lists to use.
 */
template <typename T, UINT D, UINT N>
struct GridTraits<SharedTrackedGrid<T, D, N> > : DefaultGridTraits<T, D>
{
	using ThisType = SharedTrackedGrid<T, D, N>;
	/// Type of lookup grid to use.  This is what differentiates this from TrackedGrid.
	using LookupType = SharedLookupGrid<D, N>;
};


/**
 * Traits of LazySharedTrackedGrid.
 *
 * @tparam T the type of data to store in the grid.
 * @tparam D the dimension of the grid.
 * @tparam N the number of tracking lists to use.
 */
template <typename T, UINT D, UINT N>
struct GridTraits<LazySharedTrackedGrid<T, D, N> > : DefaultGridTraits<T, D>
{
	using ThisType = LazySharedTrackedGrid<T, D, N>;
	/// Type of lookup grid to use.  This is what differentiates this from TrackedGrid.
	using LookupType = LazySharedLookupGrid<D, N>;
};

}
#endif /* INCLUDE_FELT_SHAREDTRACKEDGRID_HPP_ */
