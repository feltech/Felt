#ifndef INCLUDE_FELT_SingleTRACKEDGrid_HPP_
#define INCLUDE_FELT_SingleTRACKEDGrid_HPP_

#include "TrackedGridBase.hpp"

namespace felt
{

/**
 * A tracked grid that assumes non-overlapping tracking lists.
 *
 * A grid of arbitrary data, with active positions tracked by an internal SingleLookupGrid.
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
class SingleTrackedGrid : public TrackedGridBase< SingleTrackedGrid<T, D, N> >
{
public:
	using ThisType = SingleTrackedGrid<T, D, N>;
	using Base = TrackedGridBase<ThisType>;
	using TrackedGridBase<ThisType>::TrackedGridBase;
};

/**
 * A lazy tracked grid that assumes non-overlapping tracking lists.
 *
 * @copydetails SingleTrackedGrid
 *
 * Lazy variant, that can be activated and deactivated (data array destroyed and created).
 *
 * @tparam T the type of data to store in the grid.
 * @tparam D the dimension of the grid.
 * @tparam N the number of tracking lists to use.
 */
template <typename T, UINT D, UINT N=1>
class LazySingleTrackedGrid : public TrackedGridBase<LazySingleTrackedGrid<T, D, N>, true>
{
public:
	using ThisType = LazySingleTrackedGrid<T, D, N>;
	using Base = TrackedGridBase<ThisType, true>;
	using Base::TrackedGridBase;
	using Base::Base::is_active;

	/**
	 * Create the internal data array and fill with background value.
	 *
	 * Also activates lookup grid.
	 *
	 * @snippet test_SingleTrackedGrid.cpp LazySingleTrackedGrid activate
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
	 * @snippet test_SingleTrackedGrid.cpp LazySingleTrackedGrid deactivate
	 */
	void deactivate()
	{
		Base::deactivate();
		this->m_grid_lookup.deactivate();
	}
};


/**
 * Traits of SingleTrackedGrid.
 *
 * @tparam T the type of data to store in the grid.
 * @tparam D the dimension of the grid.
 * @tparam N the number of tracking lists to use.
 */
template <typename T, UINT D, UINT N>
struct GridTraits<SingleTrackedGrid<T, D, N> > : DefaultGridTraits<T, D>
{
	using ThisType = SingleTrackedGrid<T, D, N>;
	/// Type of lookup grid to use.  This is what differentiates this from MultiTrackedGrid.
	using MultiLookupType = SingleLookupGrid<D, N>;
};


/**
 * Traits of LazySingleTrackedGrid.
 *
 * @tparam T the type of data to store in the grid.
 * @tparam D the dimension of the grid.
 * @tparam N the number of tracking lists to use.
 */
template <typename T, UINT D, UINT N>
struct GridTraits<LazySingleTrackedGrid<T, D, N> > : DefaultGridTraits<T, D>
{
	using ThisType = LazySingleTrackedGrid<T, D, N>;
	/// Type of lookup grid to use.  This is what differentiates this from MultiTrackedGrid.
	using MultiLookupType = LazySingleLookupGrid<D, N>;
};

}
#endif /* INCLUDE_FELT_SingleTRACKEDGrid_HPP_ */
