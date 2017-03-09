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
class TrackedGrid : public TrackedGridBase< TrackedGrid<T, D, N> >
{
public:
	using ThisType = TrackedGrid<T, D, N>;
	using Base = TrackedGridBase<ThisType>;
	using Base::TrackedGridBase;
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
class LazyTrackedGrid : public TrackedGridBase<LazyTrackedGrid<T, D, N>>
{
public:
	using ThisType = LazyTrackedGrid<T, D, N>;
	using Traits = GridTraits<ThisType>;
	using Base = TrackedGridBase<ThisType>;
	using typename Base::VecDi;

public:
	using Base::TrackedGridBase;

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
	/// Type of lookup grid to use.  This is what differentiates this from MultiTrackedGrid.
	using LookupType = LookupGrid<D, N>;
	/// Set as eagerly initialised.
	static const Laziness IsLazy = Laziness::EAGER;
};


/**
 * Traits of LazyTrackedGrid.
 *
 * @tparam T the type of data to store in the grid.
 * @tparam D the dimension of the grid.
 * @tparam N the number of tracking lists to use.
 */
template <typename T, UINT D, UINT N>
struct GridTraits<LazyTrackedGrid<T, D, N> > : DefaultGridTraits<T, D>
{
	using ThisType = LazyTrackedGrid<T, D, N>;
	/// Type of lookup grid to use.  This is what differentiates this from MultiTrackedGrid.
	using LookupType = LazyLookupGrid<D, N>;
	/// Set as lazily initialised.
	static const Laziness IsLazy = Laziness::LAZY;
};

}
#endif /* INCLUDE_FELT_SingleTRACKEDGrid_HPP_ */
