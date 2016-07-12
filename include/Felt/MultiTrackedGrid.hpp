#ifndef INCLUDE_FELT_MULTITRACKEDGRID_HPP_
#define INCLUDE_FELT_MULTITRACKEDGRID_HPP_

#include "TrackedGridBase.hpp"

namespace felt
{

/**
 * Standard tracked grid.
 *
 * A grid of arbitrary data, with active positions tracked by an internal MultiLookupGrid.
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
class MultiTrackedGrid : public TrackedGridBase< MultiTrackedGrid<T, D, N>, Laziness::EAGER >
{
public:
	using TrackedGridBase<MultiTrackedGrid<T, D, N>, Laziness::EAGER >::TrackedGridBase;
};


/**
 * Traits of MultiTrackedGrid.
 *
 * @tparam T the type of data to store in the grid.
 * @tparam D the dimension of the grid.
 * @tparam N the number of tracking lists to use.
 */
template <typename T, UINT D, UINT N>
struct GridTraits<MultiTrackedGrid<T, D, N> > : DefaultGridTraits<T, D>
{
	using ThisType = MultiTrackedGrid<T, D, N>;
	/// Type of lookup grid to use.  This is what differentiates this from SingleTrackedGrid.
	using LookupType = MultiLookupGrid<D, N>;
};

}

#endif /* INCLUDE_FELT_MULTITRACKEDGRID_HPP_ */
