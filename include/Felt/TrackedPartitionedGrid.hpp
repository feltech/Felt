#ifndef INCLUDE_FELT_TRACKEDPARTITIONEDGRID_HPP_
#define INCLUDE_FELT_TRACKEDPARTITIONEDGRID_HPP_

#include "TrackingPartitionedGridBase.hpp"

namespace felt
{

/**
 * Spatially partitioned wrapper for TrackedGrid.
 *
 * @tparam T type of value stored in leaf grid nodes.
 * @tparam D dimension of the grid.
 * @tparam N number of tracking lists to use.
 */
template <typename T, UINT D, UINT N>
class TrackedPartitionedGrid
	: public TrackingPartitionedGridBase <TrackedPartitionedGrid<T, D, N>>
{
public:
	/// This class.
	using ThisType = TrackedPartitionedGrid<T, D, N>;
	/// Base class
	using Base = TrackingPartitionedGridBase<ThisType>;

	using typename Base::Child;
	using typename Base::VecDu;
	using typename Base::VecDi;
	using typename Base::ChildrenGrid;

	using PosArray = typename Child::PosArray;
public:
	using Base::TrackingPartitionedGridBase;

	/**
	 * Set value in grid at given position and add position to lookup grid.
	 *
	 * Will set value regardless whether lookup grid already set for given
	 * position + tracking list.
	 *
	 * Descend to relevant child grid to add to their tracking structure.
	 *
	 * @param pos_ position in grid.
	 * @param val_ value to set.
	 * @param arr_idx_ tracking list id.
	 * @return true if grid node set in child lookup grid and position added
	 * to tracking list, false if child grid node was already set so
	 * position already in a list.
	 */
	bool add(const VecDi& pos_, const T& val_, const UINT arr_idx_ = 0)
	{
		const VecDi& pos_child = this->pos_child(pos_);
		Base::add_child(pos_child, arr_idx_);
		return this->children().get(pos_child).add(pos_, val_, arr_idx_);
	}

	/**
	 * Thread-safely set value in grid at given position and add position to lookup grid.
	 *
	 * Will set value regardless whether lookup grid already set for given
	 * position + tracking list.
	 *
	 * Descend to relevant child grid to add to their tracking structure.
	 *
	 * @param pos_ position in grid.
	 * @param val_ value to set.
	 * @param arr_idx_ tracking list id.
	 * @return true if grid node set in child lookup grid and position added
	 * to tracking list, false if child grid node was already set so
	 * position already in a list.
	 */
	bool add_safe(const VecDi& pos_, const T& val_, const UINT arr_idx_ = 0)
	{
		const VecDi& pos_child = this->pos_child(pos_);
		Base::add_child(pos_child, arr_idx_);
		Child& child = this->children().get(pos_child);
		std::lock_guard<std::mutex> lock(child.mutex());
		return child.add(pos_, val_, arr_idx_);
	}

	/**
	 * Set every active grid node (i.e. those referenced by lookup grid)
	 * to given value and reset the lookup grid.
	 *
	 * Lookup grid will then be full of NULL indices for given tracking list
	 * and the relevant tracking list(s) will be empty.
	 *
	 * @param val_ value to set in main grid.
	 * @param arr_idx_ tracking list id to cycle over and clear.
	 */
	void reset(const T& val_, const UINT arr_idx_ = 0)
	{
		ChildrenGrid& children = this->children();
		for (const VecDi& pos_child : children.list(arr_idx_))
			children(pos_child).reset(val_, arr_idx_);

		Base::Base::reset(arr_idx_);
	}
};


/**
 * Spatially partitioned wrapper for SharedTrackedGrid.
 *
 * @tparam T type of value stored in leaf grid nodes.
 * @tparam D dimension of the grid.
 * @tparam N number of tracking lists to use.
 */
template <typename T, UINT D, UINT N>
class SharedTrackedPartitionedGrid
	: public TrackingPartitionedGridBase<SharedTrackedPartitionedGrid<T, D, N> >
{
public:
	using ThisType = SharedTrackedPartitionedGrid<T, D, N>;
	using Base = TrackingPartitionedGridBase<ThisType>;
	using typename Base::Child;
	using typename Base::VecDu;
	using typename Base::VecDi;
	using typename Base::ChildrenGrid;
	using PosArray = typename Child::PosArray;

public:
	using Base::TrackingPartitionedGridBase;
	using Base::add;

	/**
	 * Set value in grid at given position and add position to lookup grid.
	 *
	 * Will set value regardless whether lookup grid already set for given
	 * position + tracking list.
	 *
	 * Descend to relevant child grid to add to their tracking structure.
	 *
	 * @param pos_ position in grid.
	 * @param val_ value to set.
	 * @param arr_idx_ tracking list id.
	 * @return true if grid node set in child lookup grid and position added
	 * to tracking list, false if child grid node was already set so
	 * position already in a list.
	 */
	bool add(const VecDi& pos_, const T& val_, const UINT arr_idx_)
	{
		const VecDi& pos_child = this->pos_child(pos_);
		Base::add_child(pos_child, arr_idx_);
		return this->children().get(pos_child).add(pos_, val_, arr_idx_);
	}

	/**
	 * Set every active grid node (i.e. those referenced by lookup grid)
	 * to given value and reset the lookup grid.
	 *
	 * Lookup grid will then be full of NULL indices for given tracking list
	 * and the relevant tracking list(s) will be empty.
	 *
	 * @param val_ value to set in main grid.
	 * @param arr_idx_ tracking list id to cycle over and clear.
	 */
	void reset(const T& val_, const UINT arr_idx_ = 0)
	{
		ChildrenGrid& children = this->children();
		for (const VecDi& pos_child : children.list(arr_idx_))
			children(pos_child).reset(val_, arr_idx_);

		Base::Base::reset(arr_idx_);
	}
};


/**
 * Traits for TrackedPartitionedGrid.
 *
 * @tparam T type of value stored in leaf grid nodes.
 * @tparam D dimension of the grid.
 * @tparam N number of tracking lists to use.
 */
template <typename T, UINT D, UINT N>
struct GridTraits<TrackedPartitionedGrid<T, D, N> > : DefaultGridTraits<T, D>
{
	/// The class inheriting from the base.
	using ThisType = TrackedPartitionedGrid<T, D, N>;
	/// The type of lookup grid to use for tracking active grid nodes.
	using LookupType = LookupGrid<D, N>;
	/// Child grid class, in this case a TrackedGrid.
	using ChildType = TrackedGrid<T, D, N>;
	/// Base grid class to partition, in this case TrackedGridBase.
	using MixinType = TrackedGridBase<ThisType>;
	/// Number of tracking lists, from template parameter.
	static const UINT NumLists = N;
};


/**
 * Traits for SharedTrackedPartitionedGrid.
 *
 * @tparam T type of value stored in leaf grid nodes.
 * @tparam D dimension of the grid.
 * @tparam N number of tracking lists to use.
 */
template <typename T, UINT D, UINT N>
struct GridTraits<SharedTrackedPartitionedGrid<T, D, N> > : DefaultGridTraits<T, D>
{
	/// The class inheriting from TrackedGridBase.
	using ThisType = SharedTrackedPartitionedGrid<T, D, N>;
	/// The type of lookup grid to use for tracking active grid nodes.
	using LookupType = SharedLookupGrid<D, N>;
	/// Child grid class, in this case SharedTrackedGrid.
	using ChildType = SharedTrackedGrid<T, D, N>;
	/// Base grid class to partition, in this case TrackedGridBase.
	using MixinType = TrackedGridBase<ThisType>;
	/// Number of tracking lists, from template parameter.
	static const UINT NumLists = N;
};

} // End namespace felt.
#endif /* INCLUDE_FELT_TRACKEDPARTITIONEDGRID_HPP_ */
