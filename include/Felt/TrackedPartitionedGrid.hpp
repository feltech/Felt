#ifndef INCLUDE_FELT_TRACKEDPARTITIONEDGRID_HPP_
#define INCLUDE_FELT_TRACKEDPARTITIONEDGRID_HPP_

#include "TrackingPartitionedGridBase.hpp"
#include "SingleTrackedGrid.hpp"

namespace felt
{

template <class Derived, bool IsLazy=false>
class TrackedPartitionedGridBase
	: public TrackingPartitionedGridBase <TrackedPartitionedGridBase<Derived>, IsLazy>
{
public:
	/// This class.
	using ThisType = TrackedPartitionedGridBase<Derived>;
	/// Base class
	using Base = TrackingPartitionedGridBase<ThisType, IsLazy>;
	using DerivedType = Derived;

	using typename Base::Child;
	using typename Base::VecDu;
	using typename Base::VecDi;
	using typename Base::ChildrenGrid;
	using LeafType = typename GridTraits<DerivedType>::LeafType;

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
	 * @param list_idx_ tracking list id.
	 * @return true if grid node set in child lookup grid and position added
	 * to tracking list, false if child grid node was already set so
	 * position already in a list.
	 */
	bool add(const VecDi& pos_, const LeafType& val_, const UINT list_idx_)
	{
		const VecDi& pos_child = this->pos_child(pos_);
		nself->add_child(pos_child, list_idx_);
		return this->children().get(pos_child).add(pos_, val_, list_idx_);
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
	 * @param list_idx_ tracking list id.
	 * @return true if grid node set in child lookup grid and position added
	 * to tracking list, false if child grid node was already set so
	 * position already in a list.
	 */
	bool add_safe(const VecDi& pos_, const LeafType& val_, const UINT list_idx_)
	{
		const VecDi& pos_child = this->pos_child(pos_);
		nself->add_child(pos_child, list_idx_);
		Child& child = this->children().get(pos_child);
		std::lock_guard<std::mutex> lock(child.mutex());
		return child.add(pos_, val_, list_idx_);
	}

	/**
	 * Set every active grid node (i.e. those referenced by lookup grid)
	 * to given value and reset the lookup grid.
	 *
	 * MultiLookup grid will then be full of NULL indices for given tracking list
	 * and the relevant tracking list(s) will be empty.
	 *
	 * @param val_ value to set in main grid.
	 * @param list_idx_ tracking list id to cycle over and clear.
	 */
	void reset(const LeafType& val_, const UINT list_idx_ = 0)
	{
		ChildrenGrid& children = this->children();
		for (const VecDi& pos_child : children.list(list_idx_))
			children(pos_child).reset(val_, list_idx_);

		Base::Base::reset(list_idx_);
	}
};


/**
 * Spatially partitioned wrapper for MultiTrackedGrid.
 *
 * @tparam T type of value stored in leaf grid nodes.
 * @tparam D dimension of the grid.
 * @tparam N number of tracking lists to use.
 */
template <typename T, UINT D, UINT N>
class MultiTrackedPartitionedGrid
	: public TrackedPartitionedGridBase <MultiTrackedPartitionedGrid<T, D, N>>
{
public:
	/// This class.
	using ThisType = MultiTrackedPartitionedGrid<T, D, N>;
	/// Base class
	using Base = TrackedPartitionedGridBase<ThisType>;
	using Traits = GridTraits<ThisType>;
	using LeafType = typename Traits::LeafType;
	using typename Base::Child;
	using typename Base::VecDi;
	using typename Base::ChildrenGrid;

	using PosArray = typename Child::PosArray;
public:
	using Base::TrackedPartitionedGridBase;
};


/**
 * Spatially partitioned wrapper for SingleTrackedGrid.
 *
 * @tparam T type of value stored in leaf grid nodes.
 * @tparam D dimension of the grid.
 * @tparam N number of tracking lists to use.
 */
template <typename T, UINT D, UINT N>
class SingleTrackedPartitionedGrid
	: public TrackedPartitionedGridBase<SingleTrackedPartitionedGrid<T, D, N> >
{
public:
	using ThisType = SingleTrackedPartitionedGrid<T, D, N>;
	using Base = TrackedPartitionedGridBase<ThisType>;
	using Traits = GridTraits<ThisType>;

	using Base::TrackedPartitionedGridBase;
 };


/**
 * Spatially partitioned wrapper for LazySingleTrackedGrid.
 *
 * @tparam T type of value stored in leaf grid nodes.
 * @tparam D dimension of the grid.
 * @tparam N number of tracking lists to use.
 */
template <typename T, UINT D, UINT N>
class LazySingleTrackedPartitionedGrid
	: public TrackedPartitionedGridBase<LazySingleTrackedPartitionedGrid<T, D, N>, true>
{
public:
	using ThisType = LazySingleTrackedPartitionedGrid<T, D, N>;
	using Base = TrackedPartitionedGridBase<ThisType, true>;
	using Traits = GridTraits<ThisType>;
	using LeafType = typename Traits::LeafType;
	using DerivedType = typename Traits::ThisType;
	using typename Base::VecDi;
	using typename Base::Child;
public:
	using Base::TrackedPartitionedGridBase;
	using Base::Base::reset;
	using Base::remove;
	using Base::remove_child;

	/**
	 * Remove a leaf position from relevant child tracking structure and
	 * remove child from tracking list if child's list is now empty.
	 *
	 * @param pos_ leaf position to remove.
	 * @param list_idx_ tracking list id.
	 */
	void remove(const VecDi& pos_, const UINT list_idx_, const LeafType& background_)
	{
		const VecDi& pos_child = this->pos_child(pos_);
		Child& child = this->children().get(pos_child);
		child.remove(pos_, list_idx_);
		if (child.list(list_idx_).size() == 0)
			remove_child(pos_child, list_idx_, background_);
	}

	/**
	 * Move a tracked point from one tracking list to another.
	 *
	 * @param pos_ leaf position to remove.
	 * @param list_idx_ tracking list id.
	 */
	void move(
		const VecDi& pos_, const UINT list_idx_from_, const UINT list_idx_to_
	) {
		const VecDi& pos_child = this->pos_child(pos_);

		Child& child = this->children().get(pos_child);
		child.remove(pos_, list_idx_from_);
		child.add(pos_, list_idx_to_);

		std::lock_guard<std::mutex> lock(this->m_mutex_update_branch);
		nself->children().add(pos_child, list_idx_to_);
		if (child.list(list_idx_from_).size() == 0)
			nself->children().remove(pos_child, list_idx_from_);
	}


	/**
	 * Remove a spatial partition from children grid's tracking subgrid.
	 *
	 * Uses mutex for thread safety. Deactivates the child grid if not being tracked by any list.
	 *
	 * @param pos_ position of spatial partition to stop tracking.
	 * @param list_idx_ index of tracking list used to track position.
	 */
	void remove_child(
		const VecDi& pos_child_, const UINT list_idx_, const LeafType& background_
	) {
		if (!this->m_grid_children.is_active(pos_child_, list_idx_))
			return;
		std::lock_guard<std::mutex> lock(this->m_mutex_update_branch);
		this->m_grid_children.remove(pos_child_, list_idx_);
		if (!this->is_child_active(pos_child_))
		{
			Child& child = this->m_grid_children.get(pos_child_);
			child.background() = background_;
			this->m_grid_children.get(pos_child_).deactivate();
		}
	}
 };


/**
 * Traits for common base class of partitioned MultiTrackedGrids.
 */
template <class Derived>
struct GridTraits< TrackedPartitionedGridBase<Derived> > :  GridTraits<Derived> {};


/**
 * Traits for MultiTrackedPartitionedGrid.
 *
 * @tparam T type of value stored in leaf grid nodes.
 * @tparam D dimension of the grid.
 * @tparam N number of tracking lists to use.
 */
template <typename T, UINT D, UINT N>
struct GridTraits<MultiTrackedPartitionedGrid<T, D, N> > : DefaultGridTraits<T, D>
{
	/// The class inheriting from the base.
	using ThisType = MultiTrackedPartitionedGrid<T, D, N>;
	/// The type of lookup grid to use for tracking active grid nodes.
	using LookupType = MultiLookupGrid<D, N>;
	/// Child grid class, in this case a MultiTrackedGrid.
	using ChildType = MultiTrackedGrid<T, D, N>;
	/// Base grid class to partition, in this case TrackedGridBase.
	using MixinType = TrackedGridBase<ThisType>;
	/// Number of tracking lists, from template parameter.
	static const UINT NumLists = N;
};


/**
 * Traits for SingleTrackedPartitionedGrid.
 *
 * @tparam T type of value stored in leaf grid nodes.
 * @tparam D dimension of the grid.
 * @tparam N number of tracking lists to use.
 */
template <typename T, UINT D, UINT N>
struct GridTraits<SingleTrackedPartitionedGrid<T, D, N> > : DefaultGridTraits<T, D>
{
	/// The class inheriting from the base.
	using ThisType = SingleTrackedPartitionedGrid<T, D, N>;
	/// The type of lookup grid to use for tracking active grid nodes.
	using LookupType = SingleLookupGrid<D, N>;
	/// Child grid class, in this case SingleTrackedGrid.
	using ChildType = SingleTrackedGrid<T, D, N>;
	/// Base grid class to partition, in this case TrackedGridBase.
	using MixinType = TrackedGridBase<ThisType>;
	/// Number of tracking lists, from template parameter.
	static const UINT NumLists = N;
};


/**
 * Traits for LazySingleTrackedPartitionedGrid.
 *
 * @tparam T type of value stored in leaf grid nodes.
 * @tparam D dimension of the grid.
 * @tparam N number of tracking lists to use.
 */
template <typename T, UINT D, UINT N>
struct GridTraits<LazySingleTrackedPartitionedGrid<T, D, N> > : DefaultGridTraits<T, D>
{
	/// The class inheriting from the base.
	using ThisType = LazySingleTrackedPartitionedGrid<T, D, N>;
	/// Child grid class, in this case LazySingleTrackedGrid.
	using ChildType = LazySingleTrackedGrid<T, D, N>;
	/// The type of lookup grid to use for tracking active grid nodes.
	using LookupType = typename ChildType::Traits::LookupType;
	/// Base grid class to partition, in this case TrackedGridBase.
	using MixinType = TrackedGridBase<ThisType>;
	/// Number of tracking lists, from template parameter.
	static const UINT NumLists = N;
};

} // End namespace felt.
#endif /* INCLUDE_FELT_TRACKEDPARTITIONEDGRID_HPP_ */
