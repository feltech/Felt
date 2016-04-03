
#ifndef INCLUDE_FELT_LOOKUPPARTITIONEDGRID_HPP_
#define INCLUDE_FELT_LOOKUPPARTITIONEDGRID_HPP_

#include "TrackingPartitionedGridBase.hpp"

namespace felt
{

/**
 * Spatially partitioned wrapper for MultiLookupGrid.
 *
 * @tparam D dimension of the grid.
 * @tparam N number of tracking lists to use.
 */
template <UINT D, UINT N=1>
class MultiLookupPartitionedGrid
	: public TrackingPartitionedGridBase<MultiLookupPartitionedGrid<D, N> >
{
public:
	using ThisType = MultiLookupPartitionedGrid<D, N>;
	using Base = TrackingPartitionedGridBase<ThisType>;
	using typename Base::VecDu;
	using typename Base::VecDi;
	using typename Base::LeafType;
	using Base::TrackingPartitionedGridBase;

	/**
	 * Construct a spatially partitioned MultiLookupGrid.
	 *
	 * Initialises grid data to NULL_IDX.
	 *
	 * @param size_ size of the grid.
	 * @param offset_ spatial offset of the grid.
	 * @param partition_size_ size of a spatial partition.
	 */
	MultiLookupPartitionedGrid (const VecDu& size_, const VecDi& offset_, const VecDu& partition_size_) :
		Base(size_, offset_, Base::Child::Traits::NULL_IDX_DATA, partition_size_)
	{}
};


/**
 * Spatially partitioned wrapper for SingleLookupGrid.
 *
 * @snippet test_PartitionedGrid.cpp LazySingleLookupPartitionedGrid initialisation
 *
 * @tparam D dimension of the grid.
 * @tparam N number of tracking lists to use.
 */
template <UINT D, UINT N>
class SingleLookupPartitionedGrid
	: public TrackingPartitionedGridBase<SingleLookupPartitionedGrid<D, N> >
{
public:
	using ThisType = SingleLookupPartitionedGrid<D, N>;
	using Base = TrackingPartitionedGridBase<ThisType>;
	using typename Base::VecDu;
	using typename Base::VecDi;

	using Base::TrackingPartitionedGridBase;

	/**
	 * Construct a spatially partitioned SingleLookupGrid.
	 *
	 * Initialises grid data to NULL_IDX.
	 *
	 * @param size_ size of the grid.
	 * @param offset_ spatial offset of the grid.
	 * @param partition_size_ size of a spatial partition.
	 */
	SingleLookupPartitionedGrid (
		const VecDu& size_, const VecDi& offset_, const VecDu& partition_size_
	) :
		Base(size_, offset_, Base::Child::Traits::NULL_IDX_DATA, partition_size_)
	{}
};


/**
 * Spatially partitioned wrapper for LazySingleLookupGrid.
 *
 * @tparam D dimension of the grid.
 * @tparam N number of tracking lists to use.
 */
template <UINT D, UINT N>
class LazySingleLookupPartitionedGrid
	: public TrackingPartitionedGridBase<LazySingleLookupPartitionedGrid<D, N> >
{
public:
	using ThisType = LazySingleLookupPartitionedGrid<D, N>;
	using Base = TrackingPartitionedGridBase<ThisType>;
	using typename Base::VecDu;
	using typename Base::VecDi;
	using typename Base::Child;

	LazySingleLookupPartitionedGrid() : Base()
	{}

	/**
	 * Construct a spatially partitioned LazySingleLookupGrid.
	 *
	 * Initialises grid data to NULL_IDX.
	 *
	 * @param size_ size of the grid.
	 * @param offset_ spatial offset of the grid.
	 * @param partition_size_ size of a spatial partition.
	 */
	LazySingleLookupPartitionedGrid (
		const VecDu& size_, const VecDi& offset_, const VecDu& partition_size_
	) : Base()
	{
		this->init(size_, offset_, partition_size_);
	}

	/**
	 * Initialise a spatially partitioned LazySingleLookupGrid.
	 *
	 * Initialises grid data to NULL_IDX.
	 *
	 * @param size_ size of the grid.
	 * @param offset_ spatial offset of the grid.
	 * @param partition_size_ size of a spatial partition.
	 */
	void init (const VecDu& size_, const VecDi& offset_, const VecDu& partition_size_)
	{
		Base::init(size_, offset_, Base::Child::Traits::NULL_IDX_DATA, partition_size_);
	}

	/**
	 * Reset and conditionally deactivate children.
	 *
	 * All child grids will be reset, but they will not be deactivated and removed from tracking
	 * if the given master grid is currently tracking them.
	 *
	 * This is an optimisation to ensure spatial partitions that are 'paired' do not get
	 * constantly created and destroyed unnecessarily.
	 *
	 * @snippet test_PartitionedGrid.cpp LazySingleLookupPartitionedGrid reset_mixed_cases
	 * @param grid_master_ `PartitionedGrid`-like grid to use as a "master"/"mask".
	 * @param list_idx_ the tracking list id to reset.
	 * @tparam Derived child class of `PartitionedGridBase`.
	 */
	template<class Derived>
	void reset(const PartitionedGridBase<Derived>& grid_master_, const UINT list_idx_)
	{
		for (const VecDi& pos_child : this->children().list(list_idx_))
		{
			Child& child = this->children().get(pos_child);
			// If the master grid is not tracking this child, then remove it from tracking under
			// this list id, potentially destroying it.
			if (!grid_master_.is_child_active(pos_child))
			{
				this->remove_child(pos_child, list_idx_);
			}
			// If the child has not been destroyed by the above, then reset as normal (loop over
			// tracking list resetting values in grid, before resizing list to 0).
			if (child.is_active())
			{
				child.reset(list_idx_);
			}
			// If the child was destroyed above, then no need to loop over grid resetting values,
			// so just reset list.
			else
			{
				child.list(list_idx_).clear();
			}
		}
	}

	/**
	 * Reset all tracking lists and data, deactivating all children except those active in given
	 * master grid.
	 *
	 * @snippet test_PartitionedGrid.cpp LazySingleLookupPartitionedGrid reset_all

	 * @param grid_master_
	 */
	template<class Derived>
	void reset_all(const PartitionedGridBase<Derived>& grid_master_)
	{
		for (UINT idx = 0; idx < this->children().lookup().NUM_LISTS; idx++)
			reset<Derived>(grid_master_, idx);
	}

	/**
	 * Add a spatial partition to children grid's tracking subgrid.
	 *
	 * Uses mutex for thread safety. Activates the child grid.
	 *
	 * @param pos_ position in spatial partition grid to track.
	 * @param arr_idx_ index of tracking list used to track position.
	 * @return true if position was added to tracking grid, false if it was already added.
	 */
	bool add_child(const VecDi& pos_, const UINT arr_idx_ = 0)
	{
		if (this->m_grid_children.is_active(pos_, arr_idx_))
			return false;
		Child& child = this->m_grid_children.get(pos_);
		std::lock_guard<std::mutex> lock(this->m_mutex_update_branch);
		if (!child.is_active())
			this->m_grid_children.get(pos_).activate();
		return this->m_grid_children.add(pos_, arr_idx_);
	}

	/**
	 * Remove a spatial partition from children grid's tracking subgrid.
	 *
	 * Uses mutex for thread safety. Deactivates the child grid if not being tracked by any list.
	 *
	 * @param pos_ position of spatial partition to stop tracking.
	 * @param arr_idx_ index of tracking list used to track position.
	 */
	void remove_child(const VecDi& pos_child_, const UINT arr_idx_ = 0)
	{
		if (!this->m_grid_children.is_active(pos_child_, arr_idx_))
			return;
		std::lock_guard<std::mutex> lock(this->m_mutex_update_branch);
		this->m_grid_children.remove(pos_child_, arr_idx_);
		if (!this->is_child_active(pos_child_))
			this->m_grid_children.get(pos_child_).deactivate();
	}
};


/**
 * Traits for MultiLookupPartitionedGrid.
 *
 * @tparam D dimension of the grid.
 * @tparam N number of tracking lists to use.
 */
template <UINT D, UINT N>
struct GridTraits<MultiLookupPartitionedGrid<D, N> > : DefaultMultiLookupGridTraits<D, N>
{
	/// The class inheriting from TrackingPartitionedGrid.
	using ThisType = MultiLookupPartitionedGrid<D, N>;
	/// Child grid class, in this case MultiLookupGrid.
	using ChildType = MultiLookupGrid<D, N>;
	/// Base grid class to partition, in this case LookupGridBase.
	using MixinType = LookupGridBase<ThisType>;
};


/**
 * Traits for SingleLookupPartitionedGrid.
 *
 * @tparam D dimension of the grid.
 * @tparam N number of tracking lists to use.
 */
template <UINT D, UINT N>
struct GridTraits<SingleLookupPartitionedGrid<D, N> > : DefaultSingleLookupGridTraits<D, N>
{
	/// The class inheriting from TrackingPartitionedGrid.
	using ThisType = SingleLookupPartitionedGrid<D, N>;
	/// Child grid class, in this case SingleLookupGrid.
	using ChildType = SingleLookupGrid<D, N>;
	/// Grid class whose interface to copy via CRTP mixin.
	using MixinType = SingleLookupGridBase<ThisType>;
};


/**
 * Traits for LazySingleLookupPartitionedGrid.
 *
 * @tparam D dimension of the grid.
 * @tparam N number of tracking lists to use.
 */
template <UINT D, UINT N>
struct GridTraits<LazySingleLookupPartitionedGrid<D, N> > : DefaultSingleLookupGridTraits<D, N>
{
	/// The class inheriting from TrackingPartitionedGrid.
	using ThisType = LazySingleLookupPartitionedGrid<D, N>;
	/// Child grid class, in this case LazySingleLookupGrid.
	using ChildType = LazySingleLookupGrid<D, N>;
	/// Grid class whose interface to copy via CRTP mixin.
	using MixinType = SingleLookupGridBase<ThisType>;
};

} // End namespace felt
#endif
