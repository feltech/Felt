#ifndef INCLUDE_FELT_TRACKINGPARTITIONEDGRIDBASE_HPP_
#define INCLUDE_FELT_TRACKINGPARTITIONEDGRIDBASE_HPP_

#include "PartitionedGrid.hpp"

namespace felt
{

/**
 * Container wrapping iterator through leafs of a partitioned grid.
 *
 * @tparam G grid type to iterate over.
 */
template <typename G>
class LeafsContainer
{
private:
	/// Partitioned grid type to iterate over.
	using GridTree = G;

	using VecDi = typename GridTree::VecDi;
	using PosArray = typename GridTree::PosArray;

	/// Pointer to partitioned grid to iterate over
	const GridTree* m_pgrid;
	/// Index of tracking list of active points to iterate over.
	const UINT	m_list_idx;

public:
	/**
	 * Iterator class for range-based for loops across partitioned grid leafs.
	 */
	class iterator : public boost::iterator_facade<
		LeafsContainer::iterator,
		const VecDi, boost::forward_traversal_tag
	>
	{
	private:
		using Iter = typename PosArray::const_iterator;
	public:
		iterator() : m_pgrid(NULL), m_it_child_end()
		{}

		/**
		 * Construct an iterator over leafs of a partitioned grid.
		 *
		 * @param pgrid grid to iterate through
		 * @param listIdx tracking list id within grid
		 * @param it_child iterator over child grids
		 * @param it_leaf iterator over lowest level leaf grid point.
		 */
		iterator(
			const GridTree* pgrid, const UINT list_idx,
			const Iter& it_child, const Iter& it_leaf
		)
		: m_pgrid(pgrid), m_list_idx(list_idx), m_it_child(it_child),
		  m_it_leaf(it_leaf),
		  m_it_child_end(pgrid->children().list(list_idx).end())
		{}

	private:
		friend class boost::iterator_core_access;

		/// Pointer to partitioned grid to iterate over.
		const GridTree* m_pgrid;
		/// Index of tracking list to iterate through.
		const UINT	m_list_idx;
		/// Iterator pointing to current spatial partition position.
		Iter		m_it_child;
		/// Iterator pointing to current leaf grid node position.
		Iter		m_it_leaf;
		/// Iterator pointing past end of final spatial partition grid position.
		const Iter	m_it_child_end;

		/**
		 * Override specifying how to move on to next leaf, jumping from (active) child to child.
		 */
		void increment()
		{
			if (m_it_child == m_it_child_end)
				return;

			if (m_it_leaf == m_it_child)
			{
				m_it_leaf = m_pgrid->children().get(*m_it_child).list(m_list_idx).begin();
				return;
			}

			m_it_leaf++;

			if (m_it_leaf == m_pgrid->children().get(*m_it_child).list(m_list_idx).end())
			{
				m_it_child++;
				m_it_leaf = m_it_child;
				increment();
				return;
			}
		}

		/**
		 * Check for equality between this iterator and another.
		 *
		 * @param other_ iterator to compare against.
		 * @return true if equal, false if not equal.
		 */
		bool equal(const iterator& other_) const
		{
			return (
				m_it_child == other_.m_it_child
				&& m_it_leaf == other_.m_it_leaf
			);
		}

		/**
		 * Dereference iterator into grid position.
		 *
		 * @return leaf grid node position currently pointed to.
		 */
		const VecDi& dereference() const {
			const VecDi& pos = *m_it_leaf;
			return *m_it_leaf;
		}
	};

	/**
	 * Construct a wrapper for range-based for loops over active partitioned grid nodes.
	 *
	 * @param pgrid grid to iterate over
	 * @param list_idx tracking list id identifying leafs
	 */
	LeafsContainer(const GridTree* pgrid, const UINT list_idx)
	: m_pgrid(pgrid), m_list_idx(list_idx)
	{}

	/**
	 * Get first iterator for leafs identified within list.
	 *
	 * @return an iterator to the first leaf in the tracking list.
	 */
	const iterator begin() const
	{
		const typename PosArray::const_iterator& it_child_begin
			= m_pgrid->children().list(m_list_idx).cbegin();
		const typename PosArray::const_iterator& it_child_end
			= m_pgrid->children().list(m_list_idx).cend();

		typename PosArray::const_iterator it_leaf_begin = it_child_end;

		if (it_child_begin != it_child_end)
		{
			it_leaf_begin = m_pgrid->children().get(
				*it_child_begin
			).list(m_list_idx).begin();
		}

		return iterator(
			m_pgrid, m_list_idx,
			it_child_begin,
			it_leaf_begin
		);
	}

	/**
	 * Iterator representing one element past the end of the list of leafs.
	 *
	 * @return an iterator pointing past the end of tracking lists.
	 */
	const iterator end() const
	{
		return iterator(
			m_pgrid, m_list_idx,
			m_pgrid->children().list(m_list_idx).cend(),
			m_pgrid->children().list(m_list_idx).cend()
		);
	}

	/**
	 * Calculate length of list by summing lists in all partitions.
	 *
	 * @return sum of size of lists in all spatial partitions.
	 */
	UINT size() const
	{
		UINT sum = 0;
		for (const VecDi& pos_child : m_pgrid->children().list(m_list_idx))
			sum += m_pgrid->children().get(pos_child).list(m_list_idx).size();
		return sum;
	}
};


template <class Derived, bool IsLazy=false>
class TrackingPartitionedGridBase
{};


/**
 * Base class for spatially partitioned wrappers for MultiLookupGrid and MultiTrackedGrid.
 *
 * @tparam Derived the CRTP derived class.
 */
template <class Derived>
class TrackingPartitionedGridBase<Derived, false>
	: public PartitionedGridBase<TrackingPartitionedGridBase<Derived> >
{
public:
	using ThisType = TrackingPartitionedGridBase<Derived>;
	using DerivedType = typename GridTraits<Derived>::ThisType;
	/// Base class.
	using Base = PartitionedGridBase<ThisType>;
	/// Allow full access to iterator over leaf grid nodes.
	friend class LeafsContainer<ThisType>;

	using typename Base::Child;
	using typename Base::VecDu;
	using typename Base::VecDi;
	using typename Base::ChildrenGrid;

	using PosArray = typename Child::PosArray;

public:

	using Base::PartitionedGridBase;

	/**
	 * Reset the grid nodes referenced in tracking list.
	 *
	 * Descend to children to reset their tracking list.
	 *
	 * @param arr_idx_ tracking list id.
	 */
	void reset(const UINT arr_idx_ = 0)
	{
		ChildrenGrid& children = this->children();
		for (const VecDi& pos_child : children.list(arr_idx_))
			children(pos_child).reset(arr_idx_);

		Base::reset(arr_idx_);
	}

	/**
	 * Add a leaf position to be tracked to given tracking list.
	 *
	 * Descend to relevant child grid to add to their tracking structure.
	 *
	 * @param pos_ position to add.
	 * @param arr_idx_ tracking list id.
	 * @return true if grid node set in child lookup grid and position added
	 * to tracking list, false if child grid node was already set so
	 * position already in a list.
	 */
	bool add(const VecDi& pos_, const UINT arr_idx_ = 0)
	{
		const VecDi& pos_child = this->pos_child(pos_);
		self->add_child(pos_child, arr_idx_);
		return this->children().get(pos_child).add(pos_, arr_idx_);
	}

	/**
	 * Thread safely add a leaf position to be tracked to given tracking
	 * list.
	 *
	 * Descend to relevant child grid to add to their tracking structure.
	 *
	 * This is a safe version that uses a mutex lock on the child grid,
	 * which is necessary if threads can "cross over" to other partitions.
	 *
	 * @param pos_ position to add.
	 * @param arr_idx_ tracking list id.
	 * @return true if grid node set in child lookup grid and position added
	 * to tracking list, false if child grid node was already set so
	 * position already in a list.
	 */
	bool add_safe(const VecDi& pos_, const UINT arr_idx_ = 0)
	{
		const VecDi& pos_child = this->pos_child(pos_);
		self->add_child(pos_child, arr_idx_);
		Child& child = this->children().get(pos_child);
		if (child.is_active(pos_))
			return false;
		std::lock_guard<std::mutex> lock(child.mutex());
		return child.add(pos_, arr_idx_);
	}

	/**
	 * Remove a leaf position from relevant child tracking structure and
	 * remove child from tracking list if child's list is now empty.
	 *
	 * @param pos_ leaf position to remove.
	 * @param arr_idx_ tracking list id.
	 */
	void remove(const VecDi& pos_, const UINT arr_idx_ = 0)
	{
		const VecDi& pos_child = this->pos_child(pos_);
		Child& child = this->children().get(pos_child);
		child.remove(pos_, arr_idx_);
		if (child.list(arr_idx_).size() == 0)
			self->remove_child(pos_child, arr_idx_);
	}

	/**
	 * Return structure for range based for loops over leaf nodes.
	 *
	 * @param list_idx_ tracking list id.
	 * @return
	 */
	const LeafsContainer<DerivedType> leafs(const UINT list_idx_ = 0) const
	{
		return LeafsContainer<DerivedType>(
			static_cast<const DerivedType*>(this), list_idx_
		);
	}

private:
	/**
	 * Override spoofed (non-partitioned) base class's list method to make
	 * private, since it will always be empty (must use children or child).
	 *
	 * @param arr_idx_
	 * @return list of this grid (will always be empty).
	 */
	PosArray& list(const UINT arr_idx_ = 0)
	{
		return Child::list(arr_idx_);
	}
};


/**
 * Lazy variant of TrackingPartitionedGridBase.
 */
template <class Derived>
class TrackingPartitionedGridBase<Derived, true>
	: public TrackingPartitionedGridBase< TrackingPartitionedGridBase<Derived, true>, false >
{
public:
	using ThisType = TrackingPartitionedGridBase<Derived, true>;
	using Base = TrackingPartitionedGridBase<ThisType, false>;
	using typename Base::Child;
	using typename Base::VecDi;
	using Base::TrackingPartitionedGridBase;

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
	template<class MasterDerived>
	void reset(const PartitionedGridBase<MasterDerived>& grid_master_, const UINT list_idx_)
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
	template<class MasterDerived>
	void reset_all(const PartitionedGridBase<MasterDerived>& grid_master_)
	{
		for (UINT idx = 0; idx < this->children().lookup().NUM_LISTS; idx++)
			this->reset(grid_master_, idx);
	}
};


/**
 * Traits for TrackingPartitionedGridBase.
 *
 * Just forward the traits defined for TrackingPartitionedGridBase subclasses.
 *
 * @tparam Derived the CRTP derived class.
 */
template <class Derived, bool IsLazy>
struct GridTraits<TrackingPartitionedGridBase<Derived, IsLazy> > : GridTraits<Derived>
{};


} // End namespace felt.
#endif /* INCLUDE_FELT_TRACKINGPARTITIONEDGRIDBASE_HPP_ */
