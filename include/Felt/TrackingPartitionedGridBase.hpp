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


/**
 * Base class for spatially partitioned wrappers for LookupGrid and TrackedGrid.
 *
 * @tparam Derived the CRTP derived class.
 */
template <class Derived>
class TrackingPartitionedGridBase
	: public PartitionedGridBase<TrackingPartitionedGridBase<Derived> >
{
public:
	using ThisType = TrackingPartitionedGridBase<Derived>;
	using DerivedType = typename GridTraits<Derived>::ThisType;
	/// Base class.
	using Base = PartitionedGridBase<ThisType>;

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
	const LeafsContainer<Derived> leafs(const UINT list_idx_ = 0) const
	{
		return LeafsContainer<Derived>(
			static_cast<const Derived*>(this), list_idx_
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
 * Traits for TrackingPartitionedGridBase.
 *
 * Just forward the traits defined for TrackingPartitionedGridBase subclasses.
 *
 * @tparam Derived the CRTP derived class.
 */
template <class Derived>
struct GridTraits<TrackingPartitionedGridBase<Derived> > : GridTraits<Derived>
{};

} // End namespace felt.
#endif /* INCLUDE_FELT_TRACKINGPARTITIONEDGRIDBASE_HPP_ */
