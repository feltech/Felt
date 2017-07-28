#ifndef INCLUDE_FELT_IMPL_MIXIN_PARTITIONEDMIXIN_HPP_
#define INCLUDE_FELT_IMPL_MIXIN_PARTITIONEDMIXIN_HPP_
#include <memory>
#include <mutex>

#include <Felt/Util.hpp>
#include <Felt/Impl/Grid.hpp>
#include <Felt/Impl/Tracked.hpp>

namespace Felt
{
namespace Impl
{
namespace Mixin
{
namespace Partitioned
{

template <class Derived>
class Children
{
private:
	/// Traits of derived class.
	using TraitsType = Traits<Derived>;
	/// Child grid type.
	using ChildType = typename TraitsType::ChildType;
	/// Number of tracking lists.
	static constexpr TupleIdx t_num_lists = TraitsType::t_num_lists;
	/// Dimension of grid.
	static constexpr UINT t_dims = TraitsType::t_dims;
	/// D-dimensional integer vector.
	using VecDi = Felt::VecDi<t_dims>;

protected:
	/// Grid of partitions with tracking list(s) of active partitions.
	using ChildrenGrid = Impl::Tracked::MultiByRef<ChildType, t_dims, t_num_lists>;
private:
	/// Size of a child sub-grid.
	VecDi m_child_size;
	/// Grid of child grids.
	ChildrenGrid m_children;
	/// Mutex used to synchrnonise the adding/removing of elements from the tracking list(s).
	std::mutex	m_mutex;

protected:
	/// Deleted default constructor.
	Children() = delete;

	/**
	 * Construct and initialise children grid to hold child sub-grids.
	 */
	Children(
		const VecDi& size_, const VecDi& offset_, const VecDi& child_size_,
		const ChildType& background_
	) :
		m_child_size(child_size_),
		m_children(
			calc_children_size(size_, child_size_),
			(offset_.array() / child_size_.array()).matrix(),
			background_
		)
	{
		// Set each child sub-grid's size and offset.
		for (UINT idx = 0; idx < m_children.data().size(); idx++)
		{
			// Position of child in children grid.
			const VecDi& pos_child = m_children.index(idx);
			// Position of child in children grid, without offset.
			const VecDi& pos_child_offset = pos_child - m_children.offset();
			// Scaled position of child == position in world space, without offset.
			const VecDi& offset_child_offset = (
				pos_child_offset.array() * m_child_size.array()
			).matrix();
			// Position of child in world space, including offset.
			const VecDi& offset_child = offset_child_offset + offset_;

			m_children.data()[idx].resize(m_child_size, offset_child);
		}
	}

	/**
	 * Get children grid - the spatial partition grid that stores the child sub-grids.
	 *
	 * @return multi-index multi-list tracked grid storing/tracking Child objects.
	 */
	ChildrenGrid& children()
	{
		return m_children;
	}

	/**
	 * @copydoc Children::children()
	 */
	const ChildrenGrid& children() const
	{
		return m_children;
	}

	/**
	 * Get size of child sub-grids.
	 *
	 * @return size of child sub-grid.
	 */
	const VecDi& child_size() const
	{
		return m_child_size;
	}

	/**
	 * Add a spatial partition to children grid's tracking sub-grid.
	 *
	 * Uses mutex for thread safety. Activates the child grid.
	 *
	 * @param pos_idx_child_ position index of child sub-grid in parent grid.
	 * @param list_idx_ index of tracking list used to track position.
	 */
	void track_child(const PosIdx pos_idx_child_, const TupleIdx list_idx_)
	{
		FELT_DEBUG(m_children.assert_pos_idx_bounds(pos_idx_child_, "track:"));
		if (m_children.lookup().is_tracked(pos_idx_child_, list_idx_))
			return;
		std::lock_guard<std::mutex> lock(m_mutex);
		if (m_children.lookup().is_tracked(pos_idx_child_, list_idx_))
			return;
		ChildType& child = m_children.get(pos_idx_child_);
		if (!child.is_active())
			child.activate();
		m_children.lookup().track(pos_idx_child_, list_idx_);
	}

	/**
	 * Bulk add children to tracking list, activating if not already active.
	 *
	 * Not thread-safe.
	 *
	 * @param grid_mask_ grid to match partition activation with.
	 */
	template <class MaskGrid>
	void track_children(const MaskGrid& grid_mask_)
	{
		for (TupleIdx list_idx = 0; list_idx < grid_mask_.children().t_num_lists; list_idx++)
		{
			for (const PosIdx pos_idx_child : grid_mask_.children().lookup().list(list_idx))
			{
				if (m_children.lookup().is_tracked(pos_idx_child, list_idx))
					continue;
				ChildType& child = m_children.get(pos_idx_child);
				if (!child.is_active())
					child.activate();
				m_children.lookup().track(pos_idx_child, list_idx);
			}
		}
	}

	template <class MaskGrid>
	void reset(const MaskGrid& grid_mask_)
	{
		for (TupleIdx list_idx = 0; list_idx < t_num_lists; list_idx++)
		{
			for (const PosIdx pos_idx_child : m_children.lookup().list(list_idx))
			{
				ChildType& child = m_children.get(pos_idx_child);
				m_children.lookup().untrack(pos_idx_child, list_idx);

				// If the master grid is not tracking this child, then untrack it from tracking
				// under this list id, potentially destroying it.
				if (
					!grid_mask_.children().lookup().is_tracked(pos_idx_child) &&
					!m_children.lookup().is_tracked(pos_idx_child)
				) {
					child.deactivate();
				}

				// If the child has not been destroyed by the above, then reset as normal (loop
				// over tracking list resetting values in grid, before resizing list to 0).
				if (child.is_active())
				{
					child.reset();
				}
				// If the child was destroyed above, then no need to loop over grid resetting
				// values, so just reset list.
				else
				{
					child.list(list_idx).clear();
				}
			}
		}
	}

	/**
	 * Calculate the position of a child grid (i.e. partition) given the position of leaf grid node.
	 *
	 * @param pos_leaf_
	 * @return position of spatial partition in which leaf position lies.
	 */
	PosIdx pos_idx_child (const VecDi& pos_leaf_) const
	{
		// Position of leaf, without offset.
		const VecDi& pos_leaf_offset =  pos_leaf_ - pself->offset();
		// Position of child grid containing leaf, without offset.
		const VecDi& pos_child_offset = (pos_leaf_offset.array() / m_child_size.array()).matrix();
		// Position of child grid containing leaf, including offset.
		const VecDi& pos_child = pos_child_offset + m_children.offset();
		// Encode child position as an index.
		return m_children.index(pos_child);
	}


	template<typename Fn>
	void leafs(const TupleIdx layer_idx_, Fn fn_) const
	{
		const PosArray& pos_idxs_child = m_children.lookup().list(layer_idx_);
		const ListIdx num_childs = pos_idxs_child.size();

		FELT_PARALLEL_FOR(num_childs)
		for (ListIdx list_idx = 0; list_idx < num_childs; list_idx++)
		{
			const PosIdx pos_idx_child = pos_idxs_child[list_idx];
			const ChildType& child = m_children.get(pos_idx_child);
			for (const PosIdx pos_idx_leaf : child.lookup().list(layer_idx_))
			{
				fn_(pos_idx_child, pos_idx_leaf);
			}
		}
	}


	std::mutex& mutex_children()
	{
		return m_mutex;
	}

private:
	/**
	 * Calculate required size of children grid to contain child sub-grids.
	 */
	VecDi calc_children_size(const VecDi& size_, const VecDi& child_size_)
	{
		VecDi children_size = (size_.array() / child_size_.array()).matrix();

		if ((children_size.array() * child_size_.array()).matrix() != size_)
		{
			children_size += VecDi::Constant(1);
		}

		return children_size;
	}
};


template <class Derived>
class Lookup
{
private:
	/// Traits of derived class.
	using TraitsType = Traits<Derived>;
	/// Child grid type.
	using ChildType = typename TraitsType::ChildType;
	/// Dimension of grid.
	static constexpr UINT t_dims = TraitsType::t_dims;
	/// D-dimensional integer vector.
	using VecDi = Felt::VecDi<t_dims>;

protected:
	/**
	 * Add a leaf position to be tracked to given tracking list.
	 *
	 * Descend to relevant child grid to track to their tracking structure.
	 *
	 * @param pos_idx_leaf_ position to track.
	 * @param list_idx_ tracking list id.
	 */
	void track(const VecDi& pos_leaf_, const TupleIdx list_idx_)
	{
		const PosIdx pos_idx_child_ = pself->pos_idx_child(pos_leaf_);
		pself->track_child(pos_idx_child_, list_idx_);
		ChildType& child = pself->children().get(pos_idx_child_);
		const PosIdx pos_idx_leaf = child.index(pos_leaf_);
		child.track(pos_idx_leaf, list_idx_);
	}

	/**
	 * Add a leaf position to be tracked to given tracking list.
	 *
	 * Descend to relevant child grid to track to their tracking structure.
	 *
	 * @param pos_idx_leaf_ position to track.
	 * @param pos_idx_child_ position of child sub-grid containing leaf position to track.
	 * @param list_idx_ tracking list id.
	 */
	void track(const PosIdx pos_idx_child_, const PosIdx pos_idx_leaf_, const TupleIdx list_idx_)
	{
		pself->track_child(pos_idx_child_, list_idx_);
		FELT_DEBUG(
			pself->children().get(pos_idx_child_).assert_pos_idx_bounds(pos_idx_leaf_, "track:")
		);
		pself->children().get(pos_idx_child_).track(pos_idx_leaf_, list_idx_);
	}
};


template <class Derived>
class Tracked
{
private:
	/// Traits of derived class.
	using TraitsType = Traits<Derived>;
	/// Child grid type.
	using ChildType = typename TraitsType::ChildType;
	/// Leaf type.
	using LeafType = typename TraitsType::LeafType;
	/// Dimension of grid.
	static constexpr UINT t_dims = TraitsType::t_dims;
	/// D-dimensional integer vector.
	using VecDi = Felt::VecDi<t_dims>;

protected:
	/**
	 * Add a leaf position to be tracked to given tracking list.
	 *
	 * Descend to relevant child grid to track to their tracking structure.
	 *
	 * @param pos_idx_leaf_ position to track.
	 * @param list_idx_ tracking list id.
	 */
	void track(const LeafType val_, const VecDi& pos_leaf_, const TupleIdx list_idx_)
	{
		const PosIdx pos_idx_child_ = pself->pos_idx_child(pos_leaf_);
		pself->track_child(pos_idx_child_, list_idx_);
		ChildType& child = pself->children().get(pos_idx_child_);
		const PosIdx pos_idx_leaf = child.index(pos_leaf_);
		child.track(val_, pos_idx_leaf, list_idx_);
	}

	/**
	 * Add a leaf position to be tracked to given tracking list.
	 *
	 * Descend to relevant child grid to track to their tracking structure.
	 *
	 * @param pos_idx_leaf_ position to track.
	 * @param pos_idx_child_ position of child sub-grid containing leaf position to track.
	 * @param list_idx_ tracking list id.
	 */
	void track(
		const LeafType val_, const PosIdx pos_idx_child_, const PosIdx pos_idx_leaf_,
		const TupleIdx list_idx_
	) {
		pself->track_child(pos_idx_child_, list_idx_);
		FELT_DEBUG(
			pself->children().get(pos_idx_child_).assert_pos_idx_bounds(pos_idx_leaf_, "track:")
		);
		ChildType& child = pself->children().get(pos_idx_child_);
		child.track(val_, pos_idx_leaf_, list_idx_);
	}
};


template <class Derived>
class Untrack
{
private:
	/// Traits of derived class.
	using TraitsType = Traits<Derived>;
	/// Child grid type.
	using ChildType = typename TraitsType::ChildType;
	/// Leaf type.
	using LeafType = typename TraitsType::LeafType;
	/// Dimension of grid.
	static constexpr UINT t_dims = TraitsType::t_dims;
	/// D-dimensional integer vector.
	using VecDi = Felt::VecDi<t_dims>;

protected:
	void untrack(
		const LeafType background_, const PosIdx pos_idx_child_, const PosIdx pos_idx_leaf_,
		const TupleIdx list_idx_
	) {
		ChildType& child = pself->children().get(pos_idx_child_);

		// Untrack position in child sub-grid.
		child.lookup().untrack(pos_idx_leaf_, list_idx_);
		child.set(pos_idx_leaf_, background_);
		const PosArray& pos_idxs = child.lookup().list(list_idx_);

		// If tracking list is empty in child, untrack parent.  No race condition here, as long as
		// we stick to one thread per child rule.
		if (pos_idxs.size() == 0)
		{
			// Scope for lock.
			{
				// Mutex lock for children grid.
				std::lock_guard<std::mutex> lock(pself->mutex_children());
				// Untrack this list in children grid.
				pself->children().lookup().untrack(pos_idx_child_, list_idx_);
			}

			// If no position is being tracked at all in any tracking list, then deactivate
			// child. Otherwise just reset to background value.
			if (not pself->children().lookup().is_tracked(pos_idx_child_))
				child.deactivate(background_);
		}
	}

	void retrack(
		const PosIdx pos_idx_child_, const PosIdx pos_idx_leaf_, const TupleIdx list_idx_from_,
		const TupleIdx list_idx_to_
	) {
		ChildType& child = pself->children().get(pos_idx_child_);

		#ifdef FELT_DEBUG_ENABLED
		if (not pself->children().lookup().is_tracked(pos_idx_child_))
		{
			std::stringstream strs;
			strs << "Attempting to move lists within an inactive child: " <<
				format(child.index(pos_idx_leaf_)) << " from list " << list_idx_from_ <<
				" to list " << list_idx_to_ << " in partition " <<
				format(pself->children().index(pos_idx_child_));
			std::string str = strs.str();
			throw std::domain_error(str);
		}
		#endif

		// Move position between tracking lists in child grid.
		child.lookup().untrack(pos_idx_leaf_, list_idx_from_);
		child.lookup().track(pos_idx_leaf_, list_idx_to_);

		// If child is not tracked by target list or child should be untracked by source list,
		// then apply mutex and track/untrack as necessary. No race condition, provided we stick to
		// one child per thread rule.
		if (
			not pself->children().lookup().is_tracked(pos_idx_child_, list_idx_to_) ||
			child.lookup().list(list_idx_from_).size() == 0
		) {
			// Mutex lock for children grid modifications.
			std::lock_guard<std::mutex> lock(pself->mutex_children());
			FELT_DEBUG(child.assert_pos_idx_bounds(pos_idx_leaf_, "retrack"));
			// Ensure parent grid tracks this child in target list.
			pself->children().lookup().track(pos_idx_child_, list_idx_to_);
			// If child's source list is now empty, stop tracking in parent grid.
			if (child.lookup().list(list_idx_from_).size() == 0)
				pself->children().lookup().untrack(pos_idx_child_, list_idx_from_);
		}
	}
};

template <class Derived>
class Accessor
{
private:
	/// Traits of derived class.
	using TraitsType = Traits<Derived>;
	/// Child grid type.
	using ChildType = typename TraitsType::ChildType;
	/// Leaf type.
	using LeafType = typename TraitsType::LeafType;
	/// Dimension of grid.
	static constexpr UINT t_dims = TraitsType::t_dims;
	/// D-dimensional integer vector.
	using VecDi = Felt::VecDi<t_dims>;

protected:
	/**
	 * Get the leaf grid node at pos by navigating to the correct partition.
	 *
	 * @param pos_ position in grid to fetch.
	 * @return value stored at grid point.
	 */
	LeafType get (const VecDi& pos_) const
	{
		const PosIdx pos_idx_child = pself->pos_idx_child(pos_);
		const ChildType& child = pself->children().get(pos_idx_child);
		return child.get(pos_);
	}

	/**
	 * Set the leaf grid node at pos by navigating to the correct partition.
	 *
	 * @param pos_ position in grid to store at.
	 * @param value_ value to store at grid point.
	 */
	void set (const VecDi& pos_, const LeafType value_)
	{
		const PosIdx pos_idx_child = pself->pos_idx_child(pos_);
		ChildType& child = pself->children().get(pos_idx_child);
		child.set(pos_, value_);
	}
};


template <class Derived>
class Snapshot
{
private:
	/// Traits of derived class.
	using TraitsType = Traits<Derived>;
	/// Child grid type.
	using ChildType = typename TraitsType::ChildType;
	/// Leaf type.
	using LeafType = typename TraitsType::LeafType;
	/// Dimension of grid.
	static constexpr UINT t_dims = TraitsType::t_dims;
	/// D-dimensional integer vector.
	using VecDi = Felt::VecDi<t_dims>;
	/// Simple non-partitioned grid type for snapshotting the partitioned data, e.g. for tests.
	using SnapshotGrid = Impl::Grid::Snapshot<LeafType, t_dims>;
protected:
	/**
	 * Non-partitioned Grid class to store snapshots of partitioned grids.
	 *
	 * Useful for serialisation or logging.
	 */
	using SnapshotGridPtr = std::unique_ptr<SnapshotGrid>;

public:
	SnapshotGridPtr snapshot() const
	{
		SnapshotGridPtr psnapshot = std::make_unique<SnapshotGrid>(
			pself->size(), pself->offset(), LeafType()
		);

		const VecDi pos_max = (
			psnapshot->size() - VecDi::Constant(1) + psnapshot->offset()
		);
		const PosIdx pos_idx_max = psnapshot->index(pos_max);

		for (PosIdx pos_idx = 0; pos_idx <= pos_idx_max; pos_idx++)
		{
			const VecDi pos = psnapshot->index(pos_idx);
			psnapshot->set(pos_idx, pself->get(pos));
		}

		return psnapshot;
	}

	void snapshot(const SnapshotGridPtr& psnapshot)
	{
		for (PosIdx pos_idx = 0; pos_idx < psnapshot->data().size(); pos_idx++)
		{
			const LeafType val = psnapshot->get(pos_idx);
			const VecDi& pos = psnapshot->index(pos_idx);

			const PosIdx pos_idx_child = pself->pos_idx_child(pos);
			ChildType& child = pself->children().get(pos_idx_child);
			const PosIdx pos_idx_leaf = child.index(pos);

			if (!child.is_active())
			{
				if (val == child.background())
					continue;
				child.activate();
			}
			child.set(pos_idx_leaf, psnapshot->get(pos_idx));
		}
	}

	void operator=(std::initializer_list<LeafType> vals_)
	{
		SnapshotGridPtr snap = snapshot();
		snap->data() = vals_;
		snapshot(snap);
	}
};

} // Partitioned.
} // Mixin.
} // Impl.
} // Felt.


#endif /* INCLUDE_FELT_IMPL_MIXIN_PARTITIONEDMIXIN_HPP_ */
