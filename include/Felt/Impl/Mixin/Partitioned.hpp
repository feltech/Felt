#ifndef INCLUDE_FELT_IMPL_MIXIN_PARTITIONED_HPP_
#define INCLUDE_FELT_IMPL_MIXIN_PARTITIONED_HPP_
#include <mutex>

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
	using Child = typename TraitsType::ChildType;
	/// Number of tracking lists.
	static constexpr UINT NumLists = TraitsType::NumLists;
	/// Dimension of grid.
	static constexpr UINT Dims = TraitsType::Dims;
	/// D-dimensional integer vector.
	using VecDi = Felt::VecDi<Dims>;

protected:
	/// Grid of partitions with tracking list(s) of active partitions.
	using ChildrenGrid = Impl::Tracked::MultiByRef<Child, Dims, NumLists>;

protected:
	/// Grid of child grids.
	ChildrenGrid m_children;

private:
	/// Mutex used to synchrnonise the adding/removing of elements from the tracking list(s).
	std::mutex	m_mutex;
	/// Size of a child sub-grid.
	VecDi m_child_size;

protected:
	/// Deleted default constructor.
	Children() = delete;

	/**
	 * Construct and initialise children grid to hold child sub-grids.
	 */
	Children(
		const VecDi& size_, const VecDi& offset_, const VecDi& child_size_,
		const Child& background_
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
	 * Add a leaf position to be tracked to given tracking list.
	 *
	 * Descend to relevant child grid to add to their tracking structure.
	 *
	 * @param pos_idx_leaf_ position to track.
	 * @param list_idx_ tracking list id.
	 */
	void add(const VecDi& pos_leaf_, const UINT list_idx_)
	{
		const Idx pos_idx_child_ = pos_idx_child(pos_leaf_);
		add_child(pos_idx_child_, list_idx_);
		const Idx pos_idx_leaf = m_children.get(pos_idx_child_).index(pos_leaf_);
		add(pos_idx_child_, pos_idx_leaf, list_idx_);
	}


	/**
	 * Add a leaf position to be tracked to given tracking list.
	 *
	 * Descend to relevant child grid to add to their tracking structure.
	 *
	 * @param pos_idx_leaf_ position to track.
	 * @param pos_idx_child_ position of child sub-grid containing leaf position to track.
	 * @param list_idx_ tracking list id.
	 */
	void add(const Idx pos_idx_child_, const Idx pos_idx_leaf_, const UINT list_idx_)
	{
		FELT_DEBUG(m_children.get(pos_idx_child_).assert_pos_idx_bounds(pos_idx_leaf_, "add:"));
		m_children.get(pos_idx_child_).add(pos_idx_leaf_, list_idx_);
	}

	/**
	 * Add a spatial partition to children grid's tracking sub-grid.
	 *
	 * Uses mutex for thread safety. Activates the child grid.
	 *
	 * @param pos_idx_child_ position index of child sub-grid in parent grid.
	 * @param list_idx_ index of tracking list used to track position.
	 */
	void add_child(const Idx pos_idx_child_, const UINT list_idx_)
	{
		FELT_DEBUG(m_children.assert_pos_idx_bounds(pos_idx_child_, "add:"));
		std::lock_guard<std::mutex> lock(m_mutex);
		if (m_children.lookup().is_active(pos_idx_child_, list_idx_))
			return;
		Child& child = m_children.get(pos_idx_child_);
		if (!child.data().size())
			child.activate();
		m_children.lookup().add(pos_idx_child_, list_idx_);
	}

private:
	/**
	 * Calculate the position of a child grid (i.e. partition) given the position of leaf grid node.
	 *
	 * @param pos_leaf_
	 * @return position of spatial partition in which leaf position lies.
	 */
	Idx pos_idx_child (const VecDi& pos_leaf_) const
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


} // Partitioned.
} // Mixin.
} // Impl.
} // Felt.


#endif /* INCLUDE_FELT_IMPL_MIXIN_PARTITIONED_HPP_ */
