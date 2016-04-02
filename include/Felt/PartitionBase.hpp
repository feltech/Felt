#ifndef INCLUDE_FELT_PARTITIONBASE_HPP_
#define INCLUDE_FELT_PARTITIONBASE_HPP_

#include <mutex>
#include "TrackedGrid.hpp"

namespace felt
{
/**
 * Base class for spatially partitioned structures.
 *
 * A TrackedGrid is used to store and track arbitrary Child structures. The grid as a whole has a
 * size, which is the size of a Child multiplied by the size of the TrackedGrid.  The Child might
 * not be another grid type (e.g. see PartitionedArray).
 *
 * @tparam Derived the CRTP subclass type.
 */
template <class Derived>
class PartitionBase
{
public:
	/// CRTP subclass type.
	using DerivedType = Derived;
	using Traits = GridTraits<Derived>;

	/// Child object to store in each partition.
	using Child = typename Traits::ChildType;

	/// Number of tracking lists of points.
	static const UINT NUM_LISTS = Traits::NumLists;

	/// Grid of partitions with tracking list(s) of active partitions.
	using ChildrenGrid = TrackedGrid<Child, Traits::Dims, NUM_LISTS>;
	/// D-dimensional unsigned int vector.
	using VecDu = typename ChildrenGrid::VecDu;
	/// D-dimensional signed int vector.
	using VecDi = typename ChildrenGrid::VecDi;


protected:
	/// Grid of partitions with tracking list(s) of active grid points.
	ChildrenGrid	m_grid_children;

	/// Mutex used to synchrnonise the the adding/removing of elements from the tracking list(s).
	std::mutex	m_mutex_update_branch;

	/// A convenience vector of the (unsigned) size of a partition.
	VecDu	m_usize_child;
	/// A convenience vector of the (signed) size of a partition.
	VecDi 	m_isize_child;

public:
	/**
	 * Empty destructor.
	 */
	~PartitionBase ()
	{}


	/**
	 * Explicitly defined default constructor.
	 */
	PartitionBase () = default;

	/**
	 * Construct a spatially partitioned data structure with given size, spatial offset and
	 * partition size.
	 *
	 * @param size_ spatial size of grid (child size * parent size).
	 * @param offset_ spatial offset of the grid.
	 * @param partition_size_ spatial size of the child data structure.
	 */
	PartitionBase (const VecDu& size_, const VecDi& offset_, const VecDu& partition_size_)
		: m_grid_children()
	{
		this->init(size_, offset_, partition_size_);
	}


	/**
	 * Initialisation method to be called by non-trivial constructor or
	 * subclasses.
	 *
	 * Similar to a Grid::init(), with the addition of setting the size of
	 * spatial partitions.
	 *
	 * @param size_ spatial size of grid (child size * parent size).
	 * @param offset_ spatial offset of the grid.
	 * @param partition_size_ spatial size of the child data structure.

	 */
	void init (const VecDu& size_, const VecDi& offset_, const VecDu& partition_size_) {
		self->init(partition_size_);
		self->size(size_);
		self->offset(offset_);
	}

	/**
	 * Initialise size of spatial partitions.
	 *
	 * @param partition_size_ spatial size of the child data structure.
	 */
	void init (const VecDu& partition_size_)
	{
		m_usize_child = partition_size_;
		m_isize_child = partition_size_.template cast<INT>();
	}

	/**
	 * Get size of a spatial partition.
	 */
	const VecDu& child_size() const
	{
		return m_usize_child;
	}

	/**
	 * Get TrackedGrid children grid - the spatial partition grid that stores the PartitionBase::Child
	 * objects.
	 *
	 * @return TrackedGrid storing/tracking Child objects.
	 */
	ChildrenGrid& children ()
	{
		return m_grid_children;
	}

	/**
	 * Get TrackedGrid children grid - the spatial partition grid that stores the PartitionBase::Child
	 * objects.
	 *
	 * @return TrackedGrid storing/tracking Child objects.
	 */
	const ChildrenGrid& children () const
	{
		return m_grid_children;
	}

	/**
	 * Reshape grid, computing the size of the children grid.
	 *
	 * The children grid will be increased in size by one, if required, to ensure all leaf nodes are
	 * completely contained.
	 *
	 * @param grid_size_ overall size of grid.
	 */
	void size (const VecDu& grid_size_)
	{
		VecDu branch_size = (
			grid_size_.array() / m_usize_child.array()
		).matrix();

		if (
			(branch_size.array() * m_usize_child.array()).matrix()
			!= grid_size_
		) {
			branch_size += VecDu::Constant(1);
		}

		m_grid_children.size(branch_size);
	}

	/**
	 * Calculate the offset of the children grid from the given offset and the size of a spatial
	 * partition (the child size).
	 *
	 * @param grid_offset_ overall spatial offset the grid
	 */
	void offset (const VecDi& grid_offset_)
	{
		const VecDi& branch_offset = (
			grid_offset_.array() / m_isize_child.array()
		).matrix();

		m_grid_children.offset(branch_offset);
	}

	/**
	 * Add a spatial partition to children grid's tracking subgrid.
	 *
	 * Uses mutex for thread safety.
	 *
	 * @param pos_ position in spatial partition grid to track.
	 * @param arr_idx_ index of tracking list used to track position.
	 * @return true if position was added to tracking grid, false if it was already added.
	 */
	bool add_child(const VecDi& pos_, const UINT arr_idx_ = 0)
	{
		if (m_grid_children.is_active(pos_, arr_idx_))
			return false;
		std::lock_guard<std::mutex> lock(m_mutex_update_branch);
		return this->children().add(pos_, arr_idx_);
	}

	/**
	 * Remove a spatial partition from children grid's tracking subgrid.
	 *
	 * Uses mutex for thread safety.
	 *
	 * @param pos_ position of spatial partition to stop tracking.
	 * @param arr_idx_ index of tracking list used to track position.
	 */
	void remove_child(const VecDi& pos_, const UINT arr_idx_ = 0)
	{
		if (!m_grid_children.is_active(pos_, arr_idx_))
			return;
		std::lock_guard<std::mutex> lock(m_mutex_update_branch);
		this->children().remove(pos_, arr_idx_);
	}


	/**
	 *
	 * @param pos_
	 * @return
	 */
	bool is_child_active(const VecDi& pos_child_) const
	{
		using NULL_IDX_TYPE = typename ChildrenGrid::Lookup::Traits::NULL_IDX_TYPE;

		const NULL_IDX_TYPE idxs = this->children().lookup().get(pos_child_);

		return idxs != ChildrenGrid::Lookup::Traits::NULL_IDX_DATA;
	}

	/**
	 * Reset the tracking list at given index in the children grid.
	 *
	 * Removes all spatial partitions from the tracking subgrid for given list index.
	 *
	 * @param arr_idx_ index of tracking list to reset.
	 */
	void reset(const UINT arr_idx_ = 0)
	{
		this->children().reset(arr_idx_);
	}
};

} // End namespace felt.


#endif /* INCLUDE_FELT_PARTITIONBASE_HPP_ */
