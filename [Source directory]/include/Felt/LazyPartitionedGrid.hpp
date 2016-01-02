/*
 * LazyPartitionedGrid.hpp
 *
 *  Created on: 2 Jan 2016
 *      Author: dave
 */

#ifndef INCLUDE_FELT_LAZYPARTITIONEDGRID_HPP_
#define INCLUDE_FELT_LAZYPARTITIONEDGRID_HPP_

#include "PartitionedGrid.hpp"

namespace felt
{

/**
 * @ingroup Classes
 * Lazy spatially partitioned grid storing arbitrary data.
 *
 * Create/destroy child grids as necessary.
 *
 * @tparam T the type of data to store in the grid.
 * @tparam D the dimensions of the grid.
 */
template <typename T, UINT D>
class LazyPartitionedGrid : public PartitionedGridBase<LazyPartitionedGrid<T, D> >
{
public:
	using LeafType = T;

protected:
	LeafType	m_default_val;

public:
	/// Base class.
	using Base = PartitionedGridBase<LazyPartitionedGrid<T, D> >;
	/**
	 * Constructor
	 *
	 * @param size_ spatial size of the whole grid.
	 * @param offset_ spatial offset of the whole grid.
	 * @param size_partition_ size of each spatial partition.
	 * @param default_val_ leaf grid node default value when child grid is not initialised.
	 * @param delta_ grid delta value used for spatial derivatives (dx).
	 */
	LazyPartitionedGrid (
		const VecDu& size_, const VecDi& offset_,
		const VecDu& size_partition_ = VecDu::Constant(DEFAULT_PARTITION),
		const LeafType& default_val_, const FLOAT& delta_ = 1
	) : Base()
	{
		this->init(size_, offset_, size_partition_, default_val_, delta_);
	}

	/**
	 * Initialisation function called by non-trivial constructor and
	 * subclasses.
	 *
	 * @param size_ the dimensions of the whole grid.
	 * @param offset_ the offset of the whole grid.
	 * @param size_partition_ the dimensions of each spatial partition.
	 * @param default_val_ leaf grid node default value when child grid is not initialised.
	 * @param delta_ the grid delta value used for spatial derivatives (dx).
	 */
	void init (
		const VecDu& size_, const VecDi& offset,
		const VecDu& size_partition_ = VecDu::Constant(DEFAULT_PARTITION),
		const LeafType& default_val_, const FLOAT& delta_ = 1.0f
	) {
		m_default_val = default_val_;
		Base::init(size_, offset, delta_, size_partition_);
	}

	/**
	 * Set grid size, but do not initialise child grids.
	 *
	 * @param size_ size of the overall grid.
	 */
	void size (const VecDu& size_)
	{
		Base::size(size_);
		this->m_size = size_;
	}

	/**
	 * Add a spatial partition to branch grid's tracking subgrid.
	 *
	 * Uses mutex for thread safety.
	 *
	 * @param pos_ position in spatial partition grid to track.
	 * @param arr_idx_ index of tracking list used to track position.
	 * @return true if position was added to tracking grid, false if it was already added.
	 */
	bool add_child(const VecDi& pos_, const UINT& arr_idx_ = 0)
	{
		if (m_grid_branch.is_active(pos_, arr_idx_))
			return false;
		std::lock_guard<std::mutex> lock(m_mutex_update_branch);
		this->branch().get(pos_).size(this->m_usize_child);
		return this->branch().add(pos_, arr_idx_);
	}
};

/**
 * @ingroup Traits
 * Traits for GridBase to understand LazyPartitionedGrid.
 *
 * @tparam T the type of data to store in the grid.
 * @tparam D the number of dimensions of the grid.
 */
template <typename T, UINT D>
struct GridBaseTraits<LazyPartitionedGrid<T, D> > : GridBaseTraits<PartitionedGrid<T, D> >
{};


/**
 * @ingroup Traits
 * Traits for PartitionedGridBase to understand PartitionedGrid.
 *
 * @tparam T the type of data to store in the grid.
 * @tparam D the number of dimensions of the grid.
 */
template <typename T, UINT D>
struct PartitionedGridBaseTraits<LazyPartitionedGrid<T, D> >
	: PartitionedGridBaseTraits<PartitionedGrid<T, D> >
{};


} // End namespace felt.

#endif /* INCLUDE_FELT_LAZYPARTITIONEDGRID_HPP_ */
