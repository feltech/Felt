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
	using RetType = LeafType;
	/// Base class.
	using Base = PartitionedGridBase<LazyPartitionedGrid<T, D> >;
	using BaseBase = typename Base::Base;
	using MixinType = typename Base::Child;
	using Child = typename Base::Child;
	using VecDu = typename Base::VecDu;
	using VecDi = typename Base::VecDi;

protected:
	LeafType	m_default_val;

public:
	/**
	 * Constructor
	 *
	 * @snippet test_LazyPartitionedGrid.cpp Lazy grid basic usage
	 *
	 * @param size_ spatial size of the whole grid.
	 * @param offset_ spatial offset of the whole grid.
	 * @param size_partition_ size of each spatial partition.
	 * @param default_val_ leaf grid node default value when child grid is not initialised.
	 * @param delta_ grid delta value used for spatial derivatives (dx).
	 */
	LazyPartitionedGrid (
		const VecDu& size_, const VecDi& offset_, const LeafType& default_val_,
		const VecDu& size_partition_ = VecDu::Constant(DEFAULT_PARTITION),
		const FLOAT& delta_ = 1
	) : Base()
	{
		this->init(size_, offset_, default_val_, size_partition_, delta_);
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
		const VecDu& size_, const VecDi& offset, const LeafType& default_val_,
		const VecDu& size_partition_ = VecDu::Constant(DEFAULT_PARTITION),
		const FLOAT& delta_ = 1.0f
	) {
		m_default_val = default_val_;
		Base::init(size_, offset, size_partition_, delta_);
	}

	/**
	 * Add a spatial partition to branch grid's tracking subgrid.
	 *
	 * Will initialise child grid if necessary. Uses mutex for thread safety.
	 *
	 * @snippet test_LazyPartitionedGrid.cpp Lazy grid basic usage
	 *
	 * @param pos_ position in spatial partition grid to track.
	 * @param arr_idx_ index of tracking list used to track position.
	 * @return true if position was added to tracking grid, false if it was already added.
	 */
	bool add_child(const VecDi& pos_, const UINT& arr_idx_ = 0)
	{
		if (this->m_grid_branch.is_active(pos_, arr_idx_))
			return false;
		Child& child = this->branch().get(pos_);
		std::lock_guard<std::mutex> lock(this->m_mutex_update_branch);
		if (child.data().size() == 0)
			child.size(this->m_usize_child);
		return this->branch().add(pos_, arr_idx_);
	}


	/**
	 * Remove a spatial partition from branch grid's tracking subgrid.
	 *
	 * Destroys child grid. Uses mutex for thread safety.
	 *
	 * @snippet test_LazyPartitionedGrid.cpp Lazy grid basic usage
	 *
	 * @param pos_ position of spatial partition to stop tracking.
	 * @param arr_idx_ index of tracking list used to track position.
	 */
	void remove_child(const VecDi& pos_, const UINT& arr_idx_ = 0)
	{
		if (!this->m_grid_branch.is_active(pos_, arr_idx_))
			return;
		Child& child = this->branch().get(pos_);
		std::lock_guard<std::mutex> lock(this->m_mutex_update_branch);
		if (child.data().size() > 0)
			child.size(VecDu(0,0,0));
		this->branch().remove(pos_, arr_idx_);
	}

	/**
	 * Get the leaf grid node at pos by navigating to the correct partition.
	 *
	 * Return default value if child grid not yet instantiated.
	 *
	 * @snippet test_LazyPartitionedGrid.cpp Lazy grid basic usage
	 *
	 * @param pos_ position in grid to fetch.
	 * @return value stored at grid point.
	 */
	RetType& get (const VecDi& pos_)
	{
		const VecDi& pos_child = this->pos_child(pos_);
		Child& child = this->branch().get(pos_child);
		if (child.data().size() == 0)
			return m_default_val;
		return child.get(pos_);
	}

	/**
	 * @copydoc LazyPartitionedGrid::get(const VecDi&)
	 */
	const RetType& get (const VecDi& pos_) const
	{
		const Child& child = this->branch()(pos_child(pos_));
		if (child.data().size() == 0) {
			return m_default_val;
		}
		return child.get(pos_);
	}

	/**
	 * Set grid size, but do not initialise child grids.
	 *
	 * @param size_ size of the overall grid.
	 */
	void size (const VecDu& size_)
	{
		BaseBase::size(size_);
		this->m_size = size_;
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
{
	/// The class inheriting from the base.
	using ThisType = LazyPartitionedGrid<T, D>;
};


/**
 * @ingroup Traits
 * Traits for PartitionedGridBase to understand LazyPartitionedGrid.
 *
 * @tparam T the type of data to store in the grid.
 * @tparam D the number of dimensions of the grid.
 */
template <typename T, UINT D>
struct PartitionedGridBaseTraits<LazyPartitionedGrid<T, D> >
	: PartitionedGridBaseTraits<PartitionedGrid<T, D> >
{
	/// The class inheriting from the base.
	using ThisType = LazyPartitionedGrid<T, D>;
	/// Mixin type whose signature to spoof - in this case a standard GridBase.
	using MixinType = GridBase<ThisType>;
};


} // End namespace felt.

#endif /* INCLUDE_FELT_LAZYPARTITIONEDGRID_HPP_ */
