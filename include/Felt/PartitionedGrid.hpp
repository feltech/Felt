#ifndef PARTITIONEDGRID_HPP_
#define PARTITIONEDGRID_HPP_

#include <memory>
#include "PartitionBase.hpp"

namespace felt
{
/**
 * Base class for spatially partitioned grid storing arbitrary values.
 *
 * Uses multiple inheritance from child grid (MixinType) class, spoofing the signature,
 * and PartitionBase for actual storage in a (single-level) tree hierarchy,
 * i.e. PartitionedGridBase -> Child -> LeafType.
 */
template <class Derived>
class PartitionedGridBase
	: public GridTraits<Derived>::MixinType,
	  public PartitionBase<PartitionedGridBase<Derived> >
{
public:
	/// CRTP derived class.
	using DerivedType = Derived;
	using ThisType = PartitionedGridBase<DerivedType>;
	/// PartitionBase base class.
	using Base = PartitionBase<ThisType>;
	/// Traits of CRTP derived class
	using Traits = GridTraits<DerivedType>;
	/// Base grid class to partition.
	using MixinType = typename Traits::MixinType;
	/// Child grid class.  Each spatial partition contains one Child grid.
	using Child = typename Traits::ChildType;
	/// Type of data stored in leaf grid nodes.
	using LeafType = typename Traits::LeafType;
	/// Dimension of overall grid.
	static const UINT Dims = Traits::Dims;
	/**
	 * MultiTrackedGrid type storing spatial partitions.
	 *
	 * Each grid node stores a Child grid and tracks active nodes using one or more tracking lists.
	 */
	using ChildrenGrid = typename Base::ChildrenGrid;

	using VecDu = typename MixinType::VecDu;
	using VecDi = typename MixinType::VecDi;

	/**
	 * Non-partitioned Grid class to store snapshots of partitioned grids.
	 *
	 * Useful for serialisation or logging.
	 */
	using SnapshotGrid = Grid<LeafType, Dims>;

protected:
	/// Non-partitioned snapshot to maintain for serialisation or logging.
	std::unique_ptr<SnapshotGrid>	m_pgrid_snapshot;

public:
	using Base::reset;
	using MixinType::size;
	using MixinType::offset;

	/// Explicitly defined default constructor.
	PartitionedGridBase () : Base(), MixinType()
	{}

	/**
	 * Constructor
	 *
	 * @param size_ spatial size of the whole grid.
	 * @param offset_ spatial offset of the whole grid.
	 * @param size_partition_ size of each spatial partition.
	 */
	PartitionedGridBase (
		const VecDu& size_, const VecDi& offset_, const LeafType& background_,
		const VecDu& size_partition_
	) {
		this->init(size_, offset_, background_, size_partition_);
	}

	/**
	 * Initialisation function called by non-trivial constructor and
	 * subclasses.
	 *
	 * @param size the dimensions of the whole grid.
	 * @param offset the offset of the whole grid.
	 * @param size_partition the dimensions of each spatial partition.
	 */
	void init (
		const VecDu& size, const VecDi& offset, const LeafType& background_,
		const VecDu& size_partition
	) {
		Base::init(size_partition);
		MixinType::init(size, offset, background_);
	}

	/**
	 * Reshape grid to given size, initialising child grids within the spatial partitions.
	 *
	 * @param size_ size of the overall grid.
	 */
	void size (const VecDu& size_)
	{
		Base::size(size_);

		this->m_size = size_;

		for (UINT idx = 0; idx < this->children().data().size(); idx++)
		{
			Child& child = this->children().data()[idx];
			child.background() = this->m_background;
			child.size(this->m_usize_child);
		}
	}

	/**
	 * Set offset of children grid and propagate offset to children, translating as appropriate.
	 *
	 * @param offset_ offset of overall grid.
	 */
	void offset (const VecDi& offset_)
	{
		Base::offset(offset_);
		MixinType::offset(offset_);

		for (UINT idx = 0; idx < this->children().data().size(); idx++)
		{
			Child& child = this->children().data()[idx];
			const VecDi& pos_child = this->children().index(idx);
			const VecDi& offset_child = (
				(
					(pos_child - this->children().offset()).array()
					* this->m_isize_child.array()
				).matrix() + offset_
			);

			child.offset(offset_child);
		}
	}

	/**
	 * Calculate the position of a child grid (i.e. partition) given the position of leaf grid node.
	 *
	 * @param pos_leaf_
	 * @return position of spatial partition in which leaf position lies.
	 */
	VecDi pos_child (const VecDi& pos_leaf_) const
	{
		return (
			(pos_leaf_ - offset()).array() / this->m_isize_child.array()
		).matrix() + this->children().offset();
	}

	/**
	 * Fill with a single value.
	 *
	 * Loop child grids calling their fill function.
	 *
	 * @param val_ value to fill grid with.
	 */
	void fill (const LeafType& val_)
	{
		for (UINT idx = 0; idx < this->children().data().size(); idx++)
		{
			Child& child = this->children().data()[idx];
			child.fill(val_);
		}
	}

	/**
	 * Get the leaf grid node at pos by navigating to the correct partition.
	 *
	 * Overrides base Grid class's get method (and thus operator()).
	 *
	 * @param pos_ position in grid to fetch.
	 * @return value stored at grid point.
	 */
	LeafType& get (const VecDi& pos_)
	{
		const VecDi& pos_child = this->pos_child(pos_);
		Child& child = this->children().get(pos_child);
		return child.get(pos_);
	}

	/**
	 * Get the leaf grid node at pos by navigating to the correct partition.
	 *
	 * Overrides base Grid class's get method (and thus operator()).
	 *
	 * @param pos_ position in grid to fetch.
	 * @return value stored at grid point pos.
	 */
	const LeafType& get (const VecDi& pos_) const
	{
		const Child& child = this->children()(pos_child(pos_));
		return child.get(pos_);
	}

	/**
	 * Get and store a snapshot of the spatially partitioned data in a contiguous grid.
	 *
	 * Useful for serialization.
	 *
	 * @return non-partitioned contiguous grid containing a copy of the data.
	 */
	SnapshotGrid& snapshot()
	{
		m_pgrid_snapshot.reset(new SnapshotGrid(this->size(), this->offset(), LeafType()));

		const VecDi pos_max = (
			this->m_size.template cast<INT>() - VecDi::Constant(1) + this->offset()
		);

		for (UINT idx = 0; idx <= this->index(pos_max); idx++)
		{
			const VecDi pos = this->index(idx);
			m_pgrid_snapshot->get(pos) = this->get(pos);
		}

		return *m_pgrid_snapshot;
	}

	/**
	 * Copy the snapshot of the grid data back into the partitioned structure.
	 *
	 * Useful for deserialisation.
	 */
	void flush_snapshot()
	{
		if (m_pgrid_snapshot.get() == NULL)
			return;

		for (UINT child_idx = 0; child_idx < this->children().data().size(); child_idx++)
		{
			Child& child = this->children().data()[child_idx];
			for (UINT leaf_idx = 0; leaf_idx < child.data().size(); leaf_idx++)
			{
				const VecDi& pos = child.index(leaf_idx);
				if (m_pgrid_snapshot->inside(pos))
					child.get(pos) = m_pgrid_snapshot->get(pos);
			}
		}
	}
};


/**
 * Standard spatially partitioned grid storing arbitrary data.
 *
 * @tparam T the type of data to store in the grid.
 * @tparam D the dimensions of the grid.
 */
template <typename T, UINT D>
class PartitionedGrid : public PartitionedGridBase<PartitionedGrid<T, D> >
{
public:
	/// Base class.
	using Base = PartitionedGridBase<PartitionedGrid<T, D> >;
	/// Inherited constructor.
	using Base::PartitionedGridBase;
};


/**
 * Traits for PartitionedGridBase.
 *
 * Just forward the traits defined for PartitionedGridBase subclasses.
 *
 * @tparam Derived the CRTP derived class
 */
template <class Derived>
struct GridTraits<PartitionedGridBase<Derived> > : GridTraits<Derived>
{};


/**
 * Traits for PartitionedGrid.
 *
 * @tparam T the type of data to store in the grid.
 * @tparam D the number of dimensions of the grid.
 */
template <typename T, UINT D>
struct GridTraits<PartitionedGrid<T, D> > : DefaultGridTraits<T, D>
{
	/// The class inheriting from the base.
	using ThisType = PartitionedGrid<T, D>;
	/// Child type to store in spatial partitions - in this case a standard Grid.
	using ChildType = Grid<T, D>;
	/// Mixin type whose signature to spoof - in this case a standard GridBase.
	using MixinType = GridBase<ThisType>;
	/// Numer of tracking lists to use, in this case only 1.
	static const UINT NumLists = 1;
};

} // End namespace felt.
#endif /* PARTITIONEDGRID_HPP_ */
