#ifndef PARTITIONEDGRID_HPP_
#define PARTITIONEDGRID_HPP_

#include <memory>
#include <array>
#include <mutex>
#include <boost/iterator/iterator_facade.hpp>

#include "LookupGrid.hpp"
#include "SharedLookupGrid.hpp"
#include "TrackedGrid.hpp"

namespace felt
{

//template <typename T> using AlignedArray = std::vector<T, Eigen::aligned_allocator<T> >;
/// A byte-aligned array of arbitrary type with a mutex member for external thread-safety.
template <typename T>
class AlignedArray :  public std::vector<T, Eigen::aligned_allocator<T> >
{
protected:
	/// Mutex for external locking.
	std::mutex m_mutex;
public:
	/// Base class.
	using Base = std::vector<T, Eigen::aligned_allocator<T> >;

	/// Use constructor of std::vector.
	using Base::vector;

	/// Get mutex member.
	std::mutex& mutex()
	{
		return m_mutex;
	}
};


/** @addtogroup PartitionedGrids
 *
 *  Spatially partitioned versions of Grid, AlignedArray, LookupGrid and TrackedGrid.
 *
 *  @{
 */

/**
 * @defgroup Classes
 * @defgroup Traits
 */
/// Default size of a spatial partition (in each dimension).
static const UINT DEFAULT_PARTITION = 4;

/** @addtogroup Classes
 *  @{
 */


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
	void init (
		const VecDu& size_, const VecDi& offset_,
		const VecDu& partition_size_ = VecDu::Constant(DEFAULT_PARTITION)
	) {
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
		const typename ChildrenGrid::Lookup::Traits::NULL_IDX_TYPE idxs
			= this->children().lookup().get(pos_child_);

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


/**
 * Base class for common features of PartitionedArray template specialisations.
 */
template <class Derived>
class PartitionedArrayBase : public PartitionBase<PartitionedArrayBase<Derived> >
{
public:
	/// This class.
	using ThisType = PartitionedArrayBase<Derived>;
	/// Base class.
	using Base = PartitionBase<ThisType>;
	using typename Base::VecDu;
	using typename Base::VecDi;
protected:
	/// Spatial offset of 'imaginary' grid containing the list.
	VecDi	m_offset;
public:
	/**
	 * Set offset of 'imaginary' grid containing the list.
	 *
	 * @param offset_
	 */
	void offset (const VecDi& offset_)
	{
		m_offset = offset_;
		Base::offset(offset_);
	}

	/**
	 * Get spatial partition from leaf grid node in 'imaginary' grid.
	 *
	 * @param pos_leaf_
	 * @return location of spatial partition in children grid.
	 */
	const VecDi pos_child (const VecDi& pos_leaf_) const
	{
		return (
			(pos_leaf_ - m_offset).array() / this->m_isize_child.array()
		).matrix() + this->children().offset();
	}
};


/**
 * Spatially partitioned expandable lists.
 *
 * A specialised partitioned grid, where the child grids are simply expandable lists.
 *
 * @tparam T the type to store in elements of the list
 * @tparam D the dimension of the 'imaginary' grid.
 * @tparam N the dimension of the array(s).
 */
template <typename T, UINT D, UINT N=0>
class PartitionedArray
	: public PartitionedArrayBase<PartitionedArray<T, D, N> >
{
protected:
	/// This class.
	using ThisType = PartitionedArray<T, D, N>;
	/// Base class
	using Base = PartitionedArrayBase<ThisType>;

	using typename Base::Child;
	using typename Base::VecDu;
	using typename Base::VecDi;

public:
	/// Explicitly defined default constructor.
	PartitionedArray () = default;

	/**
	 * Construct multiple spatially partitioned arrays contained in an 'imaginary' grid.
	 *
	 * @param size_ spatial size of 'imaginary' grid.
	 * @param offset_ spatial offset of 'imaginary' grid.
	 * @param size_partition_ spatial size of a single partition.
	 */
	PartitionedArray (
		const VecDu& size_, const VecDi& offset_,
		const VecDu& size_partition_ = VecDu::Constant(DEFAULT_PARTITION)
	) : Base()
	{
		this->init(size_, offset_, size_partition_);
	}

	/**
	 * Add val to list, placing in partition found from pos.
	 *
	 * @param pos_ position in 'imaginary' grid.
	 * @param val_ value to insert in list.
	 * @param arr_idx_ ID of list to insert into.
	 */
	void add(const VecDi& pos_, const T& val_, const UINT arr_idx_)
	{
		const VecDi& pos_child = this->pos_child(pos_);
		this->children().get(pos_child)[arr_idx_].push_back(val_);
		Base::add_child(pos_child, arr_idx_);
	}

	/**
	 * Loop all spatial partitions, resizing the given list to zero in each.
	 *
	 * @param arr_idx_ ID of list to reset.
	 */
	void reset(const UINT arr_idx_)
	{
		for (const VecDi& pos_child : this->children().list(arr_idx_))
			this->children().get(pos_child)[arr_idx_].clear();
		Base::reset(arr_idx_);
	}
};


/**
 * Spatially partitioned expandable list - 1D array specialisation.
 *
 * A specialised partitioned grid, where the child grids are simply
 * expandable lists.
 *
 * @tparam T type to store in elements of the list
 * @tparam D dimension of the 'imaginary' grid.
 */
template <typename T, UINT D>
class PartitionedArray<T, D, 0>
	: public PartitionedArrayBase<PartitionedArray<T, D, 0> >
{
protected:
	/// This class
	using ThisType = PartitionedArray<T, D, 0>;
	/// Base class.
	using Base = PartitionedArrayBase<ThisType>;

	using typename Base::Child;
	using typename Base::VecDu;
	using typename Base::VecDi;

public:
	/// Explicitly defined default constructor.
	PartitionedArray () = default;

	/**
	 * Construct a single (1D) spatially partitioned array contained in an 'imaginary' grid.
	 *
	 * @param size_ spatial size of 'imaginary' grid.
	 * @param offset_ spatial offset of 'imaginary' grid.
	 * @param size_partition_ spatial size of a single partition.
	 */
	PartitionedArray (
		const VecDu& size_, const VecDi& offset_,
		const VecDu& size_partition_ = VecDu::Constant(DEFAULT_PARTITION)
	) : Base()
	{
		this->init(size_, offset_, size_partition_);
	}

	/**
	 * Add val to list, placing in partition found from pos.
	 *
	 * @param pos position in 'imaginary' grid.
	 * @param val value to insert in list.
	 */
	void add(const VecDi& pos, const T& val)
	{
		const VecDi& pos_child = this->pos_child(pos);
		this->children().get(pos_child).push_back(val);
		Base::add_child(pos_child);
	}

	/**
	 * Thread safely add val to list, placing in partition found from pos.
	 *
	 * @param pos position in 'imaginary' grid.
	 * @param val value to insert in list.
	 */
	void add_safe(const VecDi& pos, const T& val)
	{
		const VecDi& pos_child = this->pos_child(pos);
		Child& child = this->children().get(pos_child);
		Base::add_child(pos_child);
		std::lock_guard<std::mutex> lock(child.mutex());
		child.push_back(val);
	}

	/**
	 * Loop all spatial partitions, resizing their lists to zero.
	 */
	void reset()
	{
		for (const VecDi& pos_child : this->children().list())
			this->children().get(pos_child).clear();
		Base::reset();
	}
};


/**
 * Base class for spatially partitioned grid storing arbitrary values.
 *
 * Uses multiple inheritance from child grid (MixinType) class, spoofing the signature,
 * and PartitionBase for actual storage in a (single-level) tree hierarchy,
 * i.e. PartitionedGridBase -> Branch -> Child -> leaf.
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
	 * TrackedGrid type storing spatial partitions.
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
	PartitionedGridBase () = default;

	/**
	 * Constructor
	 *
	 * @param size_ spatial size of the whole grid.
	 * @param offset_ spatial offset of the whole grid.
	 * @param size_partition_ size of each spatial partition.
	 */
	PartitionedGridBase (
		const VecDu& size_, const VecDi& offset_,
		const VecDu& size_partition_ = VecDu::Constant(DEFAULT_PARTITION)
	) : Base()
	{
		self->init(size_, offset_, size_partition_);
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
		const VecDu& size, const VecDi& offset,
		const VecDu& size_partition = VecDu::Constant(DEFAULT_PARTITION)
	) {
		Base::init(size_partition);
		MixinType::init(size, offset);
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
			Child& child = this->children().data()(idx);
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
			Child& child = this->children().data()(idx);
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
	const VecDi pos_child (const VecDi& pos_leaf_) const
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
			Child& child = this->children().data()(idx);
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
	typename SnapshotGrid::ArrayData& data()
	{
		m_pgrid_snapshot.reset(new SnapshotGrid(this->size(), this->offset()));

		for (UINT branch_idx = 0; branch_idx < this->children().data().size(); branch_idx++)
		{
			const Child& child = this->children().data()[branch_idx];
			for (UINT leaf_idx = 0; leaf_idx < child.data().size(); leaf_idx++)
			{
				const VecDi& pos = child.index(leaf_idx);
				if (m_pgrid_snapshot->inside(pos))
					m_pgrid_snapshot->get(pos) = child.get(pos);
			}
		}
		return m_pgrid_snapshot->data();
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

		for (UINT branch_idx = 0; branch_idx < this->children().data().size(); branch_idx++)
		{
			Child& child = this->children().data()[branch_idx];
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

// End group Classes.
/** @}
 *  @addtogroup Traits
 *  @{
 */

/**
 * Traits for PartitionedArrayBase.
 *
 * Just forward the traits defined for PartitionedArrayBase subclasses.
 */
template <class Derived>
struct GridTraits<PartitionedArrayBase<Derived> > : GridTraits<Derived>
{};


/**
 * Traits for PartitionedArray.
 *
 * @tparam T type to store in elements of the list
 * @tparam D dimension of the 'imaginary' grid.
 * @tparam N number of distinct arrays to partition.
 */
template <typename T, UINT D, UINT N>
struct GridTraits<PartitionedArray<T, D, N> > : DefaultGridTraits<T, D>
{
	/// The class inheriting from the base.
	using ThisType = PartitionedArray<T, D, N>;
	/// Child type to store in spatial partitions - in this case an array of lists.
	using ChildType = std::array<AlignedArray<T>, N>;
	/// Number of distinct arrays to partition - in this case a single list.
	static const UINT NumLists = N;
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
 * Traits for PartitionedArray.
 *
 * @tparam T data type stored in array.
 * @tparam D dimension of the grid.
 */
template <typename T, UINT D>
struct GridTraits<PartitionedArray<T, D, 0> > : DefaultGridTraits<T, D>
{
	/// The class inheriting from the base.
	using ThisType = PartitionedArray<T, D, 0>;
	/// Child type to store in spatial partitions - in this case a single list.
	using ChildType = AlignedArray<T>;
	/// Number of distinct arrays to partition - in this case a single list.
	static const UINT NumLists = 1;
};


/**
 * Traits for GridBase to understand PartitionedGrid.
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

/**
 * @}
 * @} */ // End group Traits  // End group ParitionedGrid.

} // End namespace felt.
#endif /* PARTITIONEDGRID_HPP_ */
