#ifndef PARTITIONEDGRID_HPP_
#define PARTITIONEDGRID_HPP_

#include <memory>
#include <array>
#include <mutex>
#include <boost/iterator/iterator_facade.hpp>

#include "MappedGrid.hpp"

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
 * @defgroup Classes Partition grid classes.
 * @defgroup Traits Traits for CRTP static inheritance.
 */
/// Default size of a spatial partition (in each dimension).
static const UINT DEFAULT_PARTITION = 4;

/**
 * @ingroup Traits
 * Base traits class for classes CRTP derived from PartitionBase.
 */
template <class Derived> struct PartitionBaseTraits {};


/**
 * @ingroup Classes
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

	/// Child object to store in each partition.
	using Child = typename PartitionBaseTraits<Derived>::ChildType;

	/// Number of tracking lists of points.
	static const UINT NUM_LISTS = PartitionBaseTraits<Derived>::NumLists;

	/// Grid of partitions with tracking list(s) of active partitions.
	using BranchGrid = TrackedGrid<Child, PartitionBaseTraits<Derived>::Dims, NUM_LISTS>;
	/// D-dimensional unsigned int vector.
	using VecDu = typename BranchGrid::VecDu;
	/// D-dimensional signed int vector.
	using VecDi = typename BranchGrid::VecDi;


protected:
	/// Grid of partitions with tracking list(s) of active grid points.
	BranchGrid	m_grid_branch;

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
		: m_grid_branch()
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
	 * Get TrackedGrid branch grid - the spatial partition grid that stores the PartitionBase::Child 
	 * objects.
	 *
	 * @return TrackedGrid storing/tracking Child objects.
	 */
	BranchGrid& branch ()
	{
		return m_grid_branch;
	}

	/**
	 * Get TrackedGrid branch grid - the spatial partition grid that stores the PartitionBase::Child 
	 * objects.
	 *
	 * @return TrackedGrid storing/tracking Child objects.
	 */
	const BranchGrid& branch () const
	{
		return m_grid_branch;
	}

	/**
	 * Get the PartitionBase::Child object at given position.
	 *
	 * Shorthand for branch().get(pos) or branch()(pos).
	 *
	 * @param pos position within the spatial partition grid.
	 * @return PartitionBase::Child at given position.
	 */
	Child& child (const VecDi& pos_)
	{
		return m_grid_branch(pos_);
	}

	/**
	 * Get the PartitionBase::Child object at given position.
	 *
	 * Shorthand for branch().get(pos) or branch()(pos).
	 *
	 * @param pos position within the spatial partition grid.
	 * @return PartitionBase::Child at given position.
	 */
	const Child& child (const VecDi& pos_) const
	{
		return m_grid_branch(pos_);
	}

	/**
	 * Reshape grid, computing the size of the branch grid.
	 *
	 * The branch grid will be increased in size by one, if required, to ensure all leaf nodes are
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

		m_grid_branch.size(branch_size);
	}

	/**
	 * Calculate the offset of the branch grid from the given offset and the size of a spatial
	 * partition (the child size).
	 *
	 * @param grid_offset_ overall spatial offset the grid
	 */
	void offset (const VecDi& grid_offset_)
	{
		const VecDi& branch_offset = (
			grid_offset_.array() / m_isize_child.array()
		).matrix();

		m_grid_branch.offset(branch_offset);
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
		return this->branch().add(pos_, arr_idx_);
	}

	/**
	 * Remove a spatial partition from branch grid's tracking subgrid.
	 *
	 * Uses mutex for thread safety.
	 *
	 * @param pos_ position of spatial partition to stop tracking.
	 * @param arr_idx_ index of tracking list used to track position.
	 */
	void remove_child(const VecDi& pos_, const UINT& arr_idx_ = 0)
	{
		if (!m_grid_branch.is_active(pos_, arr_idx_))
			return;
		std::lock_guard<std::mutex> lock(m_mutex_update_branch);
		this->branch().remove(pos_, arr_idx_);
	}

	/**
	 * Reset the tracking list at given index in the branch grid.
	 *
	 * Removes all spatial partitions from the tracking subgrid for given list index.
	 *
	 * @param arr_idx_ index of tracking list to reset.
	 */
	void reset(const UINT& arr_idx_ = 0)
	{
		this->branch().reset(arr_idx_);
	}
};


/**
 * @ingroup Traits
 * Base traits class for classes CRTP derived from PartitionedArrayBase.
 *
 * @tparam Derived the CRTP derived class
 */
template <class Derived> struct PartitionedArrayBaseTraits {};


/**
 * @ingroup Classes
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
	 * @return location of spatial partition in branch grid.
	 */
	const VecDi pos_child (const VecDi& pos_leaf_) const
	{
		return (
			(pos_leaf_ - m_offset).array() / this->m_isize_child.array()
		).matrix() + this->branch().offset();
	}
};


/**
 * @ingroup Traits
 * Traits for PartitionBase to understand PartitionedArrayBase.
 *
 * Just forward the traits defined for PartitionedArrayBase subclasses.
 */
template <class Derived>
struct PartitionBaseTraits<PartitionedArrayBase<Derived> >
{
	/// The class inheriting from the base.
	using ThisType = typename PartitionedArrayBaseTraits<Derived>::ThisType;
	/// Child type to store in spatial partitions.
	using ChildType = typename PartitionedArrayBaseTraits<Derived>::ChildType;
	/// Dimensions of the grid.
	static const UINT Dims = PartitionedArrayBaseTraits<Derived>::Dims;
	/// Number of tracking lists.
	static const UINT NumLists = PartitionedArrayBaseTraits<Derived>::NumLists;
};


/**
 * @ingroup Classes
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
	void add(const VecDi& pos_, const T& val_, const UINT& arr_idx_)
	{
		const VecDi& pos_child = this->pos_child(pos_);
		this->child(pos_child)[arr_idx_].push_back(val_);
		Base::add_child(pos_child, arr_idx_);
	}

	/**
	 * Loop all spatial partitions, resizing the given list to zero in each.
	 *
	 * @param arr_idx_ ID of list to reset.
	 */
	void reset(const UINT& arr_idx_)
	{
		for (const VecDi& pos_child : this->branch().list(arr_idx_))
			this->child(pos_child)[arr_idx_].clear();
		Base::reset(arr_idx_);
	}
};


/**
 * @ingroup Traits
 * Traits of PartitionedArray for CRTP inheritance from PartitionedArrayBase.
 *
 * @tparam T type to store in elements of the list
 * @tparam D dimension of the 'imaginary' grid.
 * @tparam N number of distinct arrays to partition.
 */
template <typename T, UINT D, UINT N>
struct PartitionedArrayBaseTraits<PartitionedArray<T, D, N> >
{
	/// The class inheriting from the base.
	using ThisType = PartitionedArray<T, D, N>;
	/// Child type to store in spatial partitions - in this case an array of lists.
	using ChildType = std::array<AlignedArray<T>, N>;
	/// Dimensions of the grid, from template parameter.
	static const UINT Dims = D;
	/// Number of distinct arrays to partition, from template parameter.
	static const UINT NumLists = N;
};


/**
 * @ingroup Classes
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
		this->child(pos_child).push_back(val);
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
		Child& child = this->child(pos_child);
		Base::add_child(pos_child);
		std::lock_guard<std::mutex> lock(child.mutex());
		child.push_back(val);
	}

	/**
	 * Loop all spatial partitions, resizing their lists to zero.
	 */
	void reset()
	{
		for (const VecDi& pos_child : this->branch().list())
			this->child(pos_child).clear();
		Base::reset();
	}
};

/**
 * @ingroup Traits
 * Traits of 1D PartitionedArray for CRTP inheritance from PartitionedArrayBase.
 *
 * @tparam T data type stored in array.
 * @tparam D dimension of the grid.
 */
template <typename T, UINT D>
struct PartitionedArrayBaseTraits<PartitionedArray<T, D, 0> >
{
	/// The class inheriting from the base.
	using ThisType = PartitionedArray<T, D, 0>;
	/// Child type to store in spatial partitions - in this case a single list.
	using ChildType = AlignedArray<T>;
	/// Size of the grid, from template parameter.
	static const UINT Dims = D;
	/// Number of distinct arrays to partition - in this case a single list.
	static const UINT NumLists = 1;
};

/**
 * @ingroup Traits
 * Base traits class for classes CRTP derived from PartitionedGridBase.
 */
template <class Derived> struct PartitionedGridBaseTraits {};


/**
 * @ingroup Classes
 * Base class for spatially partitioned grid storing arbitrary values.
 *
 * Uses multiple inheritance from child grid (MixinType) class, spoofing the signature,
 * and PartitionBase for actual storage in a (single-level) tree hierarchy,
 * i.e. PartitionedGridBase -> Branch -> Child -> leaf.
 */
template <class Derived>
class PartitionedGridBase
	: public PartitionedGridBaseTraits<Derived>::MixinType,
	  public PartitionBase<PartitionedGridBase<Derived> >
{
public:
	/// CRTP derived class.
	using DerivedType = Derived;
	using ThisType = PartitionedGridBase<DerivedType>;
	/// PartitionBase base class.
	using Base = PartitionBase<ThisType>;
	/// Base grid class to partition.
	using MixinType = typename PartitionedGridBaseTraits<Derived>::MixinType;
	/// Traits of CRTP derived class
	using Traits = PartitionedGridBaseTraits<DerivedType>;
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
	using BranchGrid = typename Base::BranchGrid;
	
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

		for (UINT idx = 0; idx < this->branch().data().size(); idx++)
		{
			Child& child = this->branch().data()(idx);
			child.size(this->m_usize_child);
		}
	}

	/**
	 * Set offset of branch grid and propagate offset to children, translating as appropriate.
	 *
	 * @param offset_ offset of overall grid.
	 */
	void offset (const VecDi& offset_)
	{
		Base::offset(offset_);
		MixinType::offset(offset_);

		for (UINT idx = 0; idx < this->branch().data().size(); idx++)
		{
			Child& child = this->branch().data()(idx);
			const VecDi& pos_child = this->branch().index(idx);
			const VecDi& offset_child = (
				(
					(pos_child - this->branch().offset()).array()
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
		).matrix() + this->branch().offset();
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
		for (UINT idx = 0; idx < this->branch().data().size(); idx++)
		{
			Child& child = this->branch().data()(idx);
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
		Child& child = this->branch().get(pos_child);
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
		const Child& child = this->branch()(pos_child(pos_));
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

		for (UINT branch_idx = 0; branch_idx < this->branch().data().size(); branch_idx++)
		{
			const Child& child = this->branch().data()[branch_idx];
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

		for (UINT branch_idx = 0; branch_idx < this->branch().data().size(); branch_idx++)
		{
			Child& child = this->branch().data()[branch_idx];
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
 * @ingroup Traits
 * Traits for PartitionBase to understand PartitionedGridBase.
 *
 * Just forward the traits defined for PartitionedGridBase subclasses.
 *
 * @tparam Derived the CRTP derived class
 */
template <class Derived>
struct PartitionBaseTraits<PartitionedGridBase<Derived> >
{
	/// The class inheriting from the base.
	using ThisType = typename PartitionedGridBaseTraits<Derived>::ThisType;
	/// Child type to store in spatial partitions.
	using ChildType = typename PartitionedGridBaseTraits<Derived>::ChildType;
	/// Size of the grid.
	static const UINT Dims = PartitionedGridBaseTraits<Derived>::Dims;
	/// Number of distinct tracking lists to track active spatial partitions.
	static const UINT NumLists = PartitionedGridBaseTraits<Derived>::NumLists;
};


/**
 * @ingroup Classes
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
 * @ingroup Traits
 * Traits for GridBase to understand PartitionedGrid.
 *
 * @tparam T the type of data to store in the grid.
 * @tparam D the number of dimensions of the grid.
 */
template <typename T, UINT D>
struct GridBaseTraits<PartitionedGrid<T, D> >
{
	/// The class inheriting from the base.
	using ThisType = PartitionedGrid<T, D>;
	/// Dimensions of the grid, from template parameter.
	static const UINT Dims = D;
	/// The data type to store at leaf grid nodes, from template parameter.
	using LeafType = T;
};


/**
 * @ingroup Traits
 * Traits for PartitionedGridBase to understand PartitionedGrid.
 *
 * @tparam T the type of data to store in the grid.
 * @tparam D the number of dimensions of the grid.
 */
template <typename T, UINT D>
struct PartitionedGridBaseTraits<PartitionedGrid<T, D> >
{
	/// The class inheriting from the base.
	using ThisType = PartitionedGrid<T, D>;
	/// The data type to store at leaf grid nodes, from template parameter.
	using LeafType = T;
	/// Child type to store in spatial partitions - in this case a standard Grid.
	using ChildType = Grid<T, D>;
	/// Mixin type whose signature to spoof - in this case a standard GridBase.
	using MixinType = GridBase<ThisType>;
	/// Dimensions of the grid, from template parameter.
	static const UINT Dims = D;
	/// Numer of tracking lists to use, in this case only
	static const UINT NumLists = 1;
};


/**
 * @ingroup Classes
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
			const GridTree* pgrid, const UINT& list_idx,
			const Iter& it_child, const Iter& it_leaf
		)
		: m_pgrid(pgrid), m_list_idx(list_idx), m_it_child(it_child),
		  m_it_leaf(it_leaf),
		  m_it_child_end(pgrid->branch().list(list_idx).end())
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
				m_it_leaf = m_pgrid->child(*m_it_child).list(m_list_idx).begin();
				return;
			}

			m_it_leaf++;

			if (m_it_leaf == m_pgrid->child(*m_it_child).list(m_list_idx).end())
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
	LeafsContainer(const GridTree* pgrid, const UINT& list_idx)
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
			= m_pgrid->branch().list(m_list_idx).cbegin();
		const typename PosArray::const_iterator& it_child_end
			= m_pgrid->branch().list(m_list_idx).cend();

		typename PosArray::const_iterator it_leaf_begin = it_child_end;

		if (it_child_begin != it_child_end)
		{
			it_leaf_begin = m_pgrid->child(
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
			m_pgrid->branch().list(m_list_idx).cend(),
			m_pgrid->branch().list(m_list_idx).cend()
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
		for (const VecDi& pos_child : m_pgrid->branch().list(m_list_idx))
			sum += m_pgrid->child(pos_child).list(m_list_idx).size();
		return sum;
	}
};


/**
 * @ingroup Traits
 * Base traits for spatially partitioned wrapper for lookup and tracked grid.
 *
 * @tparam Derived the CRTP derived class
 */
template <class Derived> struct TrackingPartitionedGridTraits {};


/**
 * @ingroup Classes
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
	/// Base class.
	using Base = PartitionedGridBase<ThisType>;

	friend class LeafsContainer<ThisType>;

	using typename Base::Child;
	using typename Base::VecDu;
	using typename Base::VecDi;
	using typename Base::BranchGrid;

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
	void reset(const UINT& arr_idx_ = 0)
	{
		BranchGrid& branch = this->branch();
		for (const VecDi& pos_child : branch.list(arr_idx_))
			branch(pos_child).reset(arr_idx_);

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
	bool add(const VecDi& pos_, const UINT& arr_idx_ = 0)
	{
		const VecDi& pos_child = this->pos_child(pos_);
		Base::add_child(pos_child, arr_idx_);
		return this->child(pos_child).add(pos_, arr_idx_);
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
	bool add_safe(const VecDi& pos_, const UINT& arr_idx_ = 0)
	{
		const VecDi& pos_child = this->pos_child(pos_);
		Base::add_child(pos_child, arr_idx_);
		Child& child = this->child(pos_child);
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
	void remove(const VecDi& pos_, const UINT& arr_idx_ = 0)
	{
		const VecDi& pos_child = this->pos_child(pos_);
		Child& child = this->child(pos_child);
		child.remove(pos_, arr_idx_);
		if (child.list(arr_idx_).size() == 0)
			this->remove_child(pos_child, arr_idx_);
	}

	/**
	 * Return structure for range based for loops over leaf nodes.
	 *
	 * @param list_idx_ tracking list id.
	 * @return
	 */
	const LeafsContainer<Derived> leafs(const UINT& list_idx_ = 0) const
	{
		return LeafsContainer<Derived>(
			static_cast<const Derived*>(this), list_idx_
		);
	}

private:
	/**
	 * Override spoofed (non-partitioned) base class's list method to make
	 * private, since it will always be empty (must use branch or child).
	 *
	 * @param arr_idx_
	 * @return list of this grid (will always be empty).
	 */
	PosArray& list(const UINT& arr_idx_ = 0)
	{
		return Child::list(arr_idx_);
	}
};


/**
 * @ingroup Traits
 * Traits for PartitionedGridBase to understand TrackingPartitionedGridBase.
 *
 * Just forward the traits defined for TrackingPartitionedGridBase subclasses.
 *
 * @tparam Derived the CRTP derived class.
 */
template <class Derived>
struct PartitionedGridBaseTraits<TrackingPartitionedGridBase<Derived> >
{
	/// The class inheriting from the base.
	using ThisType = typename TrackingPartitionedGridTraits<Derived>::ThisType;
	/// The data type to store at leaf grid nodes.
	using LeafType = typename TrackingPartitionedGridTraits<ThisType>::LeafType;
	/// Child grid class.  Each spatial partition contains one Child grid.
	using ChildType = typename TrackingPartitionedGridTraits<ThisType>::ChildType;
	/// Base grid class to partition.
	using MixinType = typename TrackingPartitionedGridTraits<ThisType>::MixinType;
	/// Dimensions of the grid.
	static const UINT Dims = TrackingPartitionedGridTraits<ThisType>::Dims;
	/// Number of tracking lists.
	static const UINT NumLists = TrackingPartitionedGridTraits<ThisType>::NumLists;
};


/**
 * @ingroup Classes
 * Spatially partitioned wrapper for TrackedGrid.
 *
 * @tparam T type of value stored in leaf grid nodes.
 * @tparam D dimension of the grid.
 * @tparam N number of tracking lists to use.
 */
template <typename T, UINT D, UINT N>
class TrackedPartitionedGrid
	: public TrackingPartitionedGridBase <TrackedPartitionedGrid<T, D, N>>
{
public:
	/// This class.
	using ThisType = TrackedPartitionedGrid<T, D, N>;
	/// Base class
	using Base = TrackingPartitionedGridBase<ThisType>;

	using typename Base::Child;
	using typename Base::VecDu;
	using typename Base::VecDi;
	using typename Base::BranchGrid;

	using PosArray = typename Child::PosArray;
public:
	using Base::TrackingPartitionedGridBase;

	/**
	 * Set value in grid at given position and add position to lookup grid.
	 *
	 * Will set value regardless whether lookup grid already set for given
	 * position + tracking list.
	 *
	 * Descend to relevant child grid to add to their tracking structure.
	 *
	 * @param pos_ position in grid.
	 * @param val_ value to set.
	 * @param arr_idx_ tracking list id.
	 * @return true if grid node set in child lookup grid and position added
	 * to tracking list, false if child grid node was already set so
	 * position already in a list.
	 */
	bool add(const VecDi& pos_, const T& val_, const UINT& arr_idx_ = 0)
	{
		const VecDi& pos_child = this->pos_child(pos_);
		Base::add_child(pos_child, arr_idx_);
		return this->child(pos_child).add(pos_, val_, arr_idx_);
	}

	/**
	 * Thread-safely set value in grid at given position and add position to lookup grid.
	 *
	 * Will set value regardless whether lookup grid already set for given
	 * position + tracking list.
	 *
	 * Descend to relevant child grid to add to their tracking structure.
	 *
	 * @param pos_ position in grid.
	 * @param val_ value to set.
	 * @param arr_idx_ tracking list id.
	 * @return true if grid node set in child lookup grid and position added
	 * to tracking list, false if child grid node was already set so
	 * position already in a list.
	 */
	bool add_safe(const VecDi& pos_, const T& val_, const UINT& arr_idx_ = 0)
	{
		const VecDi& pos_child = this->pos_child(pos_);
		Base::add_child(pos_child, arr_idx_);
		Child& child = this->child(pos_child);
		std::lock_guard<std::mutex> lock(child.mutex());
		return child.add(pos_, val_, arr_idx_);
	}

	/**
	 * Set every active grid node (i.e. those referenced by lookup grid)
	 * to given value and reset the lookup grid.
	 *
	 * Lookup grid will then be full of NULL indices for given tracking list
	 * and the relevant tracking list(s) will be empty.
	 *
	 * @param val_ value to set in main grid.
	 * @param arr_idx_ tracking list id to cycle over and clear.
	 */
	void reset(const T& val_, const UINT& arr_idx_ = 0)
	{
		BranchGrid& branch = this->branch();
		for (const VecDi& pos_child : branch.list(arr_idx_))
			branch(pos_child).reset(val_, arr_idx_);

		Base::reset(arr_idx_);
	}

	/**
	 * Reset a tracking list on the lookup grid.
	 *
	 * @param arr_idx_ tracking list to clear.
	 */
	void reset(const UINT& arr_idx_ = 0)
	{
		BranchGrid& branch = this->branch();
		for (const VecDi& pos_child : branch.list(arr_idx_))
			branch(pos_child).reset(arr_idx_);

		Base::reset(arr_idx_);
	}
};


/**
 * @ingroup Traits
 * Traits for TrackedGridBase to understand TrackedPartitionedGrid.
 *
 * @tparam T type of value stored in leaf grid nodes.
 * @tparam D dimension of the grid.
 * @tparam N number of tracking lists to use.
 */
template <typename T, UINT D, UINT N>
struct TrackedGridBaseTraits<TrackedPartitionedGrid<T, D, N> >
{
	/// The class inheriting from the base.
	using ThisType = TrackedPartitionedGrid<T, D, N>;
	/// The type of lookup grid to use for tracking active grid nodes.
	using LookupType = LookupGrid<D, N>;
	/// The data type to store at leaf grid nodes, from template parameter.
	using LeafType = T;
	/// Dimensions of the grid, from template parameter.
	static const UINT Dims = D;
};


/**
 * @ingroup Traits
 * Traits for TrackingPartitionedGrid to understand TrackedPartitionedGrid.
 *
 * @tparam T type of value stored in leaf grid nodes.
 * @tparam D dimension of the grid.
 * @tparam N number of tracking lists to use.
 */
template <typename T, UINT D, UINT N>
struct TrackingPartitionedGridTraits<TrackedPartitionedGrid<T, D, N> >
{
	/// The class inheriting from the base.
	using ThisType = TrackedPartitionedGrid<T, D, N>;
	/// Child grid class, in this case a TrackedGrid.
	using ChildType = TrackedGrid<T, D, N>;
	/// Base grid class to partition, in this case TrackedGridBase.
	using MixinType = TrackedGridBase<ThisType>;
	/// The data type to store at leaf grid nodes, from template parameter.
	using LeafType = T;
	/// Dimensions of the grid, from template parameter.
	static const UINT Dims = D;
	/// Number of tracking lists, from template parameter.
	static const UINT NumLists = N;
};

/**
 * @ingroup Classes
 * Spatially partitioned wrapper for SharedTrackedGrid.
 *
 * @tparam T type of value stored in leaf grid nodes.
 * @tparam D dimension of the grid.
 * @tparam N number of tracking lists to use.
 */
template <typename T, UINT D, UINT N>
class SharedTrackedPartitionedGrid
	: public TrackingPartitionedGridBase<SharedTrackedPartitionedGrid<T, D, N> >
{
public:
	using ThisType = SharedTrackedPartitionedGrid<T, D, N>;
	using Base = TrackingPartitionedGridBase<ThisType>;
	using typename Base::Child;
	using typename Base::VecDu;
	using typename Base::VecDi;
	using typename Base::BranchGrid;
	using PosArray = typename Child::PosArray;

public:
	using Base::TrackingPartitionedGridBase;

	/**
	 * Set value in grid at given position and add position to lookup grid.
	 *
	 * Will set value regardless whether lookup grid already set for given
	 * position + tracking list.
	 *
	 * Descend to relevant child grid to add to their tracking structure.
	 *
	 * @param pos_ position in grid.
	 * @param val_ value to set.
	 * @param arr_idx_ tracking list id.
	 * @return true if grid node set in child lookup grid and position added
	 * to tracking list, false if child grid node was already set so
	 * position already in a list.
	 */
	bool add(const VecDi& pos_, const T& val_, const UINT& arr_idx_)
	{
		const VecDi& pos_child = this->pos_child(pos_);
		Base::add_child(pos_child, arr_idx_);
		return this->child(pos_child).add(pos_, val_, arr_idx_);
	}


	/**
	 * Add a position to the lookup grid.
	 *
	 * @param pos_ position in the grid to add.
	 * @param arr_idx_ tracking list id.
	 * @return true if grid node set in child lookup grid and position added
	 * to tracking list, false if leaf grid node was already set so position
	 * already in a list.
	 */
	bool add(const VecDi& pos_, const UINT& arr_idx_)
	{
		return Base::add(pos_, arr_idx_);
	}

	/**
	 * Set every active grid node (i.e. those referenced by lookup grid)
	 * to given value and reset the lookup grid.
	 *
	 * Lookup grid will then be full of NULL indices for given tracking list
	 * and the relevant tracking list(s) will be empty.
	 *
	 * @param val_ value to set in main grid.
	 * @param arr_idx_ tracking list id to cycle over and clear.
	 */
	void reset(const T& val_, const UINT& arr_idx_ = 0)
	{
		BranchGrid& branch = this->branch();
		for (const VecDi& pos_child : branch.list(arr_idx_))
			branch(pos_child).reset(val_, arr_idx_);

		Base::reset(arr_idx_);
	}
};


/**
 * @ingroup Traits
 * Traits class for TrackedGridBase to understand SharedTrackedPartitionedGrid.
 *
 * @tparam T type of value stored in leaf grid nodes.
 * @tparam D dimension of the grid.
 * @tparam N number of tracking lists to use.
 */
template <typename T, UINT D, UINT N>
struct TrackedGridBaseTraits<SharedTrackedPartitionedGrid<T, D, N> >
{
	/// The class inheriting from TrackedGridBase.
	using ThisType = SharedTrackedPartitionedGrid<T, D, N>;
	/// The type of lookup grid to use for tracking active grid nodes.
	using LookupType = SharedLookupGrid<D, N>;
	/// Lookup value at leaf nodes, from template parameter.
	using LeafType = T;
	/// Dimensions of the grid, from template parameter.
	static const UINT Dims = D;
};


/**
 * @ingroup Traits
 * Traits class for TrackingPartitionedGridBase to understand SharedTrackedPartitionedGrid.
 *
 * @tparam T type of value stored in leaf grid nodes.
 * @tparam D dimension of the grid.
 * @tparam N number of tracking lists to use.
 */
template <typename T, UINT D, UINT N>
struct TrackingPartitionedGridTraits<SharedTrackedPartitionedGrid<T, D, N> >
{
	/// The class inheriting from TrackingPartitionedGrid.
	using ThisType = SharedTrackedPartitionedGrid<T, D, N>;
	/// Child grid class, in this case SharedTrackedGrid.
	using ChildType = SharedTrackedGrid<T, D, N>;
	/// Base grid class to partition, in this case TrackedGridBase.
	using MixinType = TrackedGridBase<ThisType>;
	/// Lookup value at leaf nodes, from template parameter.
	using LeafType = T;
	/// Dimensions of the grid, from template parameter.
	static const UINT Dims = D;
	/// Number of tracking lists, from template parameter.
	static const UINT NumLists = N;
};


/**
 * @ingroup Classes
 * Spatially partitioned wrapper for LookupGrid.
 *
 * @tparam D dimension of the grid.
 * @tparam N number of tracking lists to use.
 */
template <UINT D, UINT N=1>
class LookupPartitionedGrid
	: public TrackingPartitionedGridBase<LookupPartitionedGrid<D, N> >
{
public:
	using ThisType = LookupPartitionedGrid<D, N>;
	using Base = TrackingPartitionedGridBase<ThisType>;
	using Base::TrackingPartitionedGridBase;
};


/**
 * @ingroup Traits
 * Traits class for TrackingPartitionedGridBase to understand LookupPartitionedGrid.
 *
 * @tparam D dimension of the grid.
 * @tparam N number of tracking lists to use.
 */
template <UINT D, UINT N>
struct TrackingPartitionedGridTraits<LookupPartitionedGrid<D, N> >
{
	/// Dimensions of the grid, from template parameter.
	static const UINT Dims = D;
	/// Number of tracking lists, from template parameter.
	static const UINT NumLists = N;
	/// The class inheriting from TrackingPartitionedGrid.
	using ThisType = LookupPartitionedGrid<Dims, NumLists>;
	/// Lookup value at leaf nodes.  N-dimensional array, one element for each tracking list.
	using LeafType = Eigen::Matrix<UINT, NumLists, 1>;
	/// Child grid class, in this case LookupGrid.
	using ChildType = LookupGrid<Dims, NumLists>;
	/// Base grid class to partition, in this case LookupGridBase.
	using MixinType = LookupGridBase<ThisType>;
};


/**
 * @ingroup Traits
 * Traits class for LookupGridBase to understand LookupPartitionedGrid.
 *
 * Just copy relevant traits from TrackingPartitionedGridTraits.
 *
 * @tparam D dimension of the grid.
 * @tparam N number of tracking lists to use.
 */
template <UINT D, UINT N>
struct LookupGridBaseTraits<LookupPartitionedGrid<D, N> >
{
	/// The class inheriting from LookupGridBase.
	using ThisType = LookupPartitionedGrid<D, N>;
	/// Lookup value at leaf nodes.
	using LeafType = typename TrackingPartitionedGridTraits<ThisType>::LeafType;
	/// Dimensions of the grid.
	static const UINT Dims = TrackingPartitionedGridTraits<ThisType>::Dims;
	/// Number of tracking lists to use.
	static const UINT NumLists = TrackingPartitionedGridTraits<ThisType>::NumLists;
};


/**
 * @ingroup Classes
 * Spatially partitioned wrapper for SharedLookupGrid.
 *
 * @tparam D dimension of the grid.
 * @tparam N number of tracking lists to use.
 */
template <UINT D, UINT N>
class SharedLookupPartitionedGrid
	: public TrackingPartitionedGridBase<SharedLookupPartitionedGrid<D, N> >
{
public:
	using ThisType = SharedLookupPartitionedGrid<D, N>;
	using Base = TrackingPartitionedGridBase<ThisType>;
	using Base::TrackingPartitionedGridBase;
};


/**
 * @ingroup Traits
 * Traits class for TrackingPartitionedGridBase to understand SharedLookupPartitionedGrid.
 *
 * @tparam D dimension of the grid.
 * @tparam N number of tracking lists to use.
 */
template <UINT D, UINT N>
struct TrackingPartitionedGridTraits<SharedLookupPartitionedGrid<D, N> >
{
	/// The class inheriting from TrackingPartitionedGrid.
	using ThisType = SharedLookupPartitionedGrid<D, N>;
	/// Lookup value at leaf nodes.
	using LeafType = UINT;
	/// Child grid class, in this case SharedLookupGrid.
	using ChildType = SharedLookupGrid<D, N>;
	using MixinType = SharedLookupGridBase<ThisType>;
	static const UINT Dims = D;
	static const UINT NumLists = N;
};


/**
 * @ingroup Traits
 * Traits class for SharedLookupGridBase to understand SharedLookupPartitionedGrid.
 *
 * Just copy relevant traits from TrackingPartitionedGridTraits.
 *
 * @tparam D dimension of the grid.
 * @tparam N number of tracking lists to use.
 */
template <UINT D, UINT N>
struct SharedLookupGridBaseTraits<SharedLookupPartitionedGrid<D, N> >
{
	/// The class inheriting from LookupGridBase.
	using ThisType = SharedLookupPartitionedGrid<D, N>;
	/// Lookup value at leaf nodes.
	using LeafType = typename TrackingPartitionedGridTraits<ThisType>::LeafType;
	/// Dimensions of the grid.
	static const UINT Dims = TrackingPartitionedGridTraits<ThisType>::Dims;
	/// Number of tracking lists to use.
	static const UINT NumLists = TrackingPartitionedGridTraits<ThisType>::NumLists;
};

/** @} */ // End group LookupGrid.

} // End namespace felt.
#endif /* PARTITIONEDGRID_HPP_ */
