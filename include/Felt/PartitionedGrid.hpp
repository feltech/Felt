#ifndef PARTITIONEDGRID_HPP_
#define PARTITIONEDGRID_HPP_

#include <memory>
#include <array>
#include <mutex>
#include <boost/iterator/iterator_facade.hpp>

#include "MappedGrid.hpp"

namespace felt
{
/// Default size of a spatial partition (in each dimension).
static const UINT DEFAULT_PARTITION = 4;


/**
 * Base traits class for classes CRTP derived from PartitionBase.
 */
template <class Derived> struct PartitionBaseTraits {};


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

	/// Child object to store in each partition.
	using Child = typename PartitionBaseTraits<Derived>::ChildType;

	/// Number of tracking lists of points.
	static const UINT NUM_LISTS = PartitionBaseTraits<Derived>::NumLists;

	/// Grid of partitions with tracking list(s) of active grid points.
	using BranchGrid = TrackedGrid<
		Child, PartitionBaseTraits<Derived>::Dims, NUM_LISTS
	>;
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
	VecDu	m_udims_child;
	/// A convenience vector of the (signed) size of a partition.
	VecDi 	m_idims_child;

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
		DerivedType* self = static_cast<DerivedType*>(this);

		self->init(partition_size_);
		self->dims(size_);
		self->offset(offset_);
	}

	/**
	 * Initialise size of spatial partitions.
	 *
	 * @param partition_size_ spatial size of the child data structure.
	 */
	void init (const VecDu& partition_size_)
	{
		m_udims_child = partition_size_;
		m_idims_child = partition_size_.template cast<INT>();
	}

	/**
	 * Get size of a spatial partition.
	 */
	const VecDu& child_dims() const
	{
		return m_udims_child;
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
	Child& child (const VecDi& pos)
	{
		return m_grid_branch(pos);
	}

	/**
	 * Get the PartitionBase::Child object at given position.
	 *
	 * Shorthand for branch().get(pos) or branch()(pos).
	 *
	 * @param pos position within the spatial partition grid.
	 * @return PartitionBase::Child at given position.
	 */
	const Child& child (const VecDi& pos) const
	{
		return m_grid_branch(pos);
	}

	/**
	 * Reshape grid, computing the size of the branch grid.
	 *
	 * The branch grid will be increased in size by one, if required, to ensure all leaf nodes are
	 * completely contained.
	 *
	 * @param grid_size_ overall size of grid.
	 */
	void dims (const VecDu& grid_size_)
	{
		VecDu branch_size = (
			grid_size_.array() / m_udims_child.array()
		).matrix();

		if (
			(branch_size.array() * m_udims_child.array()).matrix()
			!= grid_size_
		) {
			branch_size += VecDu::Constant(1);
		}

		m_grid_branch.dims(branch_size);
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
			grid_offset_.array() / m_idims_child.array()
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
 * Base traits class for classes CRTP derived from PartitionedArrayBase.
 */
template <class Derived> struct PartitionedArrayBaseTraits {};


/**
 * Base class for common features of PartitionedArray template
 * specialisations.
 */
template <class Derived>
class PartitionedArrayBase : public PartitionBase<PartitionedArrayBase<Derived> >
{
public:
	using ThisType = PartitionedArrayBase<Derived>;
	using Base = PartitionBase<ThisType>;
	using typename Base::VecDu;
	using typename Base::VecDi;
protected:
	VecDi	m_offset;
public:
	/**
	 * Set offset of 'imaginary' grid containing the list.
	 *
	 * @param offset_grid
	 */
	void offset (const VecDi& offset_grid)
	{
		m_offset = offset_grid;
		Base::offset(offset_grid);
	}

	/**
	 * Get spatial partition from leaf grid node in 'imaginary' grid.
	 *
	 * @param pos_leaf
	 * @return
	 */
	const VecDi pos_child (const VecDi& pos_leaf) const
	{
		return (
			(pos_leaf - m_offset).array() / this->m_idims_child.array()
		).matrix() + this->branch().offset();
	}
};


template <class Derived>
struct PartitionBaseTraits<PartitionedArrayBase<Derived> >
{
	using ThisType = typename PartitionedArrayBaseTraits<Derived>::ThisType;
	using ChildType = typename PartitionedArrayBaseTraits<Derived>::ChildType;
	static const UINT Dims = PartitionedArrayBaseTraits<Derived>::Dims;
	static const UINT NumLists = PartitionedArrayBaseTraits<Derived>::NumLists;
};


/**
 * Spatially partitioned expandable lists.
 *
 * A specialised partitioned grid, where the child grids are simply
 * expandable lists.
 *
 * @tparam T the type to store in elements of the list
 * @tparam D the dimension of the 'imaginary' grid.
 * @tparam N the dimension of the array(s).
 */
template <typename T, UINT D, UINT N=0>
class PartitionedArray
	: public PartitionedArrayBase<PartitionedArray<T, D, N> >
{
	static_assert(N > 0, "Number of arrays N must be greater than 0.");
protected:
	using ThisType = PartitionedArray<T, D, N>;
	using Base = PartitionedArrayBase<ThisType>;

	using typename Base::VecDu;
	using typename Base::VecDi;

public:
	PartitionedArray () = default;

	PartitionedArray (
		const VecDu& dims, const VecDi& offset,
		const VecDu& dims_partition = VecDu::Constant(DEFAULT_PARTITION)
	) : Base()
	{
		this->init(dims, offset, dims_partition);
	}

	/**
	 * Add val to list, placing in partition found from pos.
	 *
	 * @param pos position in 'imaginary' grid.
	 * @param val value to insert in list.
	 */
	void add(const VecDi& pos, const T& val, const UINT& arr_idx = 0)
	{
		const VecDi& pos_child = this->pos_child(pos);
		this->child(pos_child)[arr_idx].push_back(val);
		Base::add_child(pos_child, arr_idx);
	}

	/**
	 * Loop all spatial partitions, resizing their lists to zero.
	 */
	void reset(const UINT& arr_idx = 0)
	{
		for (const VecDi& pos_child : this->branch().list(arr_idx))
			this->child(pos_child)[arr_idx].clear();
		Base::reset(arr_idx);
	}
};


template <typename T, UINT D, UINT N>
struct PartitionedArrayBaseTraits<PartitionedArray<T, D, N> >
{
	using ThisType = PartitionedArray<T, D, N>;
	using ChildType = std::array<std::vector<T, Eigen::aligned_allocator<T> >, N>;
	static const UINT Dims = D;
	static const UINT NumLists = N;
};


/**
 * Spatially partitioned expandable list - 1D array specialisation.
 *
 * A specialised partitioned grid, where the child grids are simply
 * expandable lists.
 *
 * @tparam T the type to store in elements of the list
 * @tparam D the dimension of the 'imaginary' grid.
 */
template <typename T, UINT D>
class PartitionedArray<T, D, 0>
	: public PartitionedArrayBase<PartitionedArray<T, D, 0> >
{
protected:
	using ThisType = PartitionedArray<T, D, 0>;
	using Base = PartitionedArrayBase<ThisType>;

	using typename Base::VecDu;
	using typename Base::VecDi;

public:
	PartitionedArray () = default;

	PartitionedArray (
		const VecDu& dims, const VecDi& offset,
		const VecDu& dims_partition = VecDu::Constant(DEFAULT_PARTITION)
	) : Base()
	{
		this->init(dims, offset, dims_partition);
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
	 * Loop all spatial partitions, resizing their lists to zero.
	 */
	void reset()
	{
		for (const VecDi& pos_child : this->branch().list())
			this->child(pos_child).clear();
		Base::reset();
	}
};


template <typename T, UINT D>
struct PartitionedArrayBaseTraits<PartitionedArray<T, D, 0> >
{
	using ThisType = PartitionedArray<T, D, 0>;
	using ChildType = std::vector<T, Eigen::aligned_allocator<T> >;
	static const UINT Dims = D;
	static const UINT NumLists = 1;
};


template <class Derived> struct PartitionedGridBaseTraits {};


/**
 * A general spatially partitioned grid storing arbitrary values.
 *
 * Inherits directly from child grid class, spoofing the signature, but
 * storage is in a (single level) tree hierarchy, i.e.
 * PartitionedGridBase -> Child -> leaf.
 *
 * @tparam T type of value stored in leaf grid nodes.
 * @tparam D dimension of the grid.
 * @tparam G GridBase subclass used for child grids.
 * @tparam N number of tracking lists to use.
 */
template <class Derived>
class PartitionedGridBase
	: public PartitionedGridBaseTraits<Derived>::MixinType,
	  public PartitionBase<PartitionedGridBase<Derived> >
{
public:
	using DerivedType = Derived;
	using ThisType = PartitionedGridBase<DerivedType>;
	using Base = PartitionBase<ThisType>;
	using MixinType = typename PartitionedGridBaseTraits<Derived>::MixinType;
	using Traits = PartitionedGridBaseTraits<DerivedType>;
	using Child = typename Traits::ChildType;
	using RetType = typename Traits::RetType;
	using LeafType = typename Traits::LeafType;
	static const UINT Dims = Traits::Dims;
	using BranchGrid = TrackedGrid<Child, Dims, Traits::NumLists>;
	using VecDu = typename MixinType::VecDu;
	using VecDi = typename MixinType::VecDi;
	using SnapshotGrid = Grid<LeafType, Dims>;

protected:
	std::unique_ptr<SnapshotGrid>	m_pgrid_snapshot;

public:
	using Base::reset;
	using MixinType::dims;
	using MixinType::offset;

	PartitionedGridBase () = default;

	/**
	 * Constructor
	 * @param dims the dimensions of the whole grid.
	 * @param offset the offset of the whole grid.
	 * @param dims_partition the dimensions of each spatial partition.
	 * @param delta the grid delta value used for spatial derivatives (dx).
	 */
	PartitionedGridBase (
		const VecDu& dims, const VecDi& offset,
		const VecDu& dims_partition = VecDu::Constant(DEFAULT_PARTITION),
		const FLOAT& delta = 1
	) : Base()
	{
		DerivedType* self = static_cast<DerivedType*>(this);
		self->init(dims, offset, dims_partition, delta);
	}

	/**
	 * Initialisation function called by non-trivial constructor and
	 * subclasses.
	 *
	 * @param dims the dimensions of the whole grid.
	 * @param offset the offset of the whole grid.
	 * @param dims_partition the dimensions of each spatial partition.
	 * @param delta the grid delta value used for spatial derivatives (dx).
	 */
	void init (
		const VecDu& dims, const VecDi& offset,
		const VecDu& dims_partition = VecDu::Constant(DEFAULT_PARTITION),
		const FLOAT& delta = 1.0f
	) {
		Base::init(dims_partition);
		MixinType::init(dims, offset, delta);
	}

	/**
	 * Reshape grid.
	 *
	 * @param vec_NewDims
	 * @return
	 */
	void dims (const VecDu& dims_grid)
	{
		Base::dims(dims_grid);

		this->m_dims = dims_grid;

		for (UINT idx = 0; idx < this->branch().data().size(); idx++)
		{
			Child& child = this->branch().data()(idx);
			child.dims(this->m_udims_child);
		}
	}

	/**
	 * Set offset of branch grid and propagate offset to children,
	 * translating as appropriate.
	 *
	 * @param offset_grid
	 */
	void offset (const VecDi& offset_grid)
	{
		Base::offset(offset_grid);
		MixinType::offset(offset_grid);

		for (UINT idx = 0; idx < this->branch().data().size(); idx++)
		{
			Child& child = this->branch().data()(idx);
			const VecDi& pos_child = this->branch().index(idx);
			const VecDi& offset_child = (
				(
					(pos_child - this->branch().offset()).array()
					* this->m_idims_child.array()
				).matrix() + offset_grid
			);

			child.offset(offset_child);
		}
	}

	/**
	 * Calculate the position of a child grid (i.e. partition) given the
	 * position of leaf grid node.
	 *
	 * @param pos_leaf
	 * @return
	 */
	const VecDi pos_child (const VecDi& pos_leaf) const
	{
		return (
			(pos_leaf - offset()).array() / this->m_idims_child.array()
		).matrix() + this->branch().offset();
	}

	/**
	 * Fill with a single value.
	 *
	 * @param val
	 */
	void fill (const LeafType& val)
	{
		for (UINT idx = 0; idx < this->branch().data().size(); idx++)
		{
			Child& child = this->branch().data()(idx);
			child.fill(val);
		}
	}

	/**
	 * Get the leaf grid node at pos by navigating to the correct partition.
	 *
	 * Overrides base Grid class's get method (and thus operator()).
	 *
	 * @param pos
	 * @return value stored at grid point pos.
	 */
	RetType& get (const VecDi& pos)
	{
		const VecDi& pos_child = this->pos_child(pos);
		Child& child = this->branch().get(pos_child);
		return child.get(pos);
	}

	/**
	 * Get the leaf grid node at pos by navigating to the correct partition.
	 *
	 * Overrides base Grid class's get method (and thus operator()).
	 *
	 * @param pos
	 * @return value stored at grid point pos.
	 */
	const RetType& get (const VecDi& pos) const
	{
		const Child& child = this->branch()(pos_child(pos));
		return child.get(pos);
	}

	/**
	 * Get and store a snapshot of the spatially partitioned data in a
	 * contiguous grid.
	 *
	 * Useful for serialization.
	 *
	 * @return
	 */
	typename SnapshotGrid::ArrayData& data()
	{
		m_pgrid_snapshot.reset(new SnapshotGrid(this->dims(), this->offset()));

		for (
			UINT branch_idx = 0; branch_idx < this->branch().data().size();
			branch_idx++
		) {
			const Child& child = this->branch().data()[branch_idx];
			for (
				UINT leaf_idx = 0; leaf_idx < child.data().size();
				leaf_idx++
			) {
				const VecDi& pos = child.index(leaf_idx);
				if (m_pgrid_snapshot->inside(pos))
					m_pgrid_snapshot->get_internal(pos)
						= child.get_internal(pos);
			}
		}
		return m_pgrid_snapshot->data();
	}

	/**
	 * Copy the snapshot of the grid data back into the partitioned
	 * structure.
	 *
	 * Useful for deserialisation.
	 */
	void flush_snapshot()
	{
		if (m_pgrid_snapshot.get() == NULL)
			return;

		for (
			UINT branch_idx = 0; branch_idx < this->branch().data().size();
			branch_idx++
		) {
			Child& child = this->branch().data()[branch_idx];
			for (
				UINT leaf_idx = 0; leaf_idx < child.data().size();
				leaf_idx++
			) {
				const VecDi& pos = child.index(leaf_idx);
				if (m_pgrid_snapshot->inside(pos))
					child.get_internal(pos)
						= m_pgrid_snapshot->get_internal(pos);
			}
		}
	}
};


template <class Derived>
struct PartitionBaseTraits<PartitionedGridBase<Derived> >
{
	using ThisType = typename PartitionedGridBaseTraits<Derived>::ThisType;
	using ChildType = typename PartitionedGridBaseTraits<Derived>::ChildType;
	static const UINT Dims = PartitionedGridBaseTraits<Derived>::Dims;
	static const UINT NumLists = PartitionedGridBaseTraits<Derived>::NumLists;
};


template <typename T, UINT D>
class PartitionedGrid
	: public PartitionedGridBase<PartitionedGrid<T, D> >
{
public:
	using Base = PartitionedGridBase<PartitionedGrid<T, D> >;
	using Base::PartitionedGridBase;
};


template <typename T, UINT D>
struct GridBaseTraits<PartitionedGrid<T, D> >
{
	using ThisType = PartitionedGrid<T, D>;
	static const UINT Dims = D;
	using LeafType = T;
	using RetType = T;
};


template <typename T, UINT D>
struct PartitionedGridBaseTraits<PartitionedGrid<T, D> >
{
	using ThisType = PartitionedGrid<T, D>;
	using LeafType = typename GridBaseTraits<ThisType>::LeafType;
	using RetType = typename GridBaseTraits<ThisType>::RetType;
	using ChildType = Grid<T, D>;
	using MixinType = GridBase<ThisType>;
	static const UINT Dims = D;
	static const UINT NumLists = 1;
};



/**
 * Container wrapping iterator through leafs of partitioned grid tree.
 *
 * @tparam G grid type to iterate over.
 */
template <typename G>
class LeafsContainer
{
private:
	using GridTree = G;
	using VecDi = typename GridTree::VecDi;
	using PosArray = typename GridTree::PosArray;

	const GridTree* m_pgrid;
	const UINT	m_listIdx;

public:
	/**
	 * Iterator class for range-based for loops across partitioned grid
	 * leafs.
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
		 * Constructor.
		 *
		 * @param pgrid grid to iterate through
		 * @param listIdx tracking list id within grid
		 * @param it_child iterator over child grids
		 * @param it_leaf iterator over lowest level leaf grid point.
		 */
		iterator(
			const GridTree* pgrid, const UINT& listIdx,
			const Iter& it_child, const Iter& it_leaf
		)
		: m_pgrid(pgrid), m_listIdx(listIdx), m_it_child(it_child),
		  m_it_leaf(it_leaf),
		  m_it_child_end(pgrid->branch().list(listIdx).end())
		{}

	private:
		friend class boost::iterator_core_access;

		const GridTree* m_pgrid;
		const UINT	m_listIdx;
		Iter		m_it_child;
		Iter		m_it_leaf;
		const Iter	m_it_child_end;

		/**
		 * Override specifying how to move on to next leaf.
		 */
		void increment()
		{
			if (m_it_child == m_it_child_end)
				return;

			if (m_it_leaf == m_it_child)
			{
				m_it_leaf = m_pgrid->child(
					*m_it_child
				).list(m_listIdx).begin();
				return;
			}

			m_it_leaf++;

			if (
				m_it_leaf == m_pgrid->child(
					*m_it_child
				).list(m_listIdx).end()
			) {
				m_it_child++;
				m_it_leaf = m_it_child;
				increment();
				return;
			}
		}

		/**
		 * Check for equality between this iterator and another.
		 * @param other
		 * @return
		 */
		bool equal(const iterator& other) const
		{
			return (
				m_it_child == other.m_it_child
				&& m_it_leaf == other.m_it_leaf
			);
		}

		/**
		 * Dereference iterator into grid position.
		 */
		const VecDi& dereference() const {
			const VecDi& pos = *m_it_leaf;
			return *m_it_leaf;
		}
	};

	/**
	 * Constructor.
	 *
	 * @param pgrid grid to iterate over
	 * @param listIdx tracking list id identifying leafs
	 */
	LeafsContainer(const GridTree* pgrid, const UINT& listIdx)
	: m_pgrid(pgrid), m_listIdx(listIdx)
	{}

	/**
	 * Get first iterator for leafs identified within list.
	 *
	 * @return
	 */
	const iterator begin() const
	{
		const typename PosArray::const_iterator& it_child_begin
			= m_pgrid->branch().list(m_listIdx).cbegin();
		const typename PosArray::const_iterator& it_child_end
			= m_pgrid->branch().list(m_listIdx).cend();

		typename PosArray::const_iterator it_leaf_begin = it_child_end;

		if (it_child_begin != it_child_end)
		{
			it_leaf_begin = m_pgrid->child(
				*it_child_begin
			).list(m_listIdx).begin();
		}

		return iterator(
			m_pgrid, m_listIdx,
			it_child_begin,
			it_leaf_begin
		);
	}

	/**
	 * Iterator representing one element past the end of the list of leafs.
	 *
	 * @return
	 */
	const iterator end() const
	{
		return iterator(
			m_pgrid, m_listIdx,
			m_pgrid->branch().list(m_listIdx).cend(),
			m_pgrid->branch().list(m_listIdx).cend()
		);
	}

	/**
	 * Calculate length of list by summing lists in all partitions.
	 *
	 * @return
	 */
	UINT size() const
	{
		UINT sum = 0;
		for (const VecDi& pos_child : m_pgrid->branch().list(m_listIdx))
			sum += m_pgrid->child(pos_child).list(m_listIdx).size();
		return sum;
	}
};

/**
 * Base traits for spatially partitioned wrapper for lookup and tracked
 * grid.
 */
template <class Derived> struct TrackingPartitionedGridTraits {};

/**
 * Base class for spatially partitioned wrapper for lookup and tracked grid.
 *
 * @tparam T type of value stored in leaf grid nodes.
 * @tparam D dimension of the grid.
 * @tparam N number of tracking lists to use.
 * @tparam G lookup or tracked grid class used for child grids.
 */
template <class Derived>
class TrackingPartitionedGridBase
	: public PartitionedGridBase<TrackingPartitionedGridBase<Derived> >
{
public:
	using ThisType = TrackingPartitionedGridBase<Derived>;
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
	 * @param arr_idx tracking list id.
	 */
	void reset(const UINT& arr_idx = 0)
	{
		BranchGrid& branch = this->branch();
		for (const VecDi& pos_child : branch.list(arr_idx))
			branch(pos_child).reset(arr_idx);

		Base::reset(arr_idx);
	}

	/**
	 * Add a leaf position to be tracked to given tracking list.
	 *
	 * Descend to relevant child grid to add to their tracking structure.
	 *
	 * @param pos position to add.
	 * @param arr_idx tracking list id.
	 * @return true if grid node set in child lookup grid and position added
	 * to tracking list, false if child grid node was already set so
	 * position already in a list.
	 */
	bool add(const VecDi& pos, const UINT& arr_idx = 0)
	{
		const VecDi& pos_child = this->pos_child(pos);
		Base::add_child(pos_child, arr_idx);
		return this->child(pos_child).add(pos, arr_idx);
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
	 * @param pos position to add.
	 * @param arr_idx tracking list id.
	 * @return true if grid node set in child lookup grid and position added
	 * to tracking list, false if child grid node was already set so
	 * position already in a list.
	 */
	bool add_safe(const VecDi& pos, const UINT& arr_idx = 0)
	{
		const VecDi& pos_child = this->pos_child(pos);
		Base::add_child(pos_child, arr_idx);
		Child& child = this->child(pos_child);
		if (child.is_active(pos))
			return false;
		std::lock_guard<std::mutex> lock(child.mutex());
		return child.add(pos, arr_idx);
	}

	/**
	 * Remove a leaf position from relevant child tracking structure and
	 * remove child from tracking list if child's list is now empty.
	 *
	 * @param pos leaf position to remove.
	 * @param arr_idx tracking list id.
	 */
	void remove(const VecDi& pos, const UINT& arr_idx = 0)
	{
		const VecDi& pos_child = this->pos_child(pos);
		Child& child = this->child(pos_child);
		child.remove(pos, arr_idx);
		if (child.list(arr_idx).size() == 0)
			this->remove_child(pos_child, arr_idx);
	}

	/**
	 * Return structure for range based for loops over leaf nodes.
	 *
	 * @param listIdx tracking list id.
	 * @return
	 */
	const LeafsContainer<Derived> leafs(const UINT& listIdx = 0) const
	{
		return LeafsContainer<Derived>(
			static_cast<const Derived*>(this), listIdx
		);
	}

private:
	/**
	 * Override spoofed (non-partitioned) base class's list method to make
	 * private, since it will always be empty (must use branch or child).
	 *
	 * @param arr_idx
	 * @return
	 */
	PosArray& list(const UINT& arr_idx = 0)
	{
		return Child::list(arr_idx);
	}
};


template <class Derived>
struct PartitionedGridBaseTraits<TrackingPartitionedGridBase<Derived> >
{
	using ThisType = typename TrackingPartitionedGridTraits<Derived>::ThisType;
	using LeafType = typename TrackingPartitionedGridTraits<ThisType>::LeafType;
	using RetType = typename TrackingPartitionedGridTraits<ThisType>::RetType;
	using ChildType = typename TrackingPartitionedGridTraits<ThisType>::ChildType;
	using MixinType = typename TrackingPartitionedGridTraits<ThisType>::MixinType;
	static const UINT Dims = TrackingPartitionedGridTraits<ThisType>::Dims;
	static const UINT NumLists = TrackingPartitionedGridTraits<ThisType>::NumLists;
};


/**
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
	using ThisType = TrackedPartitionedGrid<T, D, N>;
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
	 * @param pos position in grid.
	 * @param val value to set.
	 * @param arr_idx tracking list id.
	 * @return true if grid node set in child lookup grid and position added
	 * to tracking list, false if child grid node was already set so
	 * position already in a list.
	 */
	bool add(const VecDi& pos, const T& val, const UINT& arr_idx = 0)
	{
		const VecDi& pos_child = this->pos_child(pos);
		Base::add_child(pos_child, arr_idx);
		return this->child(pos_child).add(pos, val, arr_idx);
	}

	/**
	 * Thread-safely set value in grid at given position and add position to lookup grid.
	 *
	 * Will set value regardless whether lookup grid already set for given
	 * position + tracking list.
	 *
	 * Descend to relevant child grid to add to their tracking structure.
	 *
	 * @param pos position in grid.
	 * @param val value to set.
	 * @param arr_idx tracking list id.
	 * @return true if grid node set in child lookup grid and position added
	 * to tracking list, false if child grid node was already set so
	 * position already in a list.
	 */
	bool add_safe(const VecDi& pos, const T& val, const UINT& arr_idx = 0)
	{
		const VecDi& pos_child = this->pos_child(pos);
		Base::add_child(pos_child, arr_idx);
		Child& child = this->child(pos_child);
		std::lock_guard<std::mutex> lock(child.mutex());
		return child.add(pos, val, arr_idx);
	}

	/**
	 * Set every active grid node (i.e. those referenced by lookup grid)
	 * to given value and reset the lookup grid.
	 *
	 * Lookup grid will then be full of NULL indices for given tracking list
	 * and the relevant tracking list(s) will be empty.
	 *
	 * @param val value to set in main grid.
	 * @param arr_idx tracking list id to cycle over and clear.
	 */
	void reset(const T& val, const UINT& arr_idx = 0)
	{
		BranchGrid& branch = this->branch();
		for (const VecDi& pos_child : branch.list(arr_idx))
			branch(pos_child).reset(val, arr_idx);

		Base::reset(arr_idx);
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


template <typename T, UINT D, UINT N>
struct TrackedGridBaseTraits<TrackedPartitionedGrid<T, D, N> >
{
	using ThisType = TrackedPartitionedGrid<T, D, N>;
	using LeafType = T;
	using LookupType = LookupGrid<D, N>;
	static const UINT Dims = D;
};

/**
 * Traits class for spatially partitioned wrapper for TrackedGrid.
 *
 * @tparam T type of value stored in leaf grid nodes.
 * @tparam D dimension of the grid.
 * @tparam N number of tracking lists to use.
 */
template <typename T, UINT D, UINT N>
struct TrackingPartitionedGridTraits<TrackedPartitionedGrid<T, D, N> >
{
	using ThisType = TrackedPartitionedGrid<T, D, N>;
	using ChildType = TrackedGrid<T, D, N>;
	using MixinType = TrackedGridBase<ThisType>;
	using LeafType = typename TrackedGridBaseTraits<ThisType>::LeafType;
	using RetType = LeafType;
	static const UINT Dims = TrackedGridBaseTraits<ThisType>::Dims;
	static const UINT NumLists = N;
};

/**
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
	 * @param pos position in grid.
	 * @param val value to set.
	 * @param arr_idx tracking list id.
	 * @return true if grid node set in child lookup grid and position added
	 * to tracking list, false if child grid node was already set so
	 * position already in a list.
	 */
	bool add(const VecDi& pos, const T& val, const UINT& arr_idx)
	{
		const VecDi& pos_child = this->pos_child(pos);
		Base::add_child(pos_child, arr_idx);
		return this->child(pos_child).add(pos, val, arr_idx);
	}


	/**
	 * Add a position to the lookup grid.
	 *
	 * @param pos position in the grid to add.
	 * @param arr_idx tracking list id.
	 * @return true if grid node set in child lookup grid and position added
	 * to tracking list, false if leaf grid node was already set so position
	 * already in a list.
	 */
	bool add(const VecDi& pos, const UINT& arr_idx)
	{
		return Base::add(pos, arr_idx);
	}

	/**
	 * Set every active grid node (i.e. those referenced by lookup grid)
	 * to given value and reset the lookup grid.
	 *
	 * Lookup grid will then be full of NULL indices for given tracking list
	 * and the relevant tracking list(s) will be empty.
	 *
	 * @param val value to set in main grid.
	 * @param arr_idx tracking list id to cycle over and clear.
	 */
	void reset(const T& val, const UINT& arr_idx = 0)
	{
		BranchGrid& branch = this->branch();
		for (const VecDi& pos_child : branch.list(arr_idx))
			branch(pos_child).reset(val, arr_idx);

		Base::reset(arr_idx);
	}
};


template <typename T, UINT D, UINT N>
struct TrackedGridBaseTraits<SharedTrackedPartitionedGrid<T, D, N> >
{
	using ThisType = SharedTrackedPartitionedGrid<T, D, N>;
	using LeafType = T;
	using LookupType = SharedLookupGrid<D, N>;
	static const UINT Dims = D;
};


/**
 * Traits class for spatially partitioned wrapper for SharedTrackedGrid.
 *
 * @tparam T type of value stored in leaf grid nodes.
 * @tparam D dimension of the grid.
 * @tparam N number of tracking lists to use.
 */
template <typename T, UINT D, UINT N>
struct TrackingPartitionedGridTraits<SharedTrackedPartitionedGrid<T, D, N> >
{
	using ThisType = SharedTrackedPartitionedGrid<T, D, N>;
	using ChildType = SharedTrackedGrid<T, D, N>;
	using MixinType = TrackedGridBase<ThisType>;
	using LeafType = typename TrackedGridBaseTraits<ThisType>::LeafType;
	using RetType = LeafType;
	static const UINT Dims = TrackedGridBaseTraits<ThisType>::Dims;
	static const UINT NumLists = N;
};


/**
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
 * Traits class for spatially partitioned wrapper for LookupGrid.
 *
 * @tparam D dimension of the grid.
 * @tparam N number of tracking lists to use.
 */
template <UINT D, UINT N>
struct TrackingPartitionedGridTraits<LookupPartitionedGrid<D, N> >
{
	using ThisType = LookupPartitionedGrid<D, N>;
	using LeafType = Eigen::Matrix<UINT, N, 1>;
	using RetType = LeafType;
	using ChildType = LookupGrid<D, N>;
	using MixinType = LookupGridBase<ThisType>;
	static const UINT Dims = D;
	static const UINT NumLists = N;
};


template <UINT D, UINT N>
struct LookupGridBaseTraits<LookupPartitionedGrid<D, N> >
{
	using ThisType = LookupPartitionedGrid<D, N>;
	using LeafType = typename TrackingPartitionedGridTraits<ThisType>::LeafType;
	using RetType = typename TrackingPartitionedGridTraits<ThisType>::RetType;
	static const UINT Dims = TrackingPartitionedGridTraits<ThisType>::Dims;
	static const UINT NumLists = TrackingPartitionedGridTraits<ThisType>::NumLists;
};


/**
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
 * Traits class for spatially partitioned wrapper for SharedLookupGrid.
 *
 * @tparam D dimension of the grid.
 * @tparam N number of tracking lists to use.
 */
template <UINT D, UINT N>
struct TrackingPartitionedGridTraits<SharedLookupPartitionedGrid<D, N> >
{
	using ThisType = SharedLookupPartitionedGrid<D, N>;
	using LeafType = Eigen::Matrix<UINT, N, 1>;
	using RetType = UINT;
	using ChildType = SharedLookupGrid<D, N>;
	using MixinType = SharedLookupGridBase<ThisType>;
	static const UINT Dims = D;
	static const UINT NumLists = N;
};


template <UINT D, UINT N>
struct SharedLookupGridBaseTraits<SharedLookupPartitionedGrid<D, N> >
{
	using ThisType = SharedLookupPartitionedGrid<D, N>;
	using LeafType = typename TrackingPartitionedGridTraits<ThisType>::LeafType;
	using RetType = typename TrackingPartitionedGridTraits<ThisType>::RetType;
	static const UINT Dims = TrackingPartitionedGridTraits<ThisType>::Dims;
	static const UINT NumLists = TrackingPartitionedGridTraits<ThisType>::NumLists;
};
}
#endif /* PARTITIONEDGRID_HPP_ */
