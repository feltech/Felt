#ifndef PARTITIONEDGRID_HPP_
#define PARTITIONEDGRID_HPP_

#include <boost/iterator/iterator_facade.hpp>

#include "MappedGrid.hpp"

namespace felt
{
	#define DEFAULT_PARTITION 4

	/**
	 * Base class for spatially partitioned structures.
	 *
	 * @tparam D number of dimensions of space.
	 * @tparam G type of object to store in a partition.
	 * @tparam N number of tracking lists to use (see TrackedGrid).
	 */
	template <UINT D, class G, UINT N=1>
	class PartitionBase
	{
	public:
		/// Child object to store in each partition.
		typedef G							Child;
		/// Grid of partitions with tracking list(s) of active grid points.
		typedef TrackedGrid<Child, D, N>	BranchGrid;
	protected:
		/// D-dimensional unsigned int vector.
		typedef typename BranchGrid::VecDu	VecDu;
		/// D-dimensional signed int vector.
		typedef typename BranchGrid::VecDi	VecDi;

		/// Grid of partitions with tracking list(s) of active grid points.
		BranchGrid	m_grid_branch;

		/// A convenience vector of the (unsigned) size of a partition.
		VecDu	m_udims_child;
		/// A convenience vector of the (signed) size of a partition.
		VecDi 	m_idims_child;

	public:
		virtual ~PartitionBase ()
		{}

		PartitionBase () : m_grid_branch()
		{}


		PartitionBase (
			const VecDu& dims, const VecDi& offset, const VecDu& dims_partition
		) : m_grid_branch()
		{
			this->init(dims, offset, dims_partition);
		}


		/**
		 * Initialisation method to be called by non-trivial constructor or
		 * subclasses.
		 *
		 * Similar to a Grid::init(), with the addition of setting the size of
		 * spatial partitions.
		 *
		 * @param dims
		 * @param offset
		 * @param dims_partition
		 */
		void init (
			const VecDu& dims, const VecDi& offset,
			const VecDu& dims_partition = VecDu::Constant(DEFAULT_PARTITION)
		) {
			this->init(dims_partition);
			this->dims(dims);
			this->offset(offset);
		}

		/**
		 * Initialise size of spatial partitions.
		 *
		 * @param dims_partition
		 */
		void init (const VecDu& dims_partition)
		{
			m_udims_child = dims_partition;
			m_idims_child = dims_partition.template cast<INT>();
		}

		/**
		 * Get TrackedGrid branch grid - the spatial partition grid that stores
		 * the child objects.
		 * @see PartitionedGrid
		 * @see PartitionedArray
		 *
		 * @return
		 */
		BranchGrid& branch ()
		{
			return m_grid_branch;
		}

		/**
		 * Get TrackedGrid branch grid - the spatial partition grid that stores
		 * the child objects.
		 * @see PartitionedGrid
		 * @see PartitionedArray
		 *
		 * @return
		 */
		const BranchGrid& branch () const
		{
			return m_grid_branch;
		}
		

		/**
		 * Get the child object at given position.
		 *
		 * Shorthand for branch().get(pos) or branch()(pos).
		 *
		 * @param pos
		 * @return
		 */
		Child& child (const VecDi& pos)
		{
			return m_grid_branch(pos);
		}

		/**
		 * Get the child object at given position.
		 *
		 * Shorthand for branch().get(pos) or branch()(pos).
		 *
		 * @param pos
		 * @return
		 */
		const Child& child (const VecDi& pos) const
		{
			return m_grid_branch(pos);
		}

		/**
		 * Reshape grid, computing the size of the branch grid.
		 *
		 * The branch grid will be increased in size by one, if required, to
		 * ensure dims_grid leaf nodes are completely contained.
		 *
		 * @param dims_grid
		 * @return
		 */
		virtual void dims (const VecDu& dims_grid)
		{
			VecDu dims_branch = (
				dims_grid.array() / m_udims_child.array()
			).matrix();

			if (
				(dims_branch.array() * m_udims_child.array()).matrix()
				!= dims_grid
			) {
				dims_branch += VecDu::Constant(1);
			}

			m_grid_branch.dims(dims_branch);
		}


		/**
		 * Calculate the offset of the branch grid from the given dims and the
		 * size of a spatial partition (the child size).
		 *
		 * @param offset_grid
		 */
		virtual void offset (const VecDi& offset_grid)
		{
			const VecDi& offset_parent = (
				offset_grid.array() / m_idims_child.array()
			).matrix();

			m_grid_branch.offset(offset_parent);
		}

		/**
		 * Remove a child at given position and tracking list index from
		 * branch grid's tracking subgrid.
		 *
		 * @param pos
		 * @param arr_idx
		 * @return
		 */
		bool add_child(const VecDi& pos, const UINT& arr_idx = 0)
		{
			return this->branch().add(pos, arr_idx);
		}

		/**
		 * Remove a child at given position and tracking list index from
		 * branch grid's tracking subgrid.
		 *
		 * @param pos
		 * @param arr_idx
		 * @return
		 */
		void remove_child(const VecDi& pos, const UINT& arr_idx = 0)
		{
			this->branch().remove(pos, arr_idx);
		}

		/**
		 * Reset the tracking list at given index in the branch grid.
		 *
		 * @param arr_idx
		 */
		void reset(const UINT& arr_idx = 0)
		{
			this->branch().reset(arr_idx);
		}
	};

	template <typename T, UINT D, UINT N, typename A>
	class PartitionedArrayBase : public PartitionBase<D, A, N>
	{
	protected:
		typedef PartitionBase<D, A, N> Base;
		typedef typename Base::VecDu	VecDu;
		typedef typename Base::VecDi	VecDi;

		VecDi	m_offset;

	public:
		PartitionedArrayBase () : Base()
		{}

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

		void dims (const VecDu& dims_grid)
		{
			Base::dims(dims_grid);
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
	class PartitionedArray : public PartitionedArrayBase<
		T, D, N, std::array<std::vector<T, Eigen::aligned_allocator<T> >, N>
	>
	{
		static_assert(N > 0, "Number of arrays N must be greater than 0.");
	protected:
		typedef PartitionedArrayBase<
			T, D, N, std::array<std::vector<T, Eigen::aligned_allocator<T> >, N>
		> Base;

		typedef typename Base::VecDu	VecDu;
		typedef typename Base::VecDi	VecDi;

	public:
		PartitionedArray () : Base()
		{}


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
	class PartitionedArray<T, D, 0> : public PartitionedArrayBase<
		T, D, 1, std::vector<T, Eigen::aligned_allocator<T> >
	>
	{
	protected:
		typedef PartitionedArrayBase<
			T, D, 1, std::vector<T, Eigen::aligned_allocator<T> >
		> Base;

		typedef typename Base::VecDu	VecDu;
		typedef typename Base::VecDi	VecDi;

	public:
		PartitionedArray () : Base()
		{}

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


	/**
	 * A general spatially partitioned grid storing arbitrary values.
	 *
	 * Inherits directly from child grid G class, spoofing the signature.
	 *
	 * @tparam T type of value stored in leaf grid nodes.
	 * @tparam D dimension of the grid.
	 * @tparam G GridBase subclass used for child grids.
	 * @tparam N number of tracking lists to use.
	 */
	template <
		typename T, UINT D, class G=Grid<T, D>, UINT N=1
	>
	class PartitionedGrid : public G, public PartitionBase<D, G, N>
	{
	public:
		typedef PartitionBase<D, G, N>			Base;
		typedef typename Base::Child			ChildGrid;
		typedef TrackedGrid<ChildGrid, D, N>	BranchGrid;
	protected:
		typedef typename ChildGrid::VecDu	VecDu;
		typedef typename ChildGrid::VecDi	VecDi;

		Grid<T, D>*	m_pgrid_snapshot;
	public:
		virtual ~PartitionedGrid ()
		{
			delete m_pgrid_snapshot;
		}

		PartitionedGrid () : Base(), ChildGrid(), m_pgrid_snapshot(NULL)
		{}

		/**
		 * Constructor
		 * @param dims the dimensions of the whole grid.
		 * @param offset the offset of the whole grid.
		 * @param dims_partition the dimensions of each spatial partition.
		 * @param delta the grid delta value used for spatial derivatives (dx).
		 */
		PartitionedGrid (
			const VecDu& dims, const VecDi& offset,
			const VecDu& dims_partition = VecDu::Constant(2),
			const FLOAT& delta = 1
		) : Base(), ChildGrid(), m_pgrid_snapshot(NULL)
		{
			this->init(dims, offset, dims_partition, delta);
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
			ChildGrid::init(dims, offset, delta);
		}


		/**
		 * Get grid dimensions.
		 *
		 * @return
		 */
		const VecDu& dims () const
		{
			return ChildGrid::dims();
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
				ChildGrid& child = this->branch().data()(idx);
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
			ChildGrid::offset(offset_grid);

			for (UINT idx = 0; idx < this->branch().data().size(); idx++)
			{
				ChildGrid& child = this->branch().data()(idx);
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
		 * Get the offset of the whole grid.
		 *
		 * @return
		 */
		const VecDi offset () const
		{
			return ChildGrid::offset();
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
		void fill (const T& val)
		{
			for (UINT idx = 0; idx < this->branch().data().size(); idx++)
			{
				ChildGrid& child = this->branch().data()(idx);
				child.fill(val);
			}
		}

		/**
		 * Get the leaf grid node at pos by navigating to the correct partition.
		 *
		 * Overrides base G grid class's get method (and thus operator()).
		 *
		 * @param pos
		 * @return value stored at grid point pos.
		 */
		T& get (const VecDi& pos)
		{
			ChildGrid& child = this->branch()(pos_child(pos));
			return child.get(pos);
		}

		const T& get (const VecDi& pos) const
		{
			const ChildGrid& leaf = this->branch()(pos_child(pos));
			return leaf.get(pos);
		}

		/**
		 * Delegate to PartitionBase::reset.
		 *
		 * Tracking lists are reset, but data in the grid remains unchanged.
		 *
		 * @param arr_idx
		 */
		void reset(const UINT& arr_idx = 0)
		{
			Base::reset(arr_idx);
		}

		/**
		 * Get and store a snapshot of the spatially partitioned data in a
		 * contiguous grid.
		 *
		 * Useful for serialization.
		 *
		 * @return
		 */
		typename Grid<T, D>::ArrayData& data()
		{
			delete m_pgrid_snapshot;
			m_pgrid_snapshot = new Grid<T, D>(this->dims(), this->offset());

			for (
				UINT branch_idx = 0; branch_idx < this->branch().data().size();
				branch_idx++
			) {
				const ChildGrid& child = this->branch().data()[branch_idx];
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
			if (m_pgrid_snapshot == NULL)
				return;

			for (
				UINT branch_idx = 0; branch_idx < this->branch().data().size();
				branch_idx++
			) {
				ChildGrid& child = this->branch().data()[branch_idx];
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


	/**
	 * Container wrapping iterator through leafs of partitioned grid tree.
	 *
	 * @tparam G grid type to iterate over.
	 */
	template <typename G>
	class LeafsContainer
	{
	private:
		typedef G							GridTree;
		typedef typename GridTree::VecDi	VecDi;
		typedef typename GridTree::PosArray	PosArray;

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
			typedef typename PosArray::const_iterator	Iter;
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


	/*
	 * Spatially partitioned version of MappedGrid.
	 *
	 * @tparam T type of value stored in leaf grid nodes.
	 * @tparam D dimension of the grid.
	 * @tparam N number of tracking lists to use.
	 */
//	template <typename T, UINT D, UINT N>
//	class MappedPartitionedGrid : public PartitionedGrid<
//		T, D, MappedGrid<T, D, N>, N
//	>
//	{
//	protected:
//		typedef PartitionedGrid<
//			T, D, MappedGrid<T, D, N>, N
//		> Base;
//		typedef MappedPartitionedGrid<T, D, N>			ThisType;
//		friend class LeafsContainer<ThisType>;
//	public:
//		typedef typename Base::ChildGrid				ChildGrid;
//	protected:
//		typedef typename ChildGrid::VecDu				VecDu;
//		typedef typename ChildGrid::VecDi				VecDi;
//	public:
//		typedef typename Base::BranchGrid				BranchGrid;
//		typedef typename ChildGrid::PosArray			PosArray;
//
//	public:
//		MappedPartitionedGrid () : Base()
//		{}
//
//
//		MappedPartitionedGrid (
//			const VecDu& dims, const VecDi& offset,
//			const VecDu& dims_partition = VecDu::Constant(DEFAULT_PARTITION),
//			const FLOAT& delta = 1
//		) : Base(dims, offset, dims_partition, delta)
//		{}
//
//
//		bool add(const VecDi& pos, const T& val, const UINT& arr_idx = 0)
//		{
//			BranchGrid& branch = this->branch();
//			const VecDi& pos_child = this->pos_child(pos);
//			Base::add_child(pos_child, arr_idx);
//			return branch(pos_child).add(pos, val, arr_idx);
//		}
//
//		const LeafsContainer<ThisType> leafs(const UINT& listIdx) const
//		{
//			return LeafsContainer<ThisType>(this, listIdx);
//		}
//
//		void reset(const T& val, const UINT& arr_idx = 0)
//		{
//			BranchGrid& branch = this->branch();
//			for (const VecDi& pos_child : branch.list(arr_idx))
//				branch(pos_child).reset(val, arr_idx);
//
//			Base::reset(arr_idx);
//		}
//	};


	/**
	 * Base class for spatially partitioned wrapper for lookup and tracked grid.
	 *
	 * @tparam T type of value stored in leaf grid nodes.
	 * @tparam D dimension of the grid.
	 * @tparam N number of tracking lists to use.
	 * @tparam G lookup or tracked grid class used for child grids.
	 */
	template <typename T, UINT D, UINT N, class G>
	class TrackingPartitionedGridBase : public PartitionedGrid<T, D, G, N>
	{
	protected:
		typedef PartitionedGrid<T, D, G, N> 				Base;
		typedef TrackingPartitionedGridBase<T, D, N, G>		ThisType;
		friend class LeafsContainer<ThisType>;
	public:
		typedef typename Base::ChildGrid				ChildGrid;
	protected:
		typedef typename ChildGrid::VecDu				VecDu;
		typedef typename ChildGrid::VecDi				VecDi;
	public:
		typedef typename Base::BranchGrid				BranchGrid;
		typedef typename ChildGrid::PosArray			PosArray;

	public:
		virtual ~TrackingPartitionedGridBase()
		{}

		TrackingPartitionedGridBase () : Base()
		{}

		/**
		 * Constructor
		 * @param dims the dimensions of the whole grid.
		 * @param offset the offset of the whole grid.
		 * @param dims_partition the dimensions of each spatial partition.
		 * @param delta the grid delta value used for spatial derivatives (dx).
		 */
		TrackingPartitionedGridBase (
			const VecDu& dims, const VecDi& offset,
			const VecDu& dims_partition = VecDu::Constant(DEFAULT_PARTITION),
			const FLOAT& delta = 1
		) : Base(dims, offset, dims_partition, delta)
		{}


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
		 * Remove a leaf position from relevant child tracking structure and
		 * remove child from tracking list if child's list is now empty.
		 *
		 * @param pos leaf position to remove.
		 * @param arr_idx tracking list id.
		 */
		void remove(const VecDi& pos, const UINT& arr_idx = 0)
		{
			const VecDi& pos_child = this->pos_child(pos);
			ChildGrid& child = this->child(pos_child);
			child.remove(pos, arr_idx);
			if (child.list(arr_idx).size() == 0)
				this->remove_child(pos_child, arr_idx);
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
			return ChildGrid::list(arr_idx);
		}
	};


	/**
	 * Spatially partitioned wrapper for TrackedGrid.
	 *
	 * @tparam T type of value stored in leaf grid nodes.
	 * @tparam D dimension of the grid.
	 * @tparam N number of tracking lists to use.
	 */
	template <typename T, UINT D, UINT N>
	class TrackedPartitionedGrid : public TrackingPartitionedGridBase<
		T, D, N, TrackedGrid<T, D, N>
	>
	{
	protected:
		typedef TrackingPartitionedGridBase<
			T, D, N, TrackedGrid<T, D, N>
		> Base;
		typedef TrackedPartitionedGrid<T, D, N>			ThisType;
		friend class LeafsContainer<ThisType>;
	public:
		typedef typename Base::ChildGrid				ChildGrid;
	protected:
		typedef typename ChildGrid::VecDu				VecDu;
		typedef typename ChildGrid::VecDi				VecDi;
	public:
		typedef typename Base::BranchGrid				BranchGrid;
		typedef typename ChildGrid::PosArray			PosArray;

	public:
		TrackedPartitionedGrid () : Base()
		{}

		/**
		 * Constructor
		 * @param dims the dimensions of the whole grid.
		 * @param offset the offset of the whole grid.
		 * @param dims_partition the dimensions of each spatial partition.
		 * @param delta the grid delta value used for spatial derivatives (dx).
		 */
		TrackedPartitionedGrid (
			const VecDu& dims, const VecDi& offset,
			const VecDu& dims_partition = VecDu::Constant(DEFAULT_PARTITION),
			const FLOAT& delta = 1
		) : Base(dims, offset, dims_partition, delta)
		{}

		/**
		 * Return structure for range based for loops over leaf nodes.
		 *
		 * @param listIdx tracking list id.
		 * @return
		 */
		const LeafsContainer<ThisType> leafs(const UINT& listIdx) const
		{
			return LeafsContainer<ThisType>(this, listIdx);
		}

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


	/**
	 * Spatially partitioned wrapper for SharedTrackedGrid.
	 *
	 * @tparam T type of value stored in leaf grid nodes.
	 * @tparam D dimension of the grid.
	 * @tparam N number of tracking lists to use.
	 */
	template <typename T, UINT D, UINT N>
	class SharedTrackedPartitionedGrid : public TrackingPartitionedGridBase<
		T, D, N, SharedTrackedGrid<T, D, N>
	>
	{
	protected:
		typedef TrackingPartitionedGridBase<
			T, D, N, SharedTrackedGrid<T, D, N>
		> Base;
		typedef SharedTrackedPartitionedGrid<T, D, N>	ThisType;
		friend class LeafsContainer<ThisType>;
	public:
		typedef typename Base::ChildGrid				ChildGrid;
		typedef typename ChildGrid::VecDu				VecDu;
		typedef typename ChildGrid::VecDi				VecDi;
		typedef typename Base::BranchGrid				BranchGrid;
		typedef typename ChildGrid::PosArray			PosArray;


		SharedTrackedPartitionedGrid () : Base()
		{}

		/**
		 * Constructor
		 * @param dims the dimensions of the whole grid.
		 * @param offset the offset of the whole grid.
		 * @param dims_partition the dimensions of each spatial partition.
		 * @param delta the grid delta value used for spatial derivatives (dx).
		 */
		SharedTrackedPartitionedGrid (
			const VecDu& dims, const VecDi& offset,
			const VecDu& dims_partition = VecDu::Constant(DEFAULT_PARTITION),
			const FLOAT& delta = 1
		) : Base(dims, offset, dims_partition, delta)
		{}

		/**
		 * Return structure for range based for loops over leaf nodes.
		 *
		 * @param listIdx tracking list id.
		 * @return
		 */
		const LeafsContainer<ThisType> leafs(const UINT& listIdx) const
		{
			return LeafsContainer<ThisType>(this, listIdx);
		}

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


	/**
	 * Spatially partitioned wrapper for LookupGrid.
	 *
	 * @tparam D dimension of the grid.
	 * @tparam N number of tracking lists to use.
	 */
	template <UINT D, UINT N=1>
	class LookupPartitionedGrid : public TrackingPartitionedGridBase<
		Eigen::Matrix<UINT, N, 1>, D, N, LookupGrid<D, N>
	>
	{
	protected:
		typedef LookupPartitionedGrid<D, N>				ThisType;
		typedef TrackingPartitionedGridBase<
			Eigen::Matrix<UINT, N, 1>, D, N, LookupGrid<D, N>
		> Base;
	public:
		typedef typename Base::ChildGrid				ChildGrid;
		typedef typename ChildGrid::VecDu				VecDu;
		typedef typename ChildGrid::VecDi				VecDi;
		typedef typename Base::BranchGrid				BranchGrid;

	public:
		LookupPartitionedGrid () : Base()
		{}


		LookupPartitionedGrid (
			const VecDu& dims, const VecDi& offset,
			const VecDu& dims_partition = VecDu::Constant(DEFAULT_PARTITION),
			const FLOAT& delta = 1
		) : Base(dims, offset, dims_partition, delta)
		{}

		/**
		 * Return structure for range based for loops over leaf nodes.
		 *
		 * @param listIdx tracking list id.
		 * @return
		 */
		const LeafsContainer<ThisType> leafs(const UINT& listIdx = 0) const
		{
			return LeafsContainer<ThisType>(this, listIdx);
		}
	};

	/**
	 * Spatially partitioned wrapper for SharedLookupGrid.
	 *
	 * @tparam D dimension of the grid.
	 * @tparam N number of tracking lists to use.
	 */
	template <UINT D, UINT N>
	class SharedLookupPartitionedGrid : public TrackingPartitionedGridBase<
		UINT, D, N, SharedLookupGrid<D, N>
	>
	{
	protected:
		typedef SharedLookupPartitionedGrid<D, N>		ThisType;
		typedef TrackingPartitionedGridBase<
			UINT, D, N, SharedLookupGrid<D, N>
		> Base;
		typedef typename Base::ChildGrid				ChildGrid;
		typedef typename ChildGrid::VecDu				VecDu;
		typedef typename ChildGrid::VecDi				VecDi;
	public:
		typedef typename Base::BranchGrid				BranchGrid;

	public:
		SharedLookupPartitionedGrid () : Base()
		{}


		SharedLookupPartitionedGrid (
			const VecDu& dims, const VecDi& offset,
			const VecDu& dims_partition = VecDu::Constant(DEFAULT_PARTITION),
			const FLOAT& delta = 1
		) : Base(dims, offset, dims_partition, delta)
		{}

		/**
		 * Return structure for range based for loops over leaf nodes.
		 *
		 * @param listIdx tracking list id.
		 * @return
		 */
		const LeafsContainer<ThisType> leafs(const UINT& listIdx = 0) const
		{
			return LeafsContainer<ThisType>(this, listIdx);
		}
	};

}
#endif /* PARTITIONEDGRID_HPP_ */
