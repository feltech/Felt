#ifndef PARTITIONEDGRID_HPP_
#define PARTITIONEDGRID_HPP_

#include <boost/iterator/iterator_facade.hpp>

#include "MappedGrid.hpp"

namespace felt
{
	template <UINT D, UINT P, class G, UINT N=1>
	class PartitionBase
	{
	public:
		typedef G							Child;
		typedef TrackedGrid<Child, D, N>	BranchGrid;
	protected:
		typedef typename BranchGrid::VecDu	VecDu;
		typedef typename BranchGrid::VecDi	VecDi;

		BranchGrid	m_grid_branch;

		static const VecDu	udims_child;
		static const VecDi	idims_child;

	public:
		virtual ~PartitionBase ()
		{}

		PartitionBase () : m_grid_branch()
		{}


		PartitionBase (
			const VecDu& dims, const VecDi& offset, const FLOAT& delta = 1
		) : m_grid_branch()
		{
			this->init(dims, offset);
		}


		void init (
			const VecDu& dims, const VecDi& offset = VecDi::Zero()
		) {
			this->dims(dims);
			this->offset(offset);
		}


		BranchGrid& branch ()
		{
			return m_grid_branch;
		}


		const BranchGrid& branch () const
		{
			return m_grid_branch;
		}
		

		Child& child (const VecDi& pos)
		{
			return m_grid_branch(pos);
		}


		const Child& child (const VecDi& pos) const
		{
			return m_grid_branch(pos);
		}


		/**
		 * Reshape grid.
		 *
		 * @param vec_NewDims
		 * @return
		 */
		virtual void dims (const VecDu& dims_grid)
		{
			VecDu dims_branch = (
				dims_grid.array() / udims_child.array()
			).matrix();

			if (
				(dims_branch.array() * udims_child.array()).matrix()
				!= dims_grid
			) {
				dims_branch += VecDu::Constant(1);
			}

			m_grid_branch.dims(dims_branch);
		}


		virtual void offset (const VecDi& offset_grid)
		{
			const VecDi& offset_parent = (
				offset_grid.array() / idims_child.array()
			).matrix();

			m_grid_branch.offset(offset_parent);
		}


		void add_child(const VecDi& pos, const UINT& arr_idx = 0)
		{
			this->branch().add(pos, arr_idx);
		}


		void remove_child(const VecDi& pos, const UINT& arr_idx = 0)
		{
			this->branch().remove(pos, arr_idx);
		}


		void reset(const UINT& arr_idx = 0)
		{
			this->branch().reset(arr_idx);
		}
	};


	template <UINT D, UINT P, class G, UINT N>
	const typename PartitionBase<D, P, G, N>::VecDu
	PartitionBase<D, P, G, N>::udims_child
		= PartitionBase<D, P, G, N>::VecDu::Constant(P);
	template <UINT D, UINT P, class G, UINT N>
	const typename PartitionBase<D, P, G, N>::VecDi
	PartitionBase<D, P, G, N>::idims_child
		= PartitionBase<D, P, G, N>::VecDi::Constant(P);



	template <typename T, UINT D, UINT P>
	class PartitionedArray : public PartitionBase<
		D, P, std::vector<T, Eigen::aligned_allocator<T> >, 1
	>
	{
	protected:
		typedef PartitionBase<
			D, P, std::vector<T, Eigen::aligned_allocator<T> >, 1
		> Base;

		typedef typename Base::VecDu	VecDu;
		typedef typename Base::VecDi	VecDi;

		VecDi	m_offset;

	public:
		PartitionedArray () : Base()
		{}

		PartitionedArray (
			const VecDu& dims, const VecDi& offset
		) : Base()
		{
			this->init(dims, offset);
		}

		void offset (const VecDi& offset_grid)
		{
			m_offset = offset_grid;
			Base::offset(offset_grid);
		}

		const VecDi pos_child (const VecDi& pos_leaf) const
		{
			return (
				(pos_leaf - m_offset).array() / Base::idims_child.array()
			).matrix() + this->branch().offset();
		}

		void add(const VecDi& pos, const T& val)
		{
			const VecDi& pos_child = this->pos_child(pos);
			this->child(pos_child).push_back(val);
			Base::add_child(pos_child);
		}

		void reset()
		{
			for (const VecDi& pos_child : this->branch().list())
				this->child(pos_child).clear();
			Base::reset();
		}
	};


	template <
		typename T, UINT D, UINT P, class G=Grid<T, D>, UINT N=1
	>
	class PartitionedGrid : public G, public PartitionBase<D, P, G, N>
	{
	public:
		typedef PartitionBase<D, P, G, N>		Base;
		typedef typename Base::Child			ChildGrid;
		typedef TrackedGrid<ChildGrid, D, N>	BranchGrid;
	protected:
		typedef typename ChildGrid::VecDu	VecDu;
		typedef typename ChildGrid::VecDi	VecDi;

	public:
		virtual ~PartitionedGrid ()
		{}

		PartitionedGrid () : Base(), ChildGrid()
		{}


		PartitionedGrid (
			const VecDu& dims, const VecDi& offset, const FLOAT& delta = 1
		) : Base(), ChildGrid()
		{
			this->init(dims, offset, delta);
		}


		void init (
			const VecDu& dims, const VecDi& offset = VecDi::Zero(),
			const FLOAT& delta = 1
		) {
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
				child.dims(Base::udims_child);
			}
		}


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
						* Base::idims_child.array()
					).matrix() + offset_grid
				);

				child.offset(offset_child);
			}
		}


		const VecDi offset () const
		{
			return ChildGrid::offset();
		}


		const VecDi pos_child (const VecDi& pos_leaf) const
		{
			return (
				(pos_leaf - offset()).array() / Base::idims_child.array()
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


		T& get (const VecDi& pos)
		{
			ChildGrid& leaf = this->branch()(pos_child(pos));
			return leaf.get(pos);
		}


		const T& get (const VecDi& pos) const
		{
			const ChildGrid& leaf = this->branch()(pos_child(pos));
			return leaf.get(pos);
		}

		void reset(const UINT& arr_idx = 0)
		{
			Base::reset(arr_idx);
		}
	};



	template <typename G>
	class LeafsContainer
	{
	private:
		typedef G							GridTree;
		typedef typename GridTree::VecDi	VecDi;
		typedef typename GridTree::PosArray	PosArray;

		const GridTree* m_pgrid;
		const UINT	m_listIdx;

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


		    bool equal(const iterator& other) const
		    {
		        return (
		        	m_it_child == other.m_it_child
					&& m_it_leaf == other.m_it_leaf
				);
		    }


		    const VecDi& dereference() const {
		    	const VecDi& pos = *m_it_leaf;
		    	return *m_it_leaf;
		    }
		};
	public:
		LeafsContainer(const GridTree* pgrid, const UINT& listIdx)
		: m_pgrid(pgrid), m_listIdx(listIdx)
		{}

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

		const iterator end() const
		{
			return iterator(
				m_pgrid, m_listIdx,
				m_pgrid->branch().list(m_listIdx).cend(),
				m_pgrid->branch().list(m_listIdx).cend()
			);
		}

		UINT size() const
		{
			UINT sum = 0;
			for (const VecDi& pos_child : m_pgrid->branch().list(m_listIdx))
				sum += m_pgrid->child(pos_child).list(m_listIdx).size();
			return sum;
		}
	};


	/**
	 * Specialisation for MappedGrid.
	 */
	template <typename T, UINT D, UINT P, UINT N>
	class MappedPartitionedGrid : public PartitionedGrid<
		T, D, P, MappedGrid<T, D, N>, N
	>
	{
	protected:
		typedef PartitionedGrid<
			T, D, P, MappedGrid<T, D, N>, N
		> Base;
		typedef MappedPartitionedGrid<T, D, P, N>		ThisType;
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
		MappedPartitionedGrid () : Base()
		{}


		MappedPartitionedGrid (
			const VecDu& dims, const VecDi& offset, const FLOAT& delta = 1
		) : Base(dims, offset, delta)
		{}


		void add(const VecDi& pos, const T& val, const UINT& arr_idx = 0)
		{
			BranchGrid& branch = this->branch();
			const VecDi& pos_child = this->pos_child(pos);
			branch(pos_child).add(pos, val, arr_idx);

			Base::add_child(pos_child, arr_idx);
		}

		const LeafsContainer<ThisType> leafs(const UINT& listIdx) const
		{
			return LeafsContainer<ThisType>(this, listIdx);
		}

		void reset(const T& val, const UINT& arr_idx = 0)
		{
			BranchGrid& branch = this->branch();
			for (const VecDi& pos_child : branch.list(arr_idx))
				branch(pos_child).reset(val, arr_idx);

			Base::reset(arr_idx);
		}
	};



	template <typename T, UINT D, UINT P, UINT N>
	class TrackedPartitionedGrid : public PartitionedGrid<
		T, D, P, TrackedGrid<T, D, N>, N
	>
	{
	protected:
		typedef PartitionedGrid<
			T, D, P, TrackedGrid<T, D, N>, N
		> Base;
		typedef TrackedPartitionedGrid<T, D, P, N>		ThisType;
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


		TrackedPartitionedGrid (
			const VecDu& dims, const VecDi& offset, const FLOAT& delta = 1
		) : Base(dims, offset, delta)
		{}


		void add(const VecDi& pos, const T& val, const UINT& arr_idx = 0)
		{
			BranchGrid& branch = this->branch();
			const VecDi& pos_child = this->pos_child(pos);
			branch(pos_child).add(pos, val, arr_idx);

			Base::add_child(pos_child, arr_idx);
		}


		const LeafsContainer<ThisType> leafs(const UINT& listIdx) const
		{
			return LeafsContainer<ThisType>(this, listIdx);
		}


		void reset(const T& val, const UINT& arr_idx = 0)
		{
			BranchGrid& branch = this->branch();
			for (const VecDi& pos_child : branch.list(arr_idx))
				branch(pos_child).reset(val, arr_idx);

			Base::reset(arr_idx);
		}


		void reset(const UINT& arr_idx = 0)
		{
			BranchGrid& branch = this->branch();
			for (const VecDi& pos_child : branch.list(arr_idx))
				branch(pos_child).reset(arr_idx);

			Base::reset(arr_idx);
		}


		void remove(const VecDi& pos, const UINT& arr_idx = 0)
		{
			const VecDi& pos_child = this->pos_child(pos);
			ChildGrid& child = this->child(pos_child);
			child.remove(pos, arr_idx);
			if (child.list(arr_idx).size() == 0)
				this->remove_child(pos_child, arr_idx);
		}
	};


	template <UINT D, UINT P, UINT N>
	using PosMappedPartitionedGridBase = PartitionedGrid<
		Eigen::Matrix<UINT, N, 1>, D, P, LookupGrid<D, N>, N
	>;

	template <UINT D, UINT P, UINT N>
	class PosMappedPartitionedGrid
	: public PosMappedPartitionedGridBase<D, P, N>
	{
	public:
		typedef LookupGrid<D, N>						Lookup;
		typedef typename Lookup::PosArray				PosArray;
	protected:
		typedef PosMappedPartitionedGridBase<D, P, N>	Base;
		typedef typename Base::ChildGrid				ChildGrid;
		typedef typename ChildGrid::VecDu				VecDu;
		typedef typename ChildGrid::VecDi				VecDi;
	public:
		typedef typename Base::BranchGrid				BranchGrid;

	protected:
		Lookup											m_grid_lookup;

	public:
		PosMappedPartitionedGrid () : Base(), m_grid_lookup()
		{}


		PosMappedPartitionedGrid (
			const VecDu& dims, const VecDi& offset, const FLOAT& delta = 1
		) : Base(), m_grid_lookup()
		{
			init(dims, offset, delta);
		}


		void init (
			const VecDu& dims, const VecDi& offset = VecDi::Zero(),
			const FLOAT& delta = 1
		) {
			Base::init(dims, offset, delta);
			lookup().init(dims, offset, delta);
		}


		Lookup& lookup()
		{
			return m_grid_lookup;
		}


		const Lookup& lookup() const
		{
			return m_grid_lookup;
		}


		void add(const VecDi& pos, const UINT& arr_idx = 0)
		{
			const VecDi& pos_child = this->pos_child(pos);
			BranchGrid& parent = this->branch();
			Lookup& lookup = this->lookup();
			ChildGrid& child = parent(pos_child);

			lookup.add(pos_child, arr_idx);
			child.add(pos, arr_idx);
		}


		void reset(const UINT& arr_idx = 0)
		{
			BranchGrid& parent = this->branch();
			Lookup& lookup = this->lookup();
			for (const VecDi& child_pos : lookup.list(arr_idx))
			{
				parent(child_pos).reset(arr_idx);
			}
			lookup.reset(arr_idx);
		}


		void remove(const VecDi& pos, const UINT& arr_idx = 0)
		{
			const VecDi& pos_child = this->pos_child(pos);
			BranchGrid& branch = this->branch();
			Lookup& lookup = this->lookup();
			ChildGrid& child = branch(pos_child);
			child.remove(pos, arr_idx);
			if (child.list(arr_idx).size() == 0)
				lookup.remove(pos_child, arr_idx);
		}
	};

}
#endif /* PARTITIONEDGRID_HPP_ */
