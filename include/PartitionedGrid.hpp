#ifndef PARTITIONEDGRID_HPP_
#define PARTITIONEDGRID_HPP_

#include <stdexcept>
#include "MappedGrid.hpp"

namespace felt
{
	template <
		typename T, UINT D, UINT P, class G=Grid<T, D>, UINT N=1
	>
	class PartitionedGrid : public G
	{
	public:
		typedef G								ChildGrid;
		typedef TrackedGrid<ChildGrid, D, N>	BranchGrid;
	protected:
		typedef typename ChildGrid::VecDu	VecDu;
		typedef typename ChildGrid::VecDi	VecDi;

		BranchGrid	m_grid_branch;

		static const VecDu	udims_child;
		static const VecDi	idims_child;

	public:
		virtual ~PartitionedGrid ()
		{}

		PartitionedGrid () : ChildGrid(), m_grid_branch()
		{}


		PartitionedGrid (
			const VecDu& dims, const VecDi& offset, const FLOAT& delta = 1
		) : ChildGrid(), m_grid_branch()
		{
			this->init(dims, offset, delta);
		}


		virtual void init (
			const VecDu& dims, const VecDi& offset = VecDi::Zero(),
			const FLOAT& delta = 1
		) {
			ChildGrid::init(dims, offset, delta);
		}


		const VecDi pos_child (const VecDi& pos) const
		{
			return (
				(pos - this->offset()).array() / idims_child.array()
			).matrix() + m_grid_branch.offset();
		}


		BranchGrid& branch ()
		{
			return m_grid_branch;
		}


		const BranchGrid& branch () const
		{
			return m_grid_branch;
		}
		

		ChildGrid& child (const VecDi& pos)
		{
			return m_grid_branch(pos);
		}


		const ChildGrid& child (const VecDi& pos) const
		{
			return m_grid_branch(pos);
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
			this->m_dims = dims_grid;

			VecDu dims_parent = (
				dims_grid.array() / udims_child.array()
			).matrix();

			if (
				(dims_parent.array() * udims_child.array()).matrix()
				!= dims_grid
			) {
				dims_parent += VecDu::Constant(1);
			}

			m_grid_branch.dims(dims_parent);

			for (UINT idx = 0; idx < m_grid_branch.data().size(); idx++)
			{
				ChildGrid& part = m_grid_branch.data()(idx);
				part.dims(udims_child);
			}
		}


		virtual void offset (const VecDi& offset_grid)
		{
			ChildGrid::offset(offset_grid);

			const VecDi& offset_parent = (
				offset_grid.array() / idims_child.array()
			).matrix();

			m_grid_branch.offset(offset_parent);

			for (UINT idx = 0; idx < m_grid_branch.data().size(); idx++)
			{
				ChildGrid& part = m_grid_branch.data()(idx);
				const VecDi& pos_part = m_grid_branch.index(idx);
				const VecDi& offset_part = (
					(
						(pos_part - m_grid_branch.offset()).array()
						* idims_child.array()
					).matrix() + offset_grid
				);

				part.offset(offset_part);
			}
		}


		const VecDi offset () const
		{
			return ChildGrid::offset();
		}


		/**
		 * Fill with a single value.
		 *
		 * @param val
		 */
		void fill (const T& val)
		{
			for (UINT idx = 0; idx < m_grid_branch.data().size(); idx++)
			{
				ChildGrid& child = m_grid_branch.data()(idx);
				child.fill(val);
			}
		}


		T& get (const VecDi& pos)
		{
			ChildGrid& leaf = m_grid_branch(pos_child(pos));
			return leaf.get(pos);
		}


		const T& get (const VecDi& pos) const
		{
			const ChildGrid& leaf = m_grid_branch(pos_child(pos));
			return leaf.get(pos);
		}

	};


	template <typename T, UINT D, UINT P, class G, UINT N>
	const typename PartitionedGrid<T, D, P, G, N>::VecDu
	PartitionedGrid<T, D, P, G, N>::udims_child
		= PartitionedGrid<T, D, P, G, N>::VecDu::Constant(P);
	template <typename T, UINT D, UINT P, class G, UINT N>
	const typename PartitionedGrid<T, D, P, G, N>::VecDi
	PartitionedGrid<T, D, P, G, N>::idims_child
		= PartitionedGrid<T, D, P, G, N>::VecDi::Constant(P);



	/**
	 * Specialisation for MappedGrid.
	 */

	template <typename T, UINT D, UINT P, UINT N>
	using MappedPartitionedGridBase = PartitionedGrid<
		T, D, P, MappedGrid<T, D, N>, N
	>;

	template <typename T, UINT D, UINT P, UINT N>
	class MappedPartitionedGrid : public MappedPartitionedGridBase<T, D, P, N>
	{
	protected:
		typedef MappedPartitionedGridBase<T, D, P, N>	Base;
	public:
		typedef typename Base::ChildGrid				ChildGrid;
	protected:
		typedef typename ChildGrid::VecDu				VecDu;
		typedef typename ChildGrid::VecDi				VecDi;
	public:
		typedef typename Base::BranchGrid				BranchGrid;


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
			branch.add(pos_child, arr_idx);
			branch(pos_child).add(pos, val, arr_idx);
		}


		void reset(const T& val, const UINT& arr_idx = 0)
		{
			BranchGrid& branch = this->branch();
			for (const VecDi& pos_child : branch.list(arr_idx))
				branch(pos_child).reset(val, arr_idx);

			branch.reset(arr_idx);
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
		typedef PosMappedPartitionedGridBase<D, P, N>	Base_t;
		typedef typename Base_t::ChildGrid				ChildGrid;
		typedef typename ChildGrid::VecDu				VecDu;
		typedef typename ChildGrid::VecDi				VecDi;
	public:
		typedef typename Base_t::BranchGrid				BranchGrid;

	protected:
		Lookup											m_grid_lookup;

	public:
		PosMappedPartitionedGrid () : Base_t(), m_grid_lookup()
		{}


		PosMappedPartitionedGrid (
			const VecDu& dims, const VecDi& offset, const FLOAT& delta = 1
		) : Base_t(), m_grid_lookup()
		{
			init(dims, offset, delta);
		}


		void init (
			const VecDu& dims, const VecDi& offset = VecDi::Zero(),
			const FLOAT& delta = 1
		) {
			Base_t::init(dims, offset, delta);
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
